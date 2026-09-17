# Suficiência da parametrização da torção: quantos nós são necessários?
#
# Rodar de dentro da pasta avl/:   julia suficiencia_nos.jl
#
# A pergunta é se as cinco torções livres (nas seções da asa, com a raiz
# fixa como referência) bastam. A resposta é diferente para cada objetivo:
#
#   - Para o ARRASTO da asa isolada existe um piso analítico, CL²/(π AR),
#     que corresponde à carga elíptica. Se a parametrização alcança esse
#     piso, nenhum nó a mais pode ajudar: não há o que ganhar abaixo dele.
#
#   - Para a RESTRIÇÃO DE ESTOL o critério é pontual (o cl local de cada
#     faixa), e aí mais nós podem ajudar de verdade, porque dão controle
#     mais fino do carregamento local.
#
# O estudo abaixo mede as duas coisas variando o número de nós livres.

using Printf
using JSON
using Optim

const AVL   = "./avl.exe"
const BASE  = "aft.avl"
const DADOS = JSON.parsefile("dados_designtool.json")
const PP    = DADOS["ponto_de_projeto"]
const MACH, CL_PROJ = PP["M"], round(PP["CL"], digits = 4)
const SREF, BREF = PP["Sref"], DADOS["referencia"]["Bref"]
const SEMI = BREF/2
const ETAS = [0.0, 0.1011, 0.398, 0.56, 0.90, 1.0]
const TW_MIN, TW_MAX = -12.0, 6.0
const EPS_H = 3.0
const MACH_BAIXO, ALFAS_BASE = 0.2, (8.0, 14.0)
const ETA_SEGURA, FOLGA_ALPHA = 0.50, 0.5
const ETA_LIM, CLMAX_LIM = [0.1011, 0.398, 0.90], [1.774, 1.7985, 1.7338]

function roda_avl(cmds)
    tmp = "_sn.txt"; write(tmp, cmds)
    try; return read(pipeline(`$AVL`, stdin = tmp), String)
    finally; rm(tmp, force = true); end
end
function num(t, p)
    m = collect(eachmatch(p, t)); isempty(m) && return NaN
    parse(Float64, m[end].captures[1])
end

function escreve(destino, tw; so_asa = false)
    linhas = readlines(BASE); n = length(linhas)
    inicios = [i for i in 1:n if strip(linhas[i]) in ("SURFACE", "BODY")]
    i_asa = 0
    for i in inicios
        if strip(linhas[i]) == "SURFACE"
            j = i + 1
            while j <= n && (isempty(strip(linhas[j])) ||
                             startswith(strip(linhas[j]), "#")); j += 1; end
            j <= n && strip(linhas[j]) == "Wing" && (i_asa = i)
        end
    end
    prox = findfirst(>(i_asa), inicios)
    fim = prox === nothing ? n : inicios[prox] - 1
    faixa = so_asa ? (1:fim) : (1:n)
    saida = String[]; i_sec = 0; espera = false
    for i in faixa
        ln = linhas[i]; t = strip(ln)
        dentro = (i >= i_asa && i <= fim)
        if dentro && t == "SECTION"; espera = true
        elseif dentro && espera && !isempty(t) && !startswith(t, "#")
            i_sec += 1; p = split(t); p[5] = @sprintf("%.4f", tw[i_sec])
            push!(saida, join(p, "  ")); espera = false; continue
        end
        push!(saida, ln)
    end
    write(destino, join(saida, "\n") * "\n")
end

function avalia(tw; so_asa = false, trim = false, faixas = false)
    escreve("_sn.avl", tw; so_asa = so_asa)
    c = ["load _sn.avl", "oper", "m", "mn $MACH", ""]
    trim && push!(c, "d2 pm 0")
    append!(c, ["a c $CL_PROJ", "x"])
    faixas && append!(c, ["fs", ""])
    append!(c, ["", "quit"])
    s = roda_avl(join(c, "\n") * "\n"); rm("_sn.avl", force = true)
    (CDff = num(s, r"CDff\s*=\s*([-\d.Ee+]+)"),
     CL = num(s, r"CLtot\s*=\s*([-\d.]+)"),
     e = num(s, r"\se =\s+([-\d.]+)"), saida = s)
end

function faixas_asa(saida)
    tr = split(saida, r"Surface # 1\s+Wing")[2]
    bl = String(split(tr, "Surface # 2")[1])
    L = [[parse(Float64, v) for v in split(strip(m.match))[2:end]]
         for m in eachmatch(r"^\s*\d+(?:\s+[-\dEe.+]+){12}\s*$"m, bl)]
    d = reduce(hcat, L)'
    (eta = d[:, 1]./SEMI, ccl = d[:, 4], cl_norm = d[:, 6])
end

# --------------------------------------------------------------------
# PARAMETRIZAÇÃO COM k NÓS LIVRES
#
# Os nós livres ficam em estações escolhidas da asa; as demais estações
# recebem torção por interpolação linear entre eles. Com k = 5 recai-se
# no caso de todas as seções livres.

const CONJUNTOS = Dict(
    1 => [1.0],
    2 => [0.398, 1.0],
    3 => [0.398, 0.90, 1.0],
    4 => [0.1011, 0.398, 0.90, 1.0],
    5 => [0.1011, 0.398, 0.56, 0.90, 1.0])

"Matriz que leva os k nós livres às torções das seis seções."
function mapa(k)
    nos = vcat(0.0, CONJUNTOS[k])          # a raiz entra fixa em zero
    M = zeros(length(ETAS), k)
    for (i, η) in enumerate(ETAS)
        if η <= nos[1]
            continue                        # raiz: torção nula
        elseif η >= nos[end]
            M[i, k] = 1.0
        else
            j = findlast(<=(η), nos)
            w = (η - nos[j])/(nos[j+1] - nos[j])
            j > 1 && (M[i, j-1] = 1 - w)
            M[i, j] = w
        end
    end
    return M
end

function modelo(M; so_asa, trim)
    k = size(M, 2)
    f0 = avalia(M*zeros(k); so_asa = so_asa, trim = trim).CDff
    fp, fm = zeros(k), zeros(k)
    for j in 1:k
        e = zeros(k); e[j] = EPS_H
        fp[j] = avalia(M*e; so_asa = so_asa, trim = trim).CDff
        fm[j] = avalia(M*(-e); so_asa = so_asa, trim = trim).CDff
    end
    g = (fp .- fm)./(2*EPS_H); H = zeros(k, k)
    for j in 1:k; H[j, j] = (fp[j] + fm[j] - 2*f0)/EPS_H^2; end
    for j in 1:k-1, l in j+1:k
        e = zeros(k); e[j] = e[l] = EPS_H
        H[j, l] = H[l, j] = (avalia(M*e; so_asa = so_asa, trim = trim).CDff -
                             fp[j] - fp[l] + f0)/EPS_H^2
    end
    (f0 = f0, g = g, H = H)
end
cd_mod(m, x) = m.f0 + m.g'x + 0.5*x'*m.H*x

function otimiza(m, k; pen = nothing, x0 = zeros(k))
    x = copy(x0); n_it = 0
    for w in (pen === nothing ? [0.0] : [1e2, 1e3, 1e4, 1e5, 1e6])
        f(z) = 1e4*cd_mod(m, z) + (pen === nothing ? 0.0 : w*pen(z))
        r = optimize(f, fill(TW_MIN, k), fill(TW_MAX, k), x, Fminbox(LBFGS()),
                     Optim.Options(iterations = 300, g_tol = 1e-11,
                                   store_trace = true))
        x = Optim.minimizer(r); n_it += length(Optim.trace(r))
    end
    return x, n_it
end

# --------------------------------------------------------------------
# ESTOL (para o estudo da etapa C)

function modelo_estol(M)
    k = size(M, 2)
    function amostra(x, alfa)
        escreve("_sn.avl", M*x)
        s = roda_avl(join(["load _sn.avl", "oper", "m", "mn $MACH_BAIXO", "",
                           "d2 d2 0", "a a $alfa", "x", "fs", "", "",
                           "quit"], "\n") * "\n")
        rm("_sn.avl", force = true); faixas_asa(s)
    end
    d1 = amostra(zeros(k), ALFAS_BASE[1]); d2 = amostra(zeros(k), ALFAS_BASE[2])
    b = (d2.cl_norm .- d1.cl_norm)./(ALFAS_BASE[2] - ALFAS_BASE[1])
    a = d1.cl_norm .- b.*ALFAS_BASE[1]
    C = zeros(k, length(d1.eta))
    for j in 1:k
        e = zeros(k); e[j] = 1.0
        C[j, :] = amostra(e, ALFAS_BASE[1]).cl_norm .- d1.cl_norm
    end
    pesos(η) = η <= ETA_LIM[1] ? (1.0, 0.0, 0.0) :
               η >= ETA_LIM[3] ? (0.0, 0.0, 1.0) :
               η <= ETA_LIM[2] ?
               (1 - (η-ETA_LIM[1])/(ETA_LIM[2]-ETA_LIM[1]),
                (η-ETA_LIM[1])/(ETA_LIM[2]-ETA_LIM[1]), 0.0) :
               (0.0, 1 - (η-ETA_LIM[2])/(ETA_LIM[3]-ETA_LIM[2]),
                (η-ETA_LIM[2])/(ETA_LIM[3]-ETA_LIM[2]))
    lim = [sum(CLMAX_LIM[i]*w for (i, w) in enumerate(pesos(η))) for η in d1.eta]
    (eta = d1.eta, a = a, b = b, C = C, lim = lim)
end
alfas_estol(me, x) = (me.lim .- me.a .- me.C'x)./me.b
function folga(me, x)
    al = alfas_estol(me, x); int = me.eta .<= ETA_SEGURA
    minimum(al[.!int]) - minimum(al[int])
end

# ====================================================================
println("="^76)
println("SUFICIÊNCIA DA PARAMETRIZAÇÃO DA TORÇÃO")
println("="^76)
piso = CL_PROJ^2/(pi*BREF^2/SREF)
@printf("piso analítico do arrasto da asa isolada, CL²/(π AR) = %.6f\n", piso)
println("(é a carga elíptica; nenhuma parametrização pode ficar abaixo dele)\n")

println("ASA ISOLADA -- o arrasto contra o número de nós livres:")
@printf("  %-6s %10s %10s %12s %9s %10s\n",
        "nós", "CDff", "Oswald", "acima do piso", "iterações", "desvio elíp.")
for k in 1:5
    M = mapa(k)
    m = modelo(M; so_asa = true, trim = false)
    x, nit = otimiza(m, k)
    av = avalia(M*x; so_asa = true, faixas = true)
    fx = faixas_asa(av.saida)
    elp = sqrt.(max.(0.0, 1 .- fx.eta.^2))
    dsv = maximum(abs.(fx.ccl./(4*SREF*av.CL/(pi*BREF)) .- elp))
    @printf("  %-6d %10.6f %10.4f %9.2f counts %9d %10.3f\n",
            k, av.CDff, av.e, 1e4*(av.CDff - piso), nit, dsv)
end

println("\nAERONAVE COMPLETA COM RESTRIÇÃO DE ESTOL -- aqui o critério é")
println("pontual, e mais nós podem realmente ajudar:")
@printf("  %-6s %10s %10s %12s %9s\n",
        "nós", "CDff", "Oswald", "folga estol", "iterações")
for k in 2:5
    M = mapa(k)
    m = modelo(M; so_asa = false, trim = true)
    me = modelo_estol(M)
    int = me.eta .<= ETA_SEGURA
    a0 = alfas_estol(me, zeros(k))
    i_ref = argmin([int[j] ? a0[j] : Inf for j in eachindex(a0)])
    ext = findall(.!int)
    d = [(me.lim[j]-me.a[j])/me.b[j] - (me.lim[i_ref]-me.a[i_ref])/me.b[i_ref]
         for j in ext]
    W = hcat([-me.C[:, j]./me.b[j] .+ me.C[:, i_ref]./me.b[i_ref] for j in ext]...)
    pen(z) = sum(min.(d .+ vec(W'z) .- FOLGA_ALPHA, 0.0).^2)
    melhor = nothing
    for x0 in (zeros(k), fill(-2.0, k), [j <= k÷2 ? 2.0 : -5.0 for j in 1:k])
        x, nit = otimiza(m, k; pen = pen, x0 = x0)
        if folga(me, x) >= FOLGA_ALPHA - 0.05 &&
           (melhor === nothing || cd_mod(m, x) < melhor[2])
            melhor = (x, cd_mod(m, x), nit)
        end
    end
    if melhor === nothing
        @printf("  %-6d %10s\n", k, "inviável"); continue
    end
    av = avalia(M*melhor[1]; trim = true)
    @printf("  %-6d %10.6f %10.4f %+11.2f° %9d\n",
            k, av.CDff, av.e, folga(me, melhor[1]), melhor[3])
end
