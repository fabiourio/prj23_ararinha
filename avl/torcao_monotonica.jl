# A torção precisa mesmo ser não monotônica?
#
# Rodar de dentro da pasta avl/:   julia torcao_monotonica.jl
#
# A otimização livre da etapa 3 devolveu uma torção serrilhada e com
# barriga interna. Asa real não é assim: a prática é wash-in máximo na
# raiz caindo de forma monótona até a ponta. A pergunta é se essa forma
# clássica atende a restrição de estol e a que custo de arrasto.
#
# O teste compara, no mesmo problema e com a mesma margem:
#   LIVRE       torções sem vínculo entre estações
#   MONOTÔNICA  torção que só decresce da raiz para a ponta
#   LINEAR      washout linear em η, um único grau de liberdade
#
# A raiz é a referência de gauge, então "wash-in na raiz" aparece aqui
# como torções negativas crescendo em módulo para fora.
#
# O problema é resolvido na forma exata em que ele nasce: a CL fixo o
# arrasto induzido é quadrático na torção e a margem de estol é afim,
# de modo que cada caso é um QP de restrições lineares. Ele é resolvido
# por pontos interiores (IPNewton), sem penalidade. Isso importa porque
# penalidade mal escalada é uma das explicações candidatas para o
# serrilhado, e resolvendo exato essa explicação sai de cena.

using Printf
using JSON
using Optim
using LinearAlgebra
using Plots

gr()

const AVL, BASE, SAIDA = "./avl.exe", "aft.avl", "resultados"
const DADOS = JSON.parsefile("dados_designtool.json")
const PP = DADOS["ponto_de_projeto"]
const MACH, CL_PROJ = PP["M"], round(PP["CL"], digits = 4)
const SREF, BREF = PP["Sref"], DADOS["referencia"]["Bref"]
const SEMI = BREF/2
const TW_MIN, TW_MAX = -12.0, 6.0
const EPS_H = 3.0
const MACH_BAIXO, ALFAS_BASE = 0.2, (8.0, 14.0)
const ETA_AILERON, MARGEM = 0.56, 1.0
const ETA_LIM, CLMAX_LIM = [0.1011, 0.398, 0.90], [1.774, 1.7985, 1.7338]
const ESCALA = 1e4
const PAL = ["#2a78d6", "#eb6834", "#1baf7a"]
const INK, INK2 = "#0b0b0b", "#52514e"
const ESTILO = (framestyle = :axes, gridcolor = "#e1e0d9", gridalpha = 1.0,
                gridlinewidth = 0.7, foreground_color_axis = INK2,
                foreground_color_border = "#c3c2b7",
                foreground_color_text = INK2, tickfontsize = 9,
                guidefontsize = 10, legendfontsize = 8,
                background_color = :white)

mkpath(SAIDA)

roda(c) = (write("_tm.txt", c);
           try read(pipeline(`$AVL`, stdin = "_tm.txt"), String)
           finally rm("_tm.txt", force = true) end)
function num(t, p)
    m = collect(eachmatch(p, t)); isempty(m) && return NaN
    parse(Float64, m[end].captures[1])
end

function limites_asa()
    L = readlines(BASE); n = length(L)
    ini = [i for i in 1:n if strip(L[i]) in ("SURFACE", "BODY")]
    ia = 0
    for i in ini
        strip(L[i]) == "SURFACE" || continue
        j = i + 1
        while j <= n && (isempty(strip(L[j])) || startswith(strip(L[j]), "#"))
            j += 1
        end
        j <= n && strip(L[j]) == "Wing" && (ia = i)
    end
    prox = findfirst(>(ia), ini)
    return L, ia, (prox === nothing ? n : ini[prox] - 1)
end

const LINHAS, I_ASA, FIM_ASA = limites_asa()

function le_etas()
    ys = Float64[]
    espera = false
    for i in I_ASA:FIM_ASA
        t = strip(split(LINHAS[i], "#")[1])
        if strip(LINHAS[i]) == "SECTION"
            espera = true
        elseif espera && !isempty(t)
            push!(ys, parse(Float64, split(t)[2]))
            espera = false
        end
    end
    return ys ./ ys[end]
end

const ETAS = le_etas()
const NV = length(ETAS) - 1
torcoes(x) = vcat(0.0, x)

function escreve(dest, tw)
    saida = String[]; i_sec, esp = 0, false
    for i in eachindex(LINHAS)
        ln = LINHAS[i]; t = strip(ln)
        dentro = I_ASA <= i <= FIM_ASA
        if dentro && t == "SECTION"; esp = true
        elseif dentro && esp && !isempty(t) && !startswith(t, "#")
            i_sec += 1; p = split(t); p[5] = @sprintf("%.4f", tw[i_sec])
            push!(saida, join(p, "  ")); esp = false; continue
        end
        push!(saida, ln)
    end
    write(dest, join(saida, "\n") * "\n")
end

function avalia(tw; faixas = false)
    escreve("_tm.avl", tw)
    c = ["load _tm.avl", "oper", "m", "mn $MACH", "", "d2 pm 0",
         "a c $CL_PROJ", "x"]
    faixas && append!(c, ["fs", ""])
    append!(c, ["", "quit"])
    s = roda(join(c, "\n") * "\n"); rm("_tm.avl", force = true)
    (CDff = num(s, r"CDff\s*=\s*([-\d.Ee+]+)"),
     CDvis = num(s, r"CDvis\s*=\s*([-\d.Ee+]+)"),
     CL = num(s, r"CLtot\s*=\s*([-\d.]+)"),
     e = num(s, r"\se =\s+([-\d.]+)"),
     de = num(s, r"elevator\s*=\s*([-\d.]+)"), saida = s)
end

function faixas_asa(s)
    bl = String(split(split(s, r"Surface # 1\s+Wing")[2], "Surface # 2")[1])
    L = [[parse(Float64, v) for v in split(strip(m.match))[2:end]]
         for m in eachmatch(r"^\s*\d+(?:\s+[-\dEe.+]+){12}\s*$"m, bl)]
    d = reduce(hcat, L)'
    (eta = d[:, 1]./SEMI, ccl = d[:, 4], cl_norm = d[:, 6])
end

function amostra(x, alfa)
    escreve("_tm.avl", torcoes(x))
    s = roda(join(["load _tm.avl", "oper", "m", "mn $MACH_BAIXO", "",
                   "d2 d2 0", "a a $alfa", "x", "fs", "", "", "quit"], "\n")*"\n")
    rm("_tm.avl", force = true); faixas_asa(s)
end

# --------------------------------------------------------------------
# MODELO
#
# A CL fixo a distribuição de sustentação é afim na torção, logo o
# arrasto induzido é exatamente quadrático e o cl de cada faixa é
# exatamente afim. Os dois modelos abaixo não são ajustes, são a forma
# fechada obtida de diferenças finitas. Construí-los custa cerca de 80
# rodadas de AVL, por isso ficam em cache.

const CACHE = joinpath(SAIDA, "modelo_torcao.json")

function constroi_modelo()
    println("modelo de arrasto: ", 1 + 2NV + NV*(NV-1)÷2, " rodadas de AVL...")
    f0 = avalia(torcoes(zeros(NV))).CDff
    fp, fm = zeros(NV), zeros(NV)
    for j in 1:NV
        d = zeros(NV); d[j] = EPS_H
        fp[j] = avalia(torcoes(d)).CDff
        fm[j] = avalia(torcoes(-d)).CDff
    end
    G = (fp .- fm)./(2*EPS_H)
    H = zeros(NV, NV)
    for j in 1:NV
        H[j, j] = (fp[j] + fm[j] - 2*f0)/EPS_H^2
    end
    for j in 1:NV-1, k in j+1:NV
        d = zeros(NV); d[j] = d[k] = EPS_H
        H[j, k] = H[k, j] =
            (avalia(torcoes(d)).CDff - fp[j] - fp[k] + f0)/EPS_H^2
    end

    println("modelo de estol: ", 2 + NV, " rodadas de AVL...")
    s1 = amostra(zeros(NV), ALFAS_BASE[1])
    s2 = amostra(zeros(NV), ALFAS_BASE[2])
    bb = (s2.cl_norm .- s1.cl_norm)./(ALFAS_BASE[2] - ALFAS_BASE[1])
    aa = s1.cl_norm .- bb.*ALFAS_BASE[1]
    CC = zeros(NV, length(s1.eta))
    for j in 1:NV
        d = zeros(NV); d[j] = 1.0
        CC[j, :] = amostra(d, ALFAS_BASE[1]).cl_norm .- s1.cl_norm
    end
    return (f0 = f0, G = G, H = H, aa = aa, bb = bb, CC = CC, eta = s1.eta)
end

function carrega_modelo()
    if isfile(CACHE) && mtime(CACHE) > mtime(BASE)
        d = JSON.parsefile(CACHE)
        if d["nv"] == NV
            println("modelo lido do cache ($CACHE)")
            lin(v) = reduce(hcat, [Float64.(r) for r in v])'
            return (f0 = d["f0"], G = Float64.(d["G"]), H = lin(d["H"]),
                    aa = Float64.(d["aa"]), bb = Float64.(d["bb"]),
                    CC = lin(d["CC"]), eta = Float64.(d["eta"]))
        end
    end
    m = constroi_modelo()
    open(CACHE, "w") do io
        JSON.print(io, Dict("nv" => NV, "f0" => m.f0, "G" => m.G,
                            "H" => [m.H[i, :] for i in 1:size(m.H, 1)],
                            "aa" => m.aa, "bb" => m.bb,
                            "CC" => [m.CC[i, :] for i in 1:size(m.CC, 1)],
                            "eta" => m.eta))
    end
    return m
end

const M = carrega_modelo()
const F0, G, HH = M.f0, M.G, M.H
cd_mod(x) = F0 + G'x + 0.5*x'*HH*x

# Condicionamento. Se a Hessiana for definida positiva o QP é convexo,
# cada caso tem um único mínimo global e não existe ambiguidade de
# mínimo local para explicar a forma que sair.
const AUTOV = sort(eigvals(Symmetric(HH)); rev = true)
@printf("\nautovalores da Hessiana do arrasto\n")
@printf("  maior %.4e   menor %.4e   razão %.1f\n",
        AUTOV[1], AUTOV[end], AUTOV[1]/AUTOV[end])
println(all(AUTOV .> 0) ?
    "  definida positiva: o problema é convexo e o ótimo é único" :
    "  NÃO é definida positiva: há direção de curvatura nula ou negativa")

# modelo de estol, afim
function peso(η)
    η <= ETA_LIM[1] && return (1.0, 0.0, 0.0)
    η >= ETA_LIM[3] && return (0.0, 0.0, 1.0)
    η <= ETA_LIM[2] ?
        (w = (η-ETA_LIM[1])/(ETA_LIM[2]-ETA_LIM[1]); (1-w, w, 0.0)) :
        (w = (η-ETA_LIM[2])/(ETA_LIM[3]-ETA_LIM[2]); (0.0, 1-w, w))
end
const ETA_F = M.eta
const LIM = [sum(CLMAX_LIM[i]*w for (i, w) in enumerate(peso(η))) for η in ETA_F]
alfa_estol(x) = (LIM .- M.aa .- M.CC'x)./M.bb
function margem_ail(x)
    al = alfa_estol(x); dentro = ETA_F .< ETA_AILERON
    minimum(al[.!dentro]) - minimum(al[dentro])
end
eta_crit(x) = ETA_F[argmin(alfa_estol(x))]

# --------------------------------------------------------------------
# RESTRIÇÕES, TODAS LINEARES
#
# Estol: fixada a faixa interna i_ref que estola primeiro, exigir
# α_estol(j) >= α_estol(i_ref) + MARGEM para toda faixa j do aileron.
# Varrer i_ref sobre as faixas internas e ficar com o melhor resultado
# é exato, porque a solução viável verdadeira satisfaz essa forma para
# o i_ref que de fato minimiza.

function restricao_estol(i_ref, margem = MARGEM)
    ext = findall(ETA_F .>= ETA_AILERON)
    d = [(LIM[j]-M.aa[j])/M.bb[j] - (LIM[i_ref]-M.aa[i_ref])/M.bb[i_ref]
         for j in ext]
    W = hcat([-M.CC[:, j]./M.bb[j] .+ M.CC[:, i_ref]./M.bb[i_ref] for j in ext]...)
    return Matrix(W'), margem .- d          # A*x >= lc
end

# Taxa de variação da torção entre estações vizinhas, com tw = [0; x]:
#   lo <= tw[i+1] - tw[i] <= hi
# Monotonicidade é o caso (-Inf, 0), que só proíbe subir. Suavidade é o
# caso (-Δ, +Δ), que limita o salto nos dois sentidos e é a restrição de
# fato física: o VLM só vale para variação suave de incidência ao longo
# da envergadura, e a asa ainda precisa ser fabricável.
function restricao_taxa(lo, hi)
    A = zeros(NV, NV)
    A[1, 1] = 1.0
    for i in 2:NV
        A[i, i] = 1.0; A[i, i-1] = -1.0
    end
    return A, fill(lo, NV), fill(hi, NV)
end

# Fase 1. O IPNewton exige partir de um ponto estritamente interior, e um
# chute qualquer de torção não satisfaz as restrições de estol. Aqui o
# ponto de partida é obtido minimizando a violação com uma folga δ, de
# modo que o resultado fica dentro da região viável e não sobre a borda.
function ponto_interior(A, lc, uc, lo, hi, x0; δ = 0.05)
    ml, mu = isfinite.(lc), isfinite.(uc)
    lcf = [ml[i] ? lc[i] + δ : -Inf for i in eachindex(lc)]
    ucf = [mu[i] ? uc[i] - δ :  Inf for i in eachindex(uc)]
    function residuos(x)
        c = A*x
        a = [ml[i] ? max(lcf[i] - c[i], 0.0) : 0.0 for i in eachindex(c)]
        b = [mu[i] ? max(c[i] - ucf[i], 0.0) : 0.0 for i in eachindex(c)]
        return a, b
    end
    function viol(x)
        a, b = residuos(x)
        return sum(a.^2) + sum(b.^2)
    end
    function viol_g!(g, x)
        a, b = residuos(x)
        g .= 2 .* (A' * (b .- a))
        return g
    end
    r = optimize(viol, viol_g!, lo, hi, clamp.(x0, lo, hi), Fminbox(LBFGS()),
                 Optim.Options(iterations = 600, g_tol = 1e-14))
    x = clamp.(Optim.minimizer(r), lo, hi)
    return viol(x) < 1e-9 ? x : nothing
end

function resolve_qp(Hq, Gq, A, lc, uc, lo, hi, x0)
    fun(z) = ESCALA*(F0 + Gq'z + 0.5*z'*Hq*z)
    fun_grad!(g, z) = (g .= ESCALA .* (Gq .+ Hq*z))
    fun_hess!(h, z) = (h .= ESCALA .* Hq)
    con_c!(c, z) = (c .= A*z; c)
    con_j!(J, z) = (J .= A; J)
    con_h!(h, z, λ) = h          # restrições lineares não entram na Hessiana
    df = TwiceDifferentiable(fun, fun_grad!, fun_hess!, x0)
    dfc = TwiceDifferentiableConstraints(con_c!, con_j!, con_h!, lo, hi, lc, uc)
    r = optimize(df, dfc, copy(x0), IPNewton(),
                 Optim.Options(iterations = 400))
    return Optim.minimizer(r)
end

rampa(a, b) = [a + (b-a)*(ETAS[i+1]-ETAS[2])/(1-ETAS[2]) for i in 1:NV]
const PARTIDAS = [rampa(-0.2, -4.0), rampa(-0.5, -9.0), rampa(-1.0, -11.0),
                  rampa(-2.0, -6.0), zeros(NV) .- 1.0, zeros(NV) .- 5.0]

function resolve(; taxa = nothing, margem = MARGEM)
    melhor = nothing
    for i_ref in findall(ETA_F .< ETA_AILERON)
        Ae, lce = restricao_estol(i_ref, margem)
        A, lc, uc = Ae, lce, fill(Inf, length(lce))
        if taxa !== nothing
            At, lct, uct = restricao_taxa(taxa[1], taxa[2])
            A = vcat(A, At); lc = vcat(lc, lct); uc = vcat(uc, uct)
        end
        lo, hi = fill(TW_MIN + 1e-3, NV), fill(TW_MAX - 1e-3, NV)
        for x0 in PARTIDAS
            xi = ponto_interior(A, lc, uc, lo, hi, x0)
            xi === nothing && continue
            x = try resolve_qp(HH, G, A, lc, uc, fill(TW_MIN, NV),
                               fill(TW_MAX, NV), xi) catch; xi end
            any(isnan, x) && continue
            dtw = diff(torcoes(x))
            viavel = eta_crit(x) < ETA_AILERON &&
                     margem_ail(x) >= margem - 0.05 &&
                     all(x .>= TW_MIN - 1e-6) && all(x .<= TW_MAX + 1e-6) &&
                     (taxa === nothing ||
                      (minimum(dtw) >= taxa[1] - 1e-4 &&
                       maximum(dtw) <= taxa[2] + 1e-4))
            viavel && (melhor === nothing || cd_mod(x) < melhor[2]) &&
                (melhor = (x, cd_mod(x)))
        end
    end
    return melhor
end

# Casos. Além dos dois extremos entra a família suave, que é a resposta
# de engenharia: não proíbe wash-in interno, só proíbe o degrau.
const CASOS = [("LIVRE", nothing), ("MONOTÔNICA", (-Inf, 0.0)),
               ("SUAVE 2°/est", (-2.0, 2.0)), ("SUAVE 3°/est", (-3.0, 3.0)),
               ("SUAVE 4°/est", (-4.0, 4.0)),
               ("SUAVE+MONO 3°", (-3.0, 0.0))]

println("\nresolvendo...")
res = Dict{String,Any}()
for (nome, taxa) in CASOS
    s = resolve(taxa = taxa)
    if s === nothing
        println("  $nome: inviável"); continue
    end
    res[nome] = (x = s[1], av = avalia(torcoes(s[1]); faixas = true))
    @printf("  %-14s ok\n", nome)
end

# washout linear: um grau de liberdade só
function busca_linear()
    melhor = nothing
    for t in range(-12.0, 0.0; length = 241)
        x = [t*ETAS[i+1] for i in 1:NV]      # linear em η, zero na raiz
        if eta_crit(x) < ETA_AILERON && margem_ail(x) >= MARGEM - 0.05
            (melhor === nothing || cd_mod(x) < melhor[2]) &&
                (melhor = (x, cd_mod(x)))
        end
    end
    return melhor
end
melhor_lin = busca_linear()
if melhor_lin === nothing
    println("  LINEAR: inviável, nenhum washout linear entre 0 e -12° ",
            "tira o estol de dentro do aileron com a margem exigida")
else
    res["LINEAR"] = (x = melhor_lin[1],
                     av = avalia(torcoes(melhor_lin[1]); faixas = true))
end

# --------------------------------------------------------------------
# RESULTADO

# quanto o modelo erra em relação ao AVL, para saber se a diferença
# entre as parametrizações está acima do erro do próprio modelo
for nome in keys(res)
    r = res[nome]
    @printf("modelo x AVL em %-11s previsto %.6f  AVL %.6f  erro %+.2f counts\n",
            nome, cd_mod(r.x), r.av.CDff, 1e4*(cd_mod(r.x) - r.av.CDff))
end

const ORDEM = vcat([n for (n, _) in CASOS], "LINEAR")
salto(x) = maximum(abs.(diff(torcoes(x))))

println("\n", "="^78)
println("QUE VÍNCULO A TORÇÃO PRECISA OBEDECER? (margem de $(MARGEM)° de α)")
println("="^78)
@printf("%-14s %10s %8s %8s %8s %9s %8s\n", "caso", "CDff", "counts",
        "Oswald", "estol η", "profundor", "salto")
base = avalia(torcoes(zeros(NV)))
@printf("%-14s %10.6f %8s %8.4f %8.3f %8.2f° %7.1f°\n", "sem torção",
        base.CDff, "-", base.e, eta_crit(zeros(NV)), base.de, 0.0)
for nome in ORDEM
    haskey(res, nome) || continue
    r = res[nome]
    @printf("%-14s %10.6f %+8.2f %8.4f %8.3f %8.2f° %7.1f°\n", nome, r.av.CDff,
            1e4*(r.av.CDff - base.CDff), r.av.e, eta_crit(r.x), r.av.de,
            salto(r.x))
end
println("\n(salto = maior diferença de torção entre estações vizinhas; é o que")
println(" mede se a solução é fabricável e se o VLM ainda vale para ela)")

println("\ntorções [graus], da raiz para a ponta")
for nome in ORDEM
    haskey(res, nome) || continue
    @printf("  %-14s %s\n", nome,
            join([@sprintf("%6.2f", v) for v in torcoes(res[nome].x)], ""))
end

if haskey(res, "LIVRE") && haskey(res, "MONOTÔNICA")
    d = 1e4*(res["MONOTÔNICA"].av.CDff - res["LIVRE"].av.CDff)
    @printf("\n  a monotonicidade custa %+.2f counts sobre a solução livre\n", d)
end
for nome in ORDEM
    (nome == "LIVRE" || !haskey(res, nome)) && continue
    @printf("  %-14s custa %+7.2f counts sobre a livre e %+7.2f sobre a asa sem torção\n",
            nome, 1e4*(res[nome].av.CDff - res["LIVRE"].av.CDff),
            1e4*(res[nome].av.CDff - base.CDff))
end

open(joinpath(SAIDA, "torcao_monotonica.txt"), "w") do io
    println(io, "vínculos da torção, margem $(MARGEM)° de α, restrição de estol da FAR 25.203")
    @printf(io, "Hessiana do arrasto: maior %.4e, menor %.4e, razão %.1f (definida positiva)\n",
            AUTOV[1], AUTOV[end], AUTOV[1]/AUTOV[end])
    @printf(io, "%-14s CDff %.6f  counts %8s  e %.4f  estol η %.3f\n",
            "sem torção", base.CDff, "-", base.e, eta_crit(zeros(NV)))
    for nome in ORDEM
        haskey(res, nome) || continue
        r = res[nome]
        @printf(io, "%-14s CDff %.6f  counts %+7.2f  e %.4f  estol η %.3f  salto %.1f  tw %s\n",
                nome, r.av.CDff, 1e4*(r.av.CDff - base.CDff), r.av.e,
                eta_crit(r.x), salto(r.x),
                join([@sprintf("%.2f", v) for v in torcoes(r.x)], " "))
    end
end

MOSTRAR = ["LIVRE", "MONOTÔNICA", "SUAVE 3°/est", "SUAVE+MONO 3°"]
CORES = ["#2a78d6", "#eb6834", "#1baf7a", "#8a56c4"]
p = plot(; xlabel = "η = 2y/b", ylabel = "torção [graus]",
         title = "torção que atende a restrição de estol", titlefontsize = 11,
         titlelocation = :left, legend = :bottomleft, xlims = (0, 1), ESTILO...)
vspan!(p, [ETA_AILERON, 0.90]; color = "#dcd9cd", alpha = 0.45,
       linewidth = 0, label = "aileron")
hline!(p, [0.0]; color = "#c3c2b7", linewidth = 0.8, label = "")
for (i, nome) in enumerate(MOSTRAR)
    haskey(res, nome) || continue
    plot!(p, ETAS, torcoes(res[nome].x); color = CORES[i], linewidth = 2.4,
          marker = :circle, markersize = 5, markerstrokecolor = CORES[i],
          label = @sprintf("%s  (%+.1f counts)", nome,
                           1e4*(res[nome].av.CDff - base.CDff)))
end
savefig(plot(p; size = (860, 520), dpi = 200),
        joinpath(SAIDA, "torcao_monotonica.png"))
println("\nfigura: $SAIDA/torcao_monotonica.png")

# --------------------------------------------------------------------
# O QUE A CONFORMIDADE CUSTA
#
# A asa sem torção estola em η = 0.842, dentro do aileron, logo ela não
# atende a FAR 25.203 e alguma torção é obrigatória. A pergunta de
# projeto não é se a torção compensa em arrasto, é qual o menor preço
# da conformidade. O preço depende de duas escolhas: a margem exigida e
# o quanto se admite de salto de incidência entre estações.

const DELTAS = [2.0, 3.0, 4.0, 6.0]
const MARGENS = [0.0, 0.5, 1.0, 1.5, 2.0]

println("\n", "="^78)
println("CUSTO DA CONFORMIDADE: counts sobre a asa sem torção")
println("="^78)
@printf("%-10s", "margem")
for d in DELTAS; @printf("%12s", @sprintf("Δ=%.0f°/est", d)); end
@printf("%12s\n", "sem limite")

tabela = Dict{Tuple{Float64,Float64},Float64}()
for mg in MARGENS
    @printf("%-10s", @sprintf("%.1f°", mg))
    for d in vcat(DELTAS, Inf)
        s = resolve(taxa = (isinf(d) ? nothing : (-d, d)), margem = mg)
        if s === nothing
            @printf("%12s", "inviável")
        else
            c = 1e4*(avalia(torcoes(s[1])).CDff - base.CDff)
            tabela[(mg, d)] = c
            @printf("%12s", @sprintf("%+.1f", c))
        end
    end
    println()
end
println("\n(Δ = salto máximo de torção admitido entre estações vizinhas.")
println(" A coluna sem limite é a solução livre, que chega a saltar 9.5° e")
println(" por isso está fora da validade do VLM e da fabricação.)")

open(joinpath(SAIDA, "custo_conformidade.txt"), "w") do io
    println(io, "custo da conformidade com a FAR 25.203, counts de CDff sobre a asa sem torção")
    println(io, "asa sem torção: CDff $(round(base.CDff, digits = 6)), estola em η $(round(eta_crit(zeros(NV)), digits = 3)) (dentro do aileron, não conforme)")
    println(io, "margem[graus]  delta[graus/estacao]  counts")
    for mg in MARGENS, d in vcat(DELTAS, Inf)
        haskey(tabela, (mg, d)) || continue
        @printf(io, "%6.1f %14s %12.2f\n", mg,
                isinf(d) ? "sem limite" : @sprintf("%.0f", d), tabela[(mg, d)])
    end
end

q = plot(; xlabel = "margem de estol exigida [graus de α]",
         ylabel = "custo em counts de CDff", legend = :topleft,
         title = "preço da conformidade com a FAR 25.203", titlefontsize = 11,
         titlelocation = :left, ESTILO...)
hline!(q, [0.0]; color = "#c3c2b7", linewidth = 0.8, label = "")
for (i, d) in enumerate(vcat(DELTAS, Inf))
    xs = [mg for mg in MARGENS if haskey(tabela, (mg, d))]
    isempty(xs) && continue
    plot!(q, xs, [tabela[(mg, d)] for mg in xs];
          color = ["#2a78d6", "#1baf7a", "#eb6834", "#8a56c4", "#8a8a8a"][i],
          linewidth = 2.2, marker = :circle, markersize = 4,
          linestyle = isinf(d) ? :dash : :solid,
          markerstrokecolor = :auto,
          label = isinf(d) ? "sem limite de salto" : @sprintf("Δ = %.0f°/estação", d))
end
savefig(plot(q; size = (860, 520), dpi = 200),
        joinpath(SAIDA, "custo_conformidade.png"))
println("figura: $SAIDA/custo_conformidade.png")

# --------------------------------------------------------------------
# CANDIDATOS DE PROJETO
#
# Os pontos da tabela acima que combinam ganho de arrasto com torção
# fabricável. Além do salto entre estações interessa a torção total da
# asa: transporte comercial costuma ficar na casa de 4 a 6 graus entre
# raiz e ponta, então uma solução que peça 12 graus é suspeita mesmo
# que o salto local esteja dentro do limite.

const CANDIDATOS = [(0.0, 2.0), (0.0, 3.0), (0.5, 3.0), (0.5, 4.0),
                    (1.0, 4.0)]

println("\n", "="^78)
println("CANDIDATOS DE PROJETO")
println("="^78)
@printf("%-16s %9s %8s %8s %8s %7s %7s\n", "margem / Δ", "CDff", "counts",
        "Oswald", "estol η", "salto", "faixa")
escolhas = Dict{String,Any}()
for (mg, d) in CANDIDATOS
    s = resolve(taxa = (-d, d), margem = mg)
    s === nothing && continue
    x = s[1]; tw = torcoes(x); av = avalia(tw)
    rotulo = @sprintf("%.1f° / %.0f°/est", mg, d)
    escolhas[rotulo] = Dict("margem" => mg, "delta" => d, "torcoes" => tw,
                            "CDff" => av.CDff, "CDvis" => av.CDvis,
                            "counts" => 1e4*(av.CDff - base.CDff),
                            "oswald" => av.e, "estol_eta" => eta_crit(x),
                            "profundor" => av.de)
    @printf("%-16s %9.6f %+8.2f %8.4f %8.3f %6.1f° %6.1f°\n", rotulo, av.CDff,
            1e4*(av.CDff - base.CDff), av.e, eta_crit(x), salto(x),
            maximum(tw) - minimum(tw))
    @printf("   torções  %s\n",
            join([@sprintf("%6.2f", v) for v in tw], ""))
end

# --------------------------------------------------------------------
# TORÇÃO COMO CURVA SUAVE
#
# Limitar o salto entre estações vizinhas controla a amplitude do
# degrau, mas não impede a solução de ondular dentro do limite. Para
# eliminar a ondulação de vez a torção passa a ser descrita por uma
# curva de poucos graus de liberdade:
#
#     tw(η) = c1 η + c2 η² + ... + cp η^p
#
# Não há termo constante porque a raiz é o gauge. Como a torção é
# linear nos coeficientes, o arrasto continua quadrático neles e a
# restrição de estol continua afim: é o mesmo QP convexo, resolvido num
# espaço de dimensão p. Serrilhado deixa de ser representável.

base_poli(p) = [ETAS[i+1]^k for i in 1:NV, k in 1:p]

function resolve_poli(p; margem = MARGEM, taxa = nothing)
    B = base_poli(p)
    Hq = Matrix(Symmetric(B'*HH*B)); Gq = B'*G
    lo, hi = fill(-80.0, p), fill(80.0, p)
    melhor = nothing
    for i_ref in findall(ETA_F .< ETA_AILERON)
        Ae, lce = restricao_estol(i_ref, margem)
        A = vcat(Ae*B, B, -B)                    # estol e limites de torção
        lc = vcat(lce, fill(TW_MIN, NV), fill(-TW_MAX, NV))
        uc = fill(Inf, length(lc))
        if taxa !== nothing
            At, lct, uct = restricao_taxa(taxa[1], taxa[2])
            A = vcat(A, At*B); lc = vcat(lc, lct); uc = vcat(uc, uct)
        end
        for t in (-2.0, -5.0, -9.0), k in 1:min(p, 2)
            z0 = zeros(p); z0[k] = t
            zi = ponto_interior(A, lc, uc, lo .+ 1e-3, hi .- 1e-3, z0)
            zi === nothing && continue
            z = try resolve_qp(Hq, Gq, A, lc, uc, lo, hi, zi) catch; zi end
            any(isnan, z) && continue
            x = B*z
            dtw = diff(torcoes(x))
            viavel = eta_crit(x) < ETA_AILERON &&
                     margem_ail(x) >= margem - 0.05 &&
                     all(x .>= TW_MIN - 1e-6) && all(x .<= TW_MAX + 1e-6) &&
                     (taxa === nothing ||
                      (minimum(dtw) >= taxa[1] - 1e-4 &&
                       maximum(dtw) <= taxa[2] + 1e-4))
            viavel && (melhor === nothing || cd_mod(x) < melhor[2]) &&
                (melhor = (x, cd_mod(x), z))
        end
    end
    return melhor
end

println("\n", "="^78)
println("TORÇÃO COMO CURVA SUAVE tw(η) = c1 η + c2 η² + ...")
println("="^78)
@printf("%-18s %9s %8s %8s %8s %7s %7s\n", "margem / grau", "CDff", "counts",
        "Oswald", "estol η", "salto", "faixa")
for mg in [0.0, 0.5, 1.0, 1.5], p in 2:4
    s = resolve_poli(p; margem = mg)
    s === nothing && continue
    x = s[1]; tw = torcoes(x); av = avalia(tw)
    rotulo = @sprintf("%.1f° / grau %d", mg, p)
    escolhas[rotulo] = Dict("margem" => mg, "grau" => p, "torcoes" => tw,
                            "coef" => s[3], "CDff" => av.CDff,
                            "CDvis" => av.CDvis,
                            "counts" => 1e4*(av.CDff - base.CDff),
                            "oswald" => av.e, "estol_eta" => eta_crit(x),
                            "profundor" => av.de)
    @printf("%-18s %9.6f %+8.2f %8.4f %8.3f %6.1f° %6.1f°\n", rotulo, av.CDff,
            1e4*(av.CDff - base.CDff), av.e, eta_crit(x), salto(x),
            maximum(tw) - minimum(tw))
    @printf("   torções  %s\n", join([@sprintf("%6.2f", v) for v in tw], ""))
end

# --------------------------------------------------------------------
# A MONOTONICIDADE NA MARGEM CERTA
#
# A comparação anterior entre torção livre e monotônica foi feita com
# margem de 1.0°, que é onde toda restrição custa caro. O preço da
# monotonicidade só faz sentido lido contra a margem, porque as duas
# exigências disputam o mesmo recurso: tirar carga da ponta.

println("\n", "="^78)
println("PREÇO DA MONOTONICIDADE, POR MARGEM (counts sobre a asa sem torção)")
println("="^78)
@printf("%-10s %12s %12s %12s %12s\n", "margem", "livre", "suave 3°",
        "monotônica", "cúbica mon.")
for mg in MARGENS
    @printf("%-10s", @sprintf("%.1f°", mg))
    for modo in (:livre, :suave, :mono, :poli_mono)
        s = modo === :livre     ? resolve(margem = mg) :
            modo === :suave     ? resolve(taxa = (-3.0, 3.0), margem = mg) :
            modo === :mono      ? resolve(taxa = (-Inf, 0.0), margem = mg) :
                                  resolve_poli(4; margem = mg,
                                               taxa = (-Inf, 0.0))
        if s === nothing
            @printf("%12s", "inviável")
        else
            @printf("%12s", @sprintf("%+.1f",
                    1e4*(avalia(torcoes(s[1])).CDff - base.CDff)))
        end
    end
    println()
end
println("\n(cúbica mon. = tw(η) polinomial de grau 4 e não crescente:")
println(" é a forma clássica de washout, suave e sem barriga)")

sm = resolve_poli(4; margem = 0.0, taxa = (-Inf, 0.0))
if sm !== nothing
    tw = torcoes(sm[1]); av = avalia(tw)
    println("\nmelhor torção monotônica e suave, margem 0°:")
    @printf("  torções  %s\n", join([@sprintf("%6.2f", v) for v in tw], ""))
    @printf("  CDff %.6f (%+.2f counts)  Oswald %.4f  estol η %.3f  faixa %.1f°\n",
            av.CDff, 1e4*(av.CDff - base.CDff), av.e, eta_crit(sm[1]),
            maximum(tw) - minimum(tw))
    escolhas["monotônica suave 0.0°"] =
        Dict("margem" => 0.0, "grau" => 4, "monotona" => true,
             "torcoes" => tw, "CDff" => av.CDff, "CDvis" => av.CDvis,
             "counts" => 1e4*(av.CDff - base.CDff), "oswald" => av.e,
             "estol_eta" => eta_crit(sm[1]), "profundor" => av.de)
end

println("\n", "="^78)
println("FINALISTAS")
println("="^78)
const FINAIS = [("margem 0.0, suave 3°",   0.0, () -> resolve(taxa = (-3.0, 3.0), margem = 0.0)),
                ("margem 0.0, monotônica", 0.0, () -> resolve(taxa = (-Inf, 0.0), margem = 0.0)),
                ("margem 0.0, cúbica g4",  0.0, () -> resolve_poli(4; margem = 0.0)),
                ("margem 0.0, cúbica mon", 0.0, () -> resolve_poli(4; margem = 0.0, taxa = (-Inf, 0.0))),
                ("margem 0.5, suave 3°",   0.5, () -> resolve(taxa = (-3.0, 3.0), margem = 0.5)),
                ("margem 0.5, monotônica", 0.5, () -> resolve(taxa = (-Inf, 0.0), margem = 0.5)),
                ("margem 0.5, cúbica g4",  0.5, () -> resolve_poli(4; margem = 0.5)),
                ("margem 0.0, mono+Δ3",    0.0, () -> resolve(taxa = (-3.0, 0.0), margem = 0.0)),
                ("margem 0.0, mono+Δ2",    0.0, () -> resolve(taxa = (-2.0, 0.0), margem = 0.0)),
                ("margem 0.5, mono+Δ3",    0.5, () -> resolve(taxa = (-3.0, 0.0), margem = 0.5))]
@printf("%-24s %9s %8s %8s %8s %7s %7s\n", "caso", "CDff", "counts",
        "Oswald", "estol η", "salto", "faixa")
for (nome, mg, f) in FINAIS
    s = f()
    if s === nothing
        @printf("%-24s %9s\n", nome, "inviável"); continue
    end
    x = s[1]; tw = torcoes(x); av = avalia(tw)
    escolhas[nome] = Dict("margem" => mg, "torcoes" => tw, "CDff" => av.CDff,
                          "CDvis" => av.CDvis, "oswald" => av.e,
                          "counts" => 1e4*(av.CDff - base.CDff),
                          "estol_eta" => eta_crit(x), "profundor" => av.de,
                          "salto" => salto(x),
                          "faixa" => maximum(tw) - minimum(tw))
    @printf("%-24s %9.6f %+8.2f %8.4f %8.3f %6.1f° %6.1f°\n", nome, av.CDff,
            1e4*(av.CDff - base.CDff), av.e, eta_crit(x), salto(x),
            maximum(tw) - minimum(tw))
    @printf("   torções  %s\n", join([@sprintf("%6.2f", v) for v in tw], ""))
end

open(joinpath(SAIDA, "candidatos_torcao.json"), "w") do io
    JSON.print(io, Dict("base_CDff" => base.CDff,
                        "base_estol_eta" => eta_crit(zeros(NV)),
                        "etas" => ETAS, "candidatos" => escolhas), 2)
end
println("\ngravado: $SAIDA/candidatos_torcao.json")

# --------------------------------------------------------------------
# FIGURA DE DECISÃO

const DECISAO = ["margem 0.0, livre", "margem 0.0, suave 3°",
                 "margem 0.0, mono+Δ2", "margem 0.5, monotônica"]

# A curva livre mostrada é a da MESMA margem das demais. Comparar formas
# obtidas com margens diferentes na mesma figura levaria a atribuir à
# parametrização um efeito que é da margem.
sl = resolve(margem = 0.0)
if sl !== nothing
    escolhas["margem 0.0, livre"] =
        Dict("torcoes" => torcoes(sl[1]),
             "counts" => 1e4*(avalia(torcoes(sl[1])).CDff - base.CDff))
end

r = plot(; xlabel = "η = 2y/b", ylabel = "torção [graus]", legend = :bottomleft,
         title = "torção: da solução livre à solução construível",
         titlefontsize = 11, titlelocation = :left, xlims = (0, 1), ESTILO...)
vspan!(r, [ETA_AILERON, 0.90]; color = "#dcd9cd", alpha = 0.45,
       linewidth = 0, label = "aileron")
hline!(r, [0.0]; color = "#c3c2b7", linewidth = 0.8, label = "")
for (i, nome) in enumerate(DECISAO)
    haskey(escolhas, nome) || continue
    d = escolhas[nome]
    plot!(r, ETAS, d["torcoes"];
          color = ["#8a8a8a", "#2a78d6", "#1baf7a", "#eb6834"][i],
          linewidth = i == 1 ? 1.8 : 2.6, linestyle = i == 1 ? :dash : :solid,
          marker = :circle, markersize = 4.5, markerstrokewidth = 0,
          label = @sprintf("%s  (%+.1f counts)", nome, d["counts"]))
end
savefig(plot(r; size = (900, 540), dpi = 200),
        joinpath(SAIDA, "torcao_decisao.png"))
println("figura: $SAIDA/torcao_decisao.png")

open(joinpath(SAIDA, "torcao_recomendada.txt"), "w") do io
    d = escolhas["margem 0.0, mono+Δ2"]
    println(io, "PONTO DE PROJETO RECOMENDADO PARA A TORÇÃO DA ASA")
    println(io)
    println(io, "vínculos: torção não crescente da raiz para a ponta,")
    println(io, "          salto máximo de 2 graus entre estações vizinhas,")
    println(io, "          estol iniciando fora do aileron (FAR 25.203),")
    println(io, "          margem de 0 grau de ângulo de ataque")
    println(io)
    @printf(io, "CDff   %.6f  (%+.2f counts contra a asa sem torção)\n",
            d["CDff"], d["counts"])
    @printf(io, "CDvis  %.6f\n", d["CDvis"])
    @printf(io, "Oswald %.4f\n", d["oswald"])
    @printf(io, "estol começa em η = %.3f (aileron de %.2f a 0.90)\n",
            d["estol_eta"], ETA_AILERON)
    @printf(io, "profundor de compensação %.2f graus\n", d["profundor"])
    @printf(io, "salto máximo %.1f graus, faixa total %.1f graus\n",
            d["salto"], d["faixa"])
    println(io)
    println(io, "η        torção [graus]")
    for (e, t) in zip(ETAS, d["torcoes"])
        @printf(io, "%.4f   %7.2f\n", e, t)
    end
    println(io)
    println(io, "ALTERNATIVA MAIS ROBUSTA, com 0.5 grau de margem de estol:")
    a = escolhas["margem 0.5, monotônica"]
    @printf(io, "CDff %.6f (%+.2f counts), torções %s\n", a["CDff"],
            a["counts"], join([@sprintf("%.2f", v) for v in a["torcoes"]], " "))
end
println("gravado: $SAIDA/torcao_recomendada.txt")
