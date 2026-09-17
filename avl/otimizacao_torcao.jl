# Otimização da torção da asa -- estudo próprio do Lab 04.
#
# Rodar de dentro da pasta avl/:   julia otimizacao_torcao.jl
# Saídas em resultados/.
#
# O estudo vai do problema mais simples ao de projeto, de modo que cada
# etapa valide a anterior:
#
#   ETAPA 0  linha de base: a asa sem torção, como está no modelo.
#
#   ETAPA 1  VALIDAÇÃO. Asa isolada, sem restrição, minimizando o arrasto.
#            Existe resposta conhecida: a carga elíptica, com fator de
#            Oswald igual a 1 e arrasto igual a CL²/(π AR). Se a otimização
#            chega lá, o problema está bem posto. É um teste, não projeto.
#
#   ETAPA 2  Aeronave completa e trimada, ainda sem restrição. Entram o
#            downwash da asa sobre a empenagem e a carga de compensação, e
#            o ótimo deixa de ser a asa elíptica: o mínimo é do CONJUNTO.
#
#   ETAPA 3  Projeto. Acrescenta a exigência de estol da FAR 25.203.
#
# CRITÉRIO DE ESTOL
#   A FAR 25.203(a) exige que o comando de rolamento continue eficaz até e
#   durante o estol. A tradução geométrica disso não é uma estação
#   arbitrária: é que o estol NÃO comece na região do aileron, que nesta
#   asa vai de η = 0,56 a η = 0,90. A prática dos transportes a jato é
#   garantir isso pelo perfil da raiz, com clmax menor que o da ponta, e
#   não pela torção; aqui a torção é a única variável disponível, o que
#   torna o custo dessa exigência um resultado do estudo.
#   A margem é dada em ÂNGULO DE ATAQUE: quantos graus a região do aileron
#   ainda aguenta depois que a região interna estola.
#
# ESTRUTURA DO PROBLEMA
#   Para CL fixo a distribuição de sustentação é AFIM nas torções e o
#   arrasto induzido é QUADRÁTICO. Os modelos usados aqui são, por isso,
#   exatos, e todo ótimo é verificado com rodadas do AVL fora deles.

using Printf
using Plots
using JSON
using Optim
using Statistics

gr()

# --------------------------------------------------------------------
# CONFIGURAÇÃO

const AVL   = "./avl.exe"
const BASE  = "aft.avl"
const SAIDA = "resultados"
const DADOS = JSON.parsefile("dados_designtool.json")
const PP    = DADOS["ponto_de_projeto"]

const MACH    = PP["M"]
const CL_PROJ = round(PP["CL"], digits = 4)
const SREF    = PP["Sref"]
const BREF    = DADOS["referencia"]["Bref"]
const SEMI    = BREF/2
const AR      = BREF^2/SREF

# Torções livres: uma por seção da asa, com a raiz fixa em zero. A raiz é
# referência de gauge (somar uma constante a toda a asa é girar a asa
# inteira, o que o ângulo de ataque de equilíbrio absorve).
const ETAS   = [0.0, 0.1011, 0.398, 0.56, 0.90, 1.0]
const LIVRES = 2:6
const NV     = length(LIVRES)

const TW_MIN, TW_MAX = -12.0, 6.0
const EPS_H = 3.0

# Estol
const MACH_BAIXO  = 0.2
const ALFAS_BASE  = (8.0, 14.0)
const ETA_AILERON = 0.56          # raiz do aileron: o estol tem de vir antes
const MARGEM_PROJ = 1.0           # [graus] de margem adotada no projeto
const ETA_LIM   = [0.1011, 0.398, 0.90]
const CLMAX_LIM = [1.774, 1.7985, 1.7338]

const PAL   = ["#2a78d6", "#eb6834", "#1baf7a"]
const INK   = "#0b0b0b"
const INK2  = "#52514e"
const CINZA = "#b9b8b3"
const ESTILO = (framestyle = :axes, gridcolor = "#e1e0d9", gridalpha = 1.0,
                gridlinewidth = 0.7, foreground_color_axis = INK2,
                foreground_color_border = "#c3c2b7",
                foreground_color_text = INK2, tickfontsize = 9,
                guidefontsize = 10, legendfontsize = 8,
                background_color = :white)

mkpath(SAIDA)

# --------------------------------------------------------------------
# INFRAESTRUTURA

function roda_avl(cmds::AbstractString)
    tmp = "_ot_cmds.txt"
    write(tmp, cmds)
    try
        return read(pipeline(`$AVL`, stdin = tmp), String)
    finally
        rm(tmp, force = true)
    end
end

function num(txt::AbstractString, pat::Regex)
    m = collect(eachmatch(pat, txt))
    isempty(m) && return NaN
    return parse(Float64, m[end].captures[1])
end

"Variante do modelo com as torções dadas. `so_asa` mantém só a asa."
function escreve(destino, tw; so_asa = false)
    linhas = readlines(BASE)
    n = length(linhas)
    inicios = [i for i in 1:n if strip(linhas[i]) in ("SURFACE", "BODY")]
    i_asa = 0
    for i in inicios
        if strip(linhas[i]) == "SURFACE"
            j = i + 1
            while j <= n && (isempty(strip(linhas[j])) ||
                             startswith(strip(linhas[j]), "#"))
                j += 1
            end
            j <= n && strip(linhas[j]) == "Wing" && (i_asa = i)
        end
    end
    i_asa == 0 && error("não achei a superfície Wing em $BASE")
    prox = findfirst(>(i_asa), inicios)
    fim_asa = prox === nothing ? n : inicios[prox] - 1
    faixa = so_asa ? (1:fim_asa) : (1:n)

    saida = String[]
    i_sec, espera = 0, false
    for i in faixa
        ln = linhas[i]
        t = strip(ln)
        dentro = (i >= i_asa && i <= fim_asa)
        if dentro && t == "SECTION"
            espera = true
        elseif dentro && espera && !isempty(t) && !startswith(t, "#")
            i_sec += 1
            p = split(t)
            p[5] = @sprintf("%.4f", tw[i_sec])
            push!(saida, join(p, "  "))
            espera = false
            continue
        end
        push!(saida, ln)
    end
    i_sec == length(tw) || error("apliquei $i_sec torções de $(length(tw))")
    write(destino, join(saida, "\n") * "\n")
end

function avalia(tw; so_asa = false, trim = false, faixas = false)
    arq = "_ot.avl"
    escreve(arq, tw; so_asa = so_asa)
    cmds = ["load $arq", "oper", "m", "mn $MACH", ""]
    trim && push!(cmds, "d2 pm 0")
    append!(cmds, ["a c $CL_PROJ", "x"])
    faixas && append!(cmds, ["fs", ""])
    append!(cmds, ["", "quit"])
    s = roda_avl(join(cmds, "\n") * "\n")
    rm(arq, force = true)
    return (CDff  = num(s, r"CDff\s*=\s*([-\d.Ee+]+)"),
            CDvis = num(s, r"CDvis\s*=\s*([-\d.Ee+]+)"),
            CDtot = num(s, r"CDtot\s*=\s*([-\d.Ee+]+)"),
            CL    = num(s, r"CLtot\s*=\s*([-\d.]+)"),
            e     = num(s, r"\se =\s+([-\d.]+)"),
            alpha = num(s, r"Alpha\s*=\s*([-\d.]+)"),
            de    = num(s, r"elevator\s*=\s*([-\d.]+)"),
            saida = s)
end

function faixas_asa(saida)
    tr = split(saida, r"Surface # 1\s+Wing")[2]
    bl = String(split(tr, "Surface # 2")[1])
    L = [[parse(Float64, v) for v in split(strip(m.match))[2:end]]
         for m in eachmatch(r"^\s*\d+(?:\s+[-\dEe.+]+){12}\s*$"m, bl)]
    isempty(L) && error("não achei as faixas da asa")
    d = reduce(hcat, L)'
    return (eta = d[:, 1]./SEMI, ccl = d[:, 4], cl_norm = d[:, 6])
end

torcoes(x) = (tw = zeros(length(ETAS)); tw[LIVRES] .= x; tw)

# O CDtot que o AVL imprime soma o induzido de CAMPO PRÓXIMO, que é a
# medida ruidosa. O total usado aqui é o coerente com o objetivo:
# parasita mais induzido do plano de Trefftz.
cdtot_ff(av) = av.CDvis + av.CDff
carga_norm(fx, CL) = fx.ccl ./ (4*SREF*CL/(pi*BREF))
eliptica(eta) = sqrt.(max.(0.0, 1 .- eta.^2))

# --------------------------------------------------------------------
# MODELO EXATO DO ARRASTO (quadrático nas torções, para CL fixo)

function modelo_arrasto(; so_asa, trim)
    b = avalia(torcoes(zeros(NV)); so_asa = so_asa, trim = trim)
    f0 = b.CDff
    fp, fm = zeros(NV), zeros(NV)
    for j in 1:NV
        e = zeros(NV); e[j] = EPS_H
        fp[j] = avalia(torcoes(e); so_asa = so_asa, trim = trim).CDff
        fm[j] = avalia(torcoes(-e); so_asa = so_asa, trim = trim).CDff
    end
    g = (fp .- fm)./(2*EPS_H)
    H = zeros(NV, NV)
    for j in 1:NV
        H[j, j] = (fp[j] + fm[j] - 2*f0)/EPS_H^2
    end
    for j in 1:NV-1, k in j+1:NV
        e = zeros(NV); e[j] = e[k] = EPS_H
        fjk = avalia(torcoes(e); so_asa = so_asa, trim = trim).CDff
        H[j, k] = H[k, j] = (fjk - fp[j] - fp[k] + f0)/EPS_H^2
    end
    return (f0 = f0, g = g, H = H, cdvis = b.CDvis)
end
cd_mod(m, x) = m.f0 + m.g'x + 0.5*x'*m.H*x

# --------------------------------------------------------------------
# MODELO EXATO DO ESTOL (afim nas torções e no ângulo de ataque)

function modelo_estol()
    function amostra(x, alfa)
        escreve("_ot.avl", torcoes(x))
        s = roda_avl(join(["load _ot.avl", "oper", "m", "mn $MACH_BAIXO", "",
                           "d2 d2 0", "a a $alfa", "x", "fs", "", "",
                           "quit"], "\n") * "\n")
        rm("_ot.avl", force = true)
        return faixas_asa(s)
    end
    d1 = amostra(zeros(NV), ALFAS_BASE[1])
    d2 = amostra(zeros(NV), ALFAS_BASE[2])
    b = (d2.cl_norm .- d1.cl_norm)./(ALFAS_BASE[2] - ALFAS_BASE[1])
    a = d1.cl_norm .- b.*ALFAS_BASE[1]
    C = zeros(NV, length(d1.eta))
    for j in 1:NV
        e = zeros(NV); e[j] = 1.0
        C[j, :] = amostra(e, ALFAS_BASE[1]).cl_norm .- d1.cl_norm
    end
    function peso(η)
        η <= ETA_LIM[1] && return (1.0, 0.0, 0.0)
        η >= ETA_LIM[3] && return (0.0, 0.0, 1.0)
        if η <= ETA_LIM[2]
            w = (η - ETA_LIM[1])/(ETA_LIM[2] - ETA_LIM[1]); return (1-w, w, 0.0)
        end
        w = (η - ETA_LIM[2])/(ETA_LIM[3] - ETA_LIM[2]); return (0.0, 1-w, w)
    end
    lim = [sum(CLMAX_LIM[i]*w for (i, w) in enumerate(peso(η))) for η in d1.eta]
    return (eta = d1.eta, a = a, b = b, C = C, lim = lim)
end

"Ângulo de ataque em que cada faixa alcança o clmax do seu perfil."
alfa_estol(me, x) = (me.lim .- me.a .- me.C'x)./me.b

"""
Margem de estol do aileron: quantos graus de ângulo de ataque a região do
aileron ainda aguenta depois que a região interna estola. Positiva
significa que o estol começa para dentro do aileron, que é o que a
FAR 25.203 exige para preservar o comando de rolamento.
"""
function margem_aileron(me, x)
    al = alfa_estol(me, x)
    dentro = me.eta .< ETA_AILERON
    return minimum(al[.!dentro]) - minimum(al[dentro])
end

"Estação onde o estol começa."
eta_critico(me, x) = me.eta[argmin(alfa_estol(me, x))]

"""
Restrição de estol na forma linear.

A exigência escrita como diferença de mínimos não é convexa. Elege-se
então uma faixa interna de referência e exige-se que TODA faixa do aileron
estole pelo menos `margem` graus depois dela: cada exigência é afim nas
torções. A escolha da referência é varrida e verificada, de modo que a
reformulação não restringe o resultado.
"""
function restricao_linear(me, i_ref, margem)
    ext = findall(me.eta .>= ETA_AILERON)
    d = [(me.lim[j]-me.a[j])/me.b[j] - (me.lim[i_ref]-me.a[i_ref])/me.b[i_ref]
         for j in ext]
    W = hcat([-me.C[:, j]./me.b[j] .+ me.C[:, i_ref]./me.b[i_ref] for j in ext]...)
    return z -> sum(min.(d .+ vec(W'z) .- margem, 0.0).^2)
end

# --------------------------------------------------------------------
# OTIMIZAÇÃO COM REGISTRO DO CAMINHO

function otimiza(m; pen = nothing, x0 = zeros(NV))
    caminho = [copy(x0)]
    lo, hi = fill(TW_MIN, NV), fill(TW_MAX, NV)
    x = copy(x0)
    for w in (pen === nothing ? [0.0] : [1e2, 1e3, 1e4, 1e5, 1e6])
        f(z) = 1e4*cd_mod(m, z) + (pen === nothing ? 0.0 : w*pen(z))
        r = optimize(f, lo, hi, x, Fminbox(LBFGS()),
                     Optim.Options(store_trace = true, extended_trace = true,
                                   iterations = 300, g_tol = 1e-11))
        for t in Optim.trace(r)
            haskey(t.metadata, "x") && push!(caminho, copy(t.metadata["x"]))
        end
        x = Optim.minimizer(r)
        push!(caminho, copy(x))
    end
    return x, caminho
end

"Etapa 3: varre a faixa de referência e devolve o melhor ótimo viável."
function otimiza_com_estol(m, me, margem; partidas)
    melhor = nothing
    for i_ref in findall(me.eta .< ETA_AILERON), x0 in partidas
        pen = restricao_linear(me, i_ref, margem)
        x, cam = otimiza(m; pen = pen, x0 = x0)
        ok = eta_critico(me, x) < ETA_AILERON &&
             margem_aileron(me, x) >= margem - 0.05
        f = cd_mod(m, x)
        if ok && (melhor === nothing || f < melhor.f)
            melhor = (x = x, cam = cam, f = f, eta_ref = me.eta[i_ref])
        end
    end
    return melhor
end

# --------------------------------------------------------------------
# ANIMAÇÃO DO CAMINHO

function anima(caminho, arquivo; so_asa, trim, titulo, n_quadros = 24)
    idx = unique(round.(Int, range(1, length(caminho), length = n_quadros)))
    dados = map(caminho[idx]) do x
        av = avalia(torcoes(x); so_asa = so_asa, trim = trim, faixas = true)
        fx = faixas_asa(av.saida)
        (tw = torcoes(x), eta = fx.eta, carga = carga_norm(fx, av.CL),
         CDff = av.CDff, e = av.e)
    end
    lo = minimum(minimum(d.tw) for d in dados) - 0.5
    hi = maximum(maximum(d.tw) for d in dados) + 0.5
    cd0 = dados[1].CDff
    anim = @animate for (k, d) in enumerate(dados)
        p1 = plot(; xlabel = "η = 2y/b", ylabel = "torção [graus]",
                  title = "torção", titlefontsize = 10, titlelocation = :left,
                  legend = false, xlims = (0, 1), ylims = (lo, hi), ESTILO...)
        hline!(p1, [0.0]; color = "#c3c2b7", linewidth = 0.8)
        plot!(p1, ETAS, d.tw; color = PAL[1], linewidth = 2.6,
              marker = :circle, markersize = 5, markerstrokecolor = PAL[1])
        p2 = plot(; xlabel = "η = 2y/b", ylabel = "carga normalizada",
                  title = "sustentação", titlefontsize = 10,
                  titlelocation = :left, legend = :bottomleft,
                  xlims = (0, 1), ylims = (0, 1.15), ESTILO...)
        plot!(p2, d.eta, eliptica(d.eta); color = INK, linewidth = 2,
              linestyle = :dash, label = "elíptica")
        plot!(p2, d.eta, d.carga; color = PAL[1], linewidth = 2.6, label = "asa")
        plot(p1, p2; layout = (1, 2), size = (1100, 460), dpi = 140,
             plot_title = @sprintf("%s | iteração %d/%d | CDff %.6f (%+.1f counts) | e %.4f",
                                   titulo, k, length(dados), d.CDff,
                                   1e4*(d.CDff - cd0), d.e),
             plot_titlefontsize = 10, left_margin = 6Plots.mm,
             bottom_margin = 6Plots.mm)
    end
    gif(anim, joinpath(SAIDA, arquivo), fps = 3, show_msg = false)
    println("  animação: $SAIDA/$arquivo  ($(length(dados)) quadros)")
end

# ====================================================================
# ETAPA 0 -- LINHA DE BASE
# ====================================================================
println("="^78)
println("ETAPA 0 -- LINHA DE BASE: a asa sem torção")
println("="^78)
base_asa  = avalia(torcoes(zeros(NV)); so_asa = true, faixas = true)
base_full = avalia(torcoes(zeros(NV)); trim = true, faixas = true)
piso = CL_PROJ^2/(pi*AR)
@printf("  asa isolada          CDff %.6f   Oswald %.4f\n",
        base_asa.CDff, base_asa.e)
@printf("  aeronave trimada     CDff %.6f   Oswald %.4f   profundor %.2f°\n",
        base_full.CDff, base_full.e, base_full.de)
@printf("  CDvis (parasita, constante sob torção) %.6f\n", base_full.CDvis)
@printf("  piso analítico CL²/(π AR) = %.6f  (carga elíptica)\n", piso)

# ====================================================================
# ETAPA 1 -- VALIDAÇÃO: ASA ISOLADA, IRRESTRITO
# ====================================================================
println("\n", "="^78)
println("ETAPA 1 -- VALIDAÇÃO: asa isolada, sem restrição")
println("="^78)
println("  A resposta é conhecida: carga elíptica, Oswald 1, CDff no piso.")

m1 = modelo_arrasto(so_asa = true, trim = false)
x1, cam1 = otimiza(m1)
ot1 = avalia(torcoes(x1); so_asa = true, faixas = true)
@printf("\n  torções [graus]: %s\n", join([@sprintf("%6.2f", v) for v in torcoes(x1)], " "))
@printf("  CDff   %.6f -> %.6f   (%+.2f counts)\n",
        base_asa.CDff, ot1.CDff, 1e4*(ot1.CDff - base_asa.CDff))
@printf("  Oswald %.4f   -> %.4f\n", base_asa.e, ot1.e)
@printf("  distância ao piso analítico: %+.2f counts\n", 1e4*(ot1.CDff - piso))
@printf("  o modelo previu %.6f e o AVL deu %.6f (erro %.2f counts)\n",
        cd_mod(m1, x1), ot1.CDff, 1e4*abs(cd_mod(m1, x1) - ot1.CDff))
ok1 = abs(ot1.e - 1) < 0.01 && abs(ot1.CDff - piso) < 5e-6
println(ok1 ? "\n  >> VALIDADO: Oswald ~ 1 e arrasto no piso analítico." :
              "\n  >> ATENÇÃO: não bateu a resposta conhecida.")
println("\n  gerando a animação da convergência à elíptica...")
anima(cam1, "evolucao_1_asa_isolada.gif"; so_asa = true, trim = false,
      titulo = "etapa 1: asa isolada, sem restrição")

# ====================================================================
# ETAPA 2 -- AERONAVE COMPLETA, IRRESTRITO
# ====================================================================
println("\n", "="^78)
println("ETAPA 2 -- AERONAVE COMPLETA, trimada, sem restrição")
println("="^78)

m2 = modelo_arrasto(so_asa = false, trim = true)
x2, cam2 = otimiza(m2)
ot2 = avalia(torcoes(x2); trim = true, faixas = true)
@printf("  torções [graus]: %s\n", join([@sprintf("%6.2f", v) for v in torcoes(x2)], " "))
@printf("  CDff   %.6f -> %.6f   (%+.2f counts)\n",
        base_full.CDff, ot2.CDff, 1e4*(ot2.CDff - base_full.CDff))
@printf("  CDtot (CDvis + CDff) %.6f -> %.6f   (%+.2f counts)\n",
        cdtot_ff(base_full), cdtot_ff(ot2),
        1e4*(cdtot_ff(ot2) - cdtot_ff(base_full)))
@printf("  Oswald %.4f   -> %.4f\n", base_full.e, ot2.e)
@printf("  profundor de trimagem: %.2f -> %.2f graus\n", base_full.de, ot2.de)

fx1, fx2 = faixas_asa(ot1.saida), faixas_asa(ot2.saida)
dv1 = maximum(abs.(carga_norm(fx1, ot1.CL) .- eliptica(fx1.eta)))
dv2 = maximum(abs.(carga_norm(fx2, ot2.CL) .- eliptica(fx2.eta)))
@printf("\n  desvio da elíptica: asa isolada %.3f, aeronave completa %.3f\n",
        dv1, dv2)
println("  A aeronave completa não vai para a asa elíptica, e não deveria:")
println("  a empenagem carrega para compensar a arfagem e opera no downwash")
println("  da asa, de modo que o mínimo é do conjunto, não de cada parte.")

# ====================================================================
# ETAPA 3 -- PROJETO: COM A EXIGÊNCIA DE ESTOL DA FAR 25.203
# ====================================================================
println("\n", "="^78)
println("ETAPA 3 -- PROJETO: estol antes do aileron (FAR 25.203)")
println("="^78)

me = modelo_estol()
@printf("  aileron: de η = %.2f a 0,90\n", ETA_AILERON)
@printf("  sem torção:        estol em η = %.3f, margem do aileron %+.2f°\n",
        eta_critico(me, zeros(NV)), margem_aileron(me, zeros(NV)))
@printf("  ótimo da etapa 2:  estol em η = %.3f, margem do aileron %+.2f°\n",
        eta_critico(me, x2), margem_aileron(me, x2))
println("  (margem negativa significa que o estol começa DENTRO do aileron,")
println("   o que a FAR 25.203 não admite: o rolamento perde eficácia)\n")

partidas = [zeros(NV), copy(x2), [2.0, 1.0, -2.0, -5.0, -8.0],
            [0.0, 0.0, -3.0, -6.0, -9.0], [4.0, 3.0, -1.0, -4.0, -7.0]]

println("  custo da margem de estol exigida:")
@printf("  %-10s %10s %10s %9s %10s   %s\n", "margem", "CDff", "counts",
        "Oswald", "estol η", "profundor")
sols = Dict{Float64,Any}()
for mg in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0)
    s = otimiza_com_estol(m2, me, mg; partidas = partidas)
    if s === nothing
        @printf("  %-10.1f %10s\n", mg, "inviável"); continue
    end
    av = avalia(torcoes(s.x); trim = true, faixas = true)
    sols[mg] = (s = s, av = av)
    @printf("  %-10.1f %10.6f %+10.2f %9.4f %10.3f   %8.2f°\n", mg, av.CDff,
            1e4*(av.CDff - base_full.CDff), av.e, eta_critico(me, s.x), av.de)
end

esc = sols[MARGEM_PROJ]
x3, cam3, ot3 = esc.s.x, esc.s.cam, esc.av
@printf("\n  margem adotada no projeto: %.1f grau de ângulo de ataque\n",
        MARGEM_PROJ)
@printf("  torções [graus]: %s\n", join([@sprintf("%6.2f", v) for v in torcoes(x3)], " "))
@printf("  CDff  %.6f  (%+.2f counts contra a asa sem torção)\n",
        ot3.CDff, 1e4*(ot3.CDff - base_full.CDff))
@printf("  CDtot (CDvis + CDff) %.6f  (%+.2f counts)\n", cdtot_ff(ot3),
        1e4*(cdtot_ff(ot3) - cdtot_ff(base_full)))
@printf("  estol em η = %.3f com margem de %+.2f graus\n",
        eta_critico(me, x3), margem_aileron(me, x3))
@printf("  custo da exigência de estol: %+.2f counts sobre a etapa 2\n",
        1e4*(ot3.CDff - ot2.CDff))

println("\n  gerando a animação do caminho da etapa 3...")
anima(cam3, "evolucao_3_com_estol.gif"; so_asa = false, trim = true,
      titulo = "etapa 3: aeronave completa com estol restrito")

# ====================================================================
# O QUE O PERFIL DA RAIZ RESOLVERIA
# ====================================================================
# Posicionar o estol pela torção sai caro nesta asa porque os clmax dos
# três perfis do Lab 03 são quase iguais e, pior, o MENOR deles é o da
# ponta: raiz 1,774, meio 1,799, ponta 1,734. Aquela otimização buscou
# arrasto com restrição de clmax MÍNIMO, então nada empurrou o clmax da
# raiz para baixo. É por isso que os transportes combinam três coisas para
# garantir estol de raiz: torção, perfil de raiz com clmax menor, e
# dispositivos de bordo de ataque diferentes na parte interna e externa.
#
# A conta abaixo responde de quanto precisaria ser o clmax interno para
# que, mantendo a torção ÓTIMA DE ARRASTO da etapa 2, o estol já começasse
# para dentro do aileron com a margem de projeto. Se isso for viável no
# próximo ciclo de perfis, a torção fica livre para fazer só arrasto.

println("\n", "="^78)
println("O QUE O PERFIL DA RAIZ RESOLVERIA")
println("="^78)
al2 = alfa_estol(me, x2)
dentro = me.eta .< ETA_AILERON
alfa_ail = minimum(al2[.!dentro])
@printf("  com a torção da etapa 2, a região do aileron estola em α = %.2f°\n",
        alfa_ail)
@printf("  para o estol vir antes com %.1f° de margem, o clmax interno teria\n",
        MARGEM_PROJ)
println("  de cair para:")
@printf("  %-10s %12s %12s %10s\n", "η", "clmax atual", "clmax alvo", "redução")
for k in findall(dentro)
    me.eta[k] > 0.45 && continue
    alvo = me.a[k] + me.C[:, k]'x2 + me.b[k]*(alfa_ail - MARGEM_PROJ)
    @printf("  %-10.3f %12.3f %12.3f %9.1f%%\n", me.eta[k], me.lim[k], alvo,
            100*(alvo/me.lim[k] - 1))
end
@printf("\n  ganho potencial: a torção voltaria a ser a da etapa 2 e o arrasto\n")
@printf("  cairia de %.6f para %.6f, ou seja %.1f counts.\n",
        ot3.CDff, ot2.CDff, 1e4*(ot3.CDff - ot2.CDff))

# ====================================================================
# FIGURA DE SÍNTESE
# ====================================================================
etapas = [("1: asa isolada", x1, cam1, m1, ot1, PAL[1]),
          ("2: aeronave completa", x2, cam2, m2, ot2, PAL[2]),
          ("3: com estol restrito", x3, cam3, m2, ot3, PAL[3])]

p1 = plot(; xlabel = "iteração", ylabel = "CDff em counts",
          title = "(a) caminho da otimização", titlefontsize = 10,
          titlelocation = :left, legend = :topright, ESTILO...)
for (nome, _, cam, m, _, cor) in etapas
    y = [1e4*cd_mod(m, x) for x in cam]
    plot!(p1, 0:length(y)-1, y; color = cor, linewidth = 2, label = nome)
end

p2 = plot(; xlabel = "η = 2y/b", ylabel = "torção [graus]",
          title = "(b) torção ótima de cada etapa", titlefontsize = 10,
          titlelocation = :left, legend = :bottomleft, xlims = (0, 1),
          ESTILO...)
hline!(p2, [0.0]; color = "#c3c2b7", linewidth = 0.8, label = "")
for (nome, x, _, _, _, cor) in etapas
    plot!(p2, ETAS, torcoes(x); color = cor, linewidth = 2.2, marker = :circle,
          markersize = 5, markerstrokecolor = cor, label = nome)
end

p3 = plot(; xlabel = "η = 2y/b", ylabel = "carga normalizada",
          title = "(c) distribuição de sustentação", titlefontsize = 10,
          titlelocation = :left, legend = :bottomleft, xlims = (0, 1),
          ESTILO...)
fxb = faixas_asa(base_full.saida)
plot!(p3, fxb.eta, eliptica(fxb.eta); color = INK, linewidth = 2,
      linestyle = :dash, label = "elíptica")
plot!(p3, fxb.eta, carga_norm(fxb, base_full.CL); color = CINZA,
      linewidth = 2, label = "sem torção")
for (nome, _, _, _, ot, cor) in etapas
    fx = faixas_asa(ot.saida)
    plot!(p3, fx.eta, carga_norm(fx, ot.CL); color = cor, linewidth = 2.2,
          label = nome)
end

p4 = plot(; xlabel = "margem de estol exigida [graus de α]",
          ylabel = "CDff em counts", legend = :topleft,
          title = "(d) custo da exigência de estol", titlefontsize = 10,
          titlelocation = :left, ESTILO...)
mg = sort(collect(keys(sols)))
plot!(p4, mg, [1e4*sols[k].av.CDff for k in mg]; color = PAL[3],
      linewidth = 2.4, marker = :circle, markersize = 5,
      markerstrokecolor = PAL[3], label = "ótimo com restrição")
hline!(p4, [1e4*base_full.CDff]; color = CINZA, linewidth = 2,
       linestyle = :dash, label = "asa sem torção")
hline!(p4, [1e4*ot2.CDff]; color = PAL[2], linewidth = 2, linestyle = :dot,
       label = "ótimo sem restrição (etapa 2)")
scatter!(p4, [MARGEM_PROJ], [1e4*ot3.CDff]; color = PAL[3], markersize = 10,
         marker = :star5, markerstrokecolor = INK, label = "adotado")

savefig(plot(p1, p2, p3, p4; layout = (2, 2), size = (1300, 900), dpi = 200,
             left_margin = 6Plots.mm, bottom_margin = 6Plots.mm),
        joinpath(SAIDA, "otimizacao_torcao.png"))
println("figura: $SAIDA/otimizacao_torcao.png")

# ====================================================================
# RESUMO
# ====================================================================
println("\n", "="^78)
println("RESUMO")
println("="^78)
@printf("%-24s %10s %10s %9s %9s\n", "caso", "CDff", "counts", "Oswald", "estol η")
@printf("%-24s %10.6f %10s %9.4f %9.3f\n", "sem torção", base_full.CDff, "-",
        base_full.e, eta_critico(me, zeros(NV)))
@printf("%-24s %10.6f %+10.2f %9.4f %9s\n", "1: asa isolada", ot1.CDff,
        1e4*(ot1.CDff - base_asa.CDff), ot1.e, "-")
@printf("%-24s %10.6f %+10.2f %9.4f %9.3f\n", "2: completa irrestrito",
        ot2.CDff, 1e4*(ot2.CDff - base_full.CDff), ot2.e, eta_critico(me, x2))
@printf("%-24s %10.6f %+10.2f %9.4f %9.3f\n", "3: projeto final", ot3.CDff,
        1e4*(ot3.CDff - base_full.CDff), ot3.e, eta_critico(me, x3))

open(joinpath(SAIDA, "torcao_otimizada.json"), "w") do io
    JSON.print(io, Dict(
        "etas" => ETAS,
        "criterio_estol" => Dict(
            "norma" => "FAR 25.203(a): rolamento eficaz até e durante o estol",
            "eta_aileron" => ETA_AILERON,
            "margem_adotada_graus" => MARGEM_PROJ),
        "sem_torcao" => Dict("CDff" => base_full.CDff, "e" => base_full.e,
                             "eta_estol" => eta_critico(me, zeros(NV))),
        "etapa_1_asa_isolada" => Dict("torcoes" => torcoes(x1),
                                      "CDff" => ot1.CDff, "e" => ot1.e,
                                      "piso_analitico" => piso),
        "etapa_2_completa" => Dict("torcoes" => torcoes(x2),
                                   "CDff" => ot2.CDff, "e" => ot2.e,
                                   "eta_estol" => eta_critico(me, x2)),
        "etapa_3_projeto" => Dict("torcoes" => torcoes(x3),
                                  "CDff" => ot3.CDff,
                                  "CDtot_ff" => cdtot_ff(ot3),
                                  "e" => ot3.e,
                                  "eta_estol" => eta_critico(me, x3),
                                  "margem_aileron" => margem_aileron(me, x3)),
        "custo_da_margem" => Dict(string(k) => sols[k].av.CDff for k in mg)), 2)
end
println("\ngravado: $SAIDA/torcao_otimizada.json")
