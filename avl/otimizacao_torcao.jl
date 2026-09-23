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
#   ETAPA 1  VALIDAÇÃO com a asa isolada, já feita e registrada no README:
#            a spline monotônica chegou à carga elíptica, Oswald 1,0035.
#            Não é repetida aqui.
#
#   ETAPA 2  Aeronave completa e compensada, ainda sem restrição. Entram o
#            downwash da asa sobre a empenagem e a carga de compensação, e
#            o ótimo deixa de ser a asa elíptica: o mínimo é do CONJUNTO.
#
#   ETAPA 3  Projeto. Acrescenta a exigência de estol da FAR 25.203.
#
# PARAMETRIZAÇÃO DA TORÇÃO
#   A torção NÃO é livre por estação. Ela é descrita por poucos pontos de
#   controle ao longo da envergadura e interpolada por uma spline cúbica
#   de Hermite com as inclinações de Fritsch-Carlson, a PCHIP. Essa spline
#   tem a propriedade que interessa aqui: com valores de controle
#   monótonos a curva interpolada é monótona, sem ultrapassagem entre nós.
#   Uma spline cúbica natural não garante isso e pode fazer barriga.
#
#   A monotonicidade é imposta pela própria variável de projeto. Em vez de
#   otimizar os valores de torção, otimizam-se os DECREMENTOS entre nós
#   consecutivos, com a raiz em zero por gauge:
#       t_1 = 0,   t_k = -(s_1 + ... + s_{k-1}),   s_k >= 0
#   Isso é não crescente para qualquer s >= 0, de modo que a exigência
#   vira uma caixa simples e não precisa de penalidade nem multiplicador.
#   A raiz é a estação de maior incidência, que é o washout clássico.
#
#   Torções livres por estação foram abandonadas: com elas o ótimo do
#   modelo é serrilhado, com saltos de quase dez graus entre estações
#   vizinhas, porque o arrasto induzido quase não responde aos modos de
#   alta frequência da torção enquanto o critério de estol responde muito.
#   O ótimo serrilhado é o ótimo verdadeiro DO MODELO, não ruído numérico,
#   mas explora uma região onde o próprio método de malha de vórtices
#   deixa de valer e onde a asa não é fabricável.
#
# CRITÉRIO DE ESTOL
#   A FAR 25.203(a) exige que o comando de rolamento continue eficaz até e
#   durante o estol, o que se traduz em o estol não começar na região do
#   aileron, que nesta asa vai de η = 0,56 a η = 0,90. A exigência usada é
#   de POSIÇÃO: a primeira faixa a estolar tem de estar para dentro de uma
#   estação η_lim. O projeto usa η_lim = 0,56, a raiz do aileron, que é o
#   que a FAR pede; a varredura mostra o custo de exigir posições mais
#   para dentro.
#
#   DESEMPATE. Sem nenhuma folga, o ótimo sempre cai num empate: a faixa
#   interna e uma faixa do aileron estolam no mesmo ângulo, e a exigência
#   vale só no papel. Medido direto no AVL, um projeto assim estolou
#   primeiro em η = 0,842, dentro do aileron, porque o modelo afim de estol
#   é um pouco otimista perto do estol (entre 0,11 e 0,19 grau nos casos
#   medidos). A folga de 0,2 grau exigida aqui não é margem de projeto: é
#   o erro medido do modelo, para que a posição do estol valha também na
#   medição direta, que o verifica_torcao.jl confere.
#
# ESTRUTURA DO PROBLEMA
#   Para CL fixo a distribuição de sustentação é AFIM nas torções de
#   estação e o arrasto induzido é QUADRÁTICO nelas. Os modelos abaixo são
#   construídos nesse espaço de estações e por isso são exatos; a spline
#   entra depois, como um mapa dos poucos parâmetros de projeto para as
#   torções de estação. Todo ótimo é verificado com rodadas do AVL.

using Printf
using Plots
using JSON
using Optim
using Random
using Statistics
using LinearAlgebra

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

# As estações são lidas do próprio modelo, de modo que refinar a asa
# (refina_asa.py) não exige tocar neste script.
function etas_da_asa(arquivo)
    linhas = readlines(arquivo)
    n = length(linhas)
    ini = [i for i in 1:n if strip(linhas[i]) in ("SURFACE", "BODY")]
    i_asa = 0
    for i in ini
        strip(linhas[i]) == "SURFACE" || continue
        j = i + 1
        while j <= n && (isempty(strip(linhas[j])) ||
                         startswith(strip(linhas[j]), "#")); j += 1; end
        j <= n && strip(linhas[j]) == "Wing" && (i_asa = i)
    end
    prox = findfirst(>(i_asa), ini)
    fim = prox === nothing ? n : ini[prox] - 1
    ys, espera = Float64[], false
    for i in i_asa:fim
        t = split(linhas[i], "#")[1] |> strip
        if strip(linhas[i]) == "SECTION"
            espera = true
        elseif espera && !isempty(t)
            y = parse(Float64, split(t)[2]); espera = false
            # a seção do winglet repete o y da ponta e só sobe em z: ela não
            # é estação de torção, acompanha a incidência da ponta
            (isempty(ys) || y > ys[end] + 1e-6) && push!(ys, y)
        end
    end
    return ys ./ ys[end]
end

const ETAS   = etas_da_asa(BASE)
const LIVRES = 2:length(ETAS)
const NV     = length(LIVRES)

const TW_MIN, TW_MAX = -12.0, 6.0
const TW_TOTAL = 12.0             # torção total admitida entre raiz e ponta
const DEC_MAX  = 6.0              # decremento máximo por intervalo de nós
const EPS_H = 3.0

# Nós da spline, escolhidos pelo estudo mais adiante. O par 0,48 e 0,56 é
# o que faz diferença: 0,56 é a raiz do aileron e 0,48 é a estação logo
# antes dela. Com os dois a curva pode ficar plana até 0,48 e virar
# depressa em seguida, de modo que o washout não vaza para a região que
# precisa estolar primeiro. Sem o nó em 0,48 a mesma exigência de estol
# custa vinte counts a mais, ainda que para arrasto puro a colocação dos
# nós quase não importe.
const NOS_PROJ = [0.0, 0.30, 0.48, 0.56, 0.70, 0.85, 1.0]

# Estol
const MACH_BAIXO  = 0.2
const ALFAS_BASE  = (8.0, 14.0)
const ETA_AILERON = 0.56          # raiz do aileron, só para as figuras
const ETA_ESTOL_PROJ = 0.56       # o estol tem de começar para dentro daqui
const DESEMPATE = 0.2             # [graus] erro medido do modelo de estol
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
            p[5] = @sprintf("%.4f", tw[min(i_sec, length(tw))])
            push!(saida, join(p, "  "))
            espera = false
            continue
        end
        push!(saida, ln)
    end
    i_sec >= length(tw) || error("apliquei $i_sec torções de $(length(tw))")
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
    d = d[d[:, 1] .< SEMI - 1e-3, :]   # sem as faixas do winglet, todas em y = b/2
    return (eta = d[:, 1]./SEMI, ccl = d[:, 4], cl_norm = d[:, 6])
end

"Torções de estação a partir dos valores livres, com a raiz em zero."
torcoes(x) = (tw = zeros(length(ETAS)); tw[LIVRES] .= x; tw)
livres(tw) = tw[LIVRES]

# O CDtot que o AVL imprime soma o induzido de CAMPO PRÓXIMO, que é a
# medida ruidosa. O total usado aqui é o coerente com o objetivo:
# parasita mais induzido do plano de Trefftz.
cdtot_ff(av) = av.CDvis + av.CDff
carga_norm(fx, CL) = fx.ccl ./ (4*SREF*CL/(pi*BREF))
eliptica(eta) = sqrt.(max.(0.0, 1 .- eta.^2))

# --------------------------------------------------------------------
# SPLINE INTERPOLADORA MONOTÔNICA (PCHIP)
#
# Hermite cúbica por partes com as inclinações de Fritsch-Carlson. A
# inclinação em cada nó interno é a média harmônica ponderada das
# inclinações dos segmentos vizinhos, zerada quando eles têm sinais
# opostos. É isso que impede a curva de ultrapassar os valores de
# controle e garante que dados monótonos gerem curva monótona.

function inclinacao_extremo(h1, h2, Δ1, Δ2)
    d = ((2h1 + h2)*Δ1 - h1*Δ2)/(h1 + h2)
    d*Δ1 <= 0 && return zero(d)
    (Δ1*Δ2 < 0 && abs(d) > 3abs(Δ1)) && return 3Δ1
    return d
end

function inclinacoes_pchip(xk, yk)
    n = length(xk)
    h = diff(xk)
    Δ = diff(yk) ./ h
    d = zeros(eltype(Δ), n)
    n == 2 && return fill(Δ[1], 2)
    for k in 2:n-1
        if Δ[k-1]*Δ[k] > 0
            w1, w2 = 2h[k] + h[k-1], h[k] + 2h[k-1]
            d[k] = (w1 + w2)/(w1/Δ[k-1] + w2/Δ[k])
        end
    end
    d[1] = inclinacao_extremo(h[1], h[2], Δ[1], Δ[2])
    d[n] = inclinacao_extremo(h[n-1], h[n-2], Δ[n-1], Δ[n-2])
    return d
end

function pchip(xk, yk, xq)
    dk = inclinacoes_pchip(xk, yk)
    map(xq) do q
        i = clamp(searchsortedlast(xk, q), 1, length(xk) - 1)
        h = xk[i+1] - xk[i]
        s = (q - xk[i])/h
        s2 = s*s; s3 = s2*s
        (2s3 - 3s2 + 1)*yk[i] + (s3 - 2s2 + s)*h*dk[i] +
            (-2s3 + 3s2)*yk[i+1] + (s3 - s2)*h*dk[i+1]
    end
end

"Valores de torção nos nós a partir dos decrementos: não crescentes por construção."
valores_no(s) = vcat(zero(eltype(s)), -cumsum(s))

"Torções de estação a partir dos decrementos da spline."
torcoes_spline(nos, s) = pchip(nos, valores_no(s), ETAS)

"Penalidade de torção total, para a asa não pedir mais do que TW_TOTAL."
excesso_total(s) = max(sum(s) - TW_TOTAL, 0.0)^2

"""
Maior variação de torção entre estações vizinhas da asa, em graus.

Monotonicidade e suavidade da spline não impedem um gradiente local
forte: a curva pode ficar plana e virar depressa num trecho curto. Como
essa é a medida que diz se a asa é fabricável e se o método de malha de
vórtices ainda vale, ela é reportada junto com o arrasto.
"""
salto_max(tw) = maximum(abs.(diff(tw)))

"Maior taxa de torção ao longo da envergadura, em graus por metro."
function taxa_max(tw)
    d = abs.(diff(tw)) ./ (diff(ETAS) .* SEMI)
    return maximum(d)
end

# --------------------------------------------------------------------
# MODELO EXATO DO ARRASTO (quadrático nas torções de estação, CL fixo)
#
# Construir cada variante custa 1 + 2 NV + NV(NV-1)/2 rodadas do AVL, por
# isso os modelos ficam em cache e só são refeitos se o aft.avl mudar.

const CACHE = joinpath(SAIDA, "modelo_arrasto_cache.json")

function constroi_arrasto(; so_asa, trim)
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

function cache_valido()
    isfile(CACHE) && mtime(CACHE) > mtime(BASE) || return nothing
    d = JSON.parsefile(CACHE)
    d["nv"] == NV || return nothing
    return d
end

# Só a aeronave completa e compensada. O modelo da asa isolada servia à
# etapa de validação, que já foi feita e está registrada no README: a
# spline monotônica chegou à carga elíptica, com Oswald 1,0035.
function modelo_arrasto()
    d = cache_valido()
    if d !== nothing && haskey(d, "completa")
        println("modelo de arrasto lido do cache ($CACHE)")
        mat(v) = reduce(hcat, [Float64.(r) for r in v])'
        k = d["completa"]
        return (f0 = k["f0"], g = Float64.(k["g"]), H = mat(k["H"]),
                cdvis = k["cdvis"])
    end
    n = 1 + 2NV + NV*(NV-1)÷2
    println("construindo o modelo de arrasto ($n rodadas de AVL)...")
    mc = constroi_arrasto(so_asa = false, trim = true)
    open(CACHE, "w") do io
        JSON.print(io, Dict("nv" => NV, "completa" => Dict(
            "f0" => mc.f0, "g" => mc.g, "cdvis" => mc.cdvis,
            "H" => [mc.H[i, :] for i in 1:size(mc.H, 1)])))
    end
    return mc
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
Folga de estol em relação à estação `eta_lim`: o menor ângulo de estol das
faixas externas a eta_lim menos o das internas. Positiva ou nula significa
que o estol começa para dentro de eta_lim, que é toda a exigência. Ela é
dada em graus só porque é uma diferença de ângulos; não há margem angular
imposta, o limite é zero.
"""
function folga_estol(me, x, eta_lim)
    al = alfa_estol(me, x)
    dentro = me.eta .< eta_lim
    return minimum(al[.!dentro]) - minimum(al[dentro])
end

"Estação onde o estol começa."
eta_critico(me, x) = me.eta[argmin(alfa_estol(me, x))]

"""
Restrição de estol como penalidade direta sobre a folga de estol.

A exigência é uma diferença de mínimos, portanto não é convexa nem suave.
Quando as torções eram livres por estação isso obrigava a reescrevê-la de
forma linear, elegendo uma faixa interna de referência e varrendo a
escolha. Com a spline o problema já é não linear e de poucas variáveis, de
modo que a penalidade direta é ao mesmo tempo mais simples e exata: não há
reformulação a validar, e a viabilidade é conferida no fim pelo próprio
critério.
"""
penalidade_estol(me, eta_lim) = x -> max(DESEMPATE - folga_estol(me, x, eta_lim), 0.0)^2

# --------------------------------------------------------------------
# OTIMIZAÇÃO SOBRE OS DECREMENTOS DA SPLINE
#
# As variáveis naturais são os decrementos s >= 0 entre nós, mas otimizar
# direto na caixa 0 <= s <= DEC_MAX dá problema: o ótimo costuma ter
# vários decrementos EXATAMENTE em zero, ou seja sobre a fronteira, e o
# Fminbox usa barreira logarítmica, que diverge ali. Na prática a busca
# de linha falha em achar ponto finito e o resultado fica preso.
#
# A saída é reparametrizar com uma logística, s = DEC_MAX/(1+exp(-u)):
# qualquer u real dá 0 < s < DEC_MAX, os extremos viram assíntotas em vez
# de paredes, e o problema passa a ser irrestrito. Não há barreira nem
# projeção, e um decremento nulo aparece naturalmente como u bem negativo.

s_de_u(u) = DEC_MAX ./ (1 .+ exp.(-u))
function u_de_s(s)
    f = clamp.(s ./ DEC_MAX, 1e-6, 1 - 1e-6)
    return log.(f ./ (1 .- f))
end

function otimiza(m, nos; pen = nothing, s0 = nothing)
    ns = length(nos) - 1
    s0 === nothing && (s0 = fill(0.3, ns))
    u = u_de_s(clamp.(copy(s0), 0.0, DEC_MAX))
    caminho = [s_de_u(u)]
    for w in (pen === nothing ? [0.0] : [1e2, 1e3, 1e4, 1e5, 1e6])
        function f(v)
            s = s_de_u(v)
            x = livres(torcoes_spline(nos, s))
            return 1e4*cd_mod(m, x) + 1e4*excesso_total(s) +
                   (pen === nothing ? 0.0 : w*pen(x))
        end
        r = optimize(f, u, LBFGS(),
                     Optim.Options(store_trace = true, extended_trace = true,
                                   iterations = 400, g_tol = 1e-10))
        for t in Optim.trace(r)
            haskey(t.metadata, "x") && push!(caminho, s_de_u(t.metadata["x"]))
        end
        u = Optim.minimizer(r)
        push!(caminho, s_de_u(u))
    end
    return s_de_u(u), caminho
end

"""
Partidas com alívio de ponta crescente, dimensionadas pelo número de nós.

Sem a varredura da faixa de referência a diversidade das partidas passa a
ser a única defesa contra mínimo local, então a lista cobre asa sem
torção, washout uniforme em vários tamanhos, washout concentrado na parte
externa e washout concentrado na interna.
"""
function partidas_de(ns, extras = Vector{Float64}[])
    meio = max(1, ns ÷ 2)
    p = [zeros(ns), fill(0.3, ns), fill(1.0, ns), fill(2.0, ns),
         fill(3.0, ns),
         [k <= meio ? 0.0 : 2.0 for k in 1:ns],
         [k <= meio ? 0.2 : 3.0 for k in 1:ns],
         [k <= meio ? 0.0 : 4.0 for k in 1:ns],
         [k <= meio ? 2.0 : 0.2 for k in 1:ns],
         [k*3.0/ns for k in 1:ns],
         [(ns - k + 1)*3.0/ns for k in 1:ns]]
    return vcat(p, [copy(e) for e in extras if length(e) == ns])
end

"Um candidato é viável se o estol começa para dentro de eta_lim."
function viavel(me, nos, s, eta_lim)
    x = livres(torcoes_spline(nos, s))
    return folga_estol(me, x, eta_lim) >= DESEMPATE - 1e-3 &&
           sum(s) <= TW_TOTAL + 1e-6
end

"""
Busca global no espaço dos decrementos.

São poucas variáveis e cada avaliação é apenas uma spline seguida de uma
forma quadrática, sem AVL, de modo que amostrar dezenas de milhares de
pontos custa menos de um segundo. Isso resolve de vez o mínimo local, que
com multipartida apenas já tinha dado resultado incoerente: uma margem de
1,5 grau aparecia inviável enquanto 2,0 graus, que é mais exigente, saía
viável. A amostragem mistura washout uniforme, concentrado na parte
externa e concentrado na interna, para cobrir as formas de interesse.
"""
function amostra_global(m, me, nos, eta_lim; n = 60000, semente = 20240917)
    ns = length(nos) - 1
    rng = MersenneTwister(semente)
    achados = Tuple{Float64,Vector{Float64}}[]
    for k in 1:n
        escala = TW_TOTAL*rand(rng)
        p = rand(rng, ns)
        if k % 3 == 1                      # concentrado na parte externa
            p .*= range(0.15, 1.0; length = ns)
        elseif k % 3 == 2                  # concentrado na parte interna
            p .*= range(1.0, 0.15; length = ns)
        end
        soma = sum(p)
        soma < 1e-9 && continue
        s = clamp.(escala .* p ./ soma, 0.0, DEC_MAX)
        viavel(me, nos, s, eta_lim) || continue
        push!(achados, (cd_mod(m, livres(torcoes_spline(nos, s))), s))
    end
    sort!(achados; by = first)
    return [a[2] for a in achados]
end

"Etapa 3: busca global seguida de polimento local, melhor ótimo viável."
function otimiza_com_estol(m, me, nos, eta_lim; partidas, n_polir = 8)
    pen = penalidade_estol(me, eta_lim)
    brutos = amostra_global(m, me, nos, eta_lim)
    sementes = vcat(brutos[1:min(end, n_polir)], partidas)
    melhor = nothing
    for s0 in sementes
        s, cam = otimiza(m, nos; pen = pen, s0 = s0)
        viavel(me, nos, s, eta_lim) || continue
        x = livres(torcoes_spline(nos, s))
        f = cd_mod(m, x)
        if melhor === nothing || f < melhor.f
            melhor = (s = s, x = x, cam = cam, f = f)
        end
    end
    # se o polimento estragou todos os candidatos, fica com o melhor bruto
    if melhor === nothing
        isempty(brutos) && return nothing
        s = brutos[1]
        x = livres(torcoes_spline(nos, s))
        melhor = (s = s, x = x, cam = [copy(s)], f = cd_mod(m, x))
    end
    return melhor
end

# --------------------------------------------------------------------
# ANIMAÇÃO DO CAMINHO

function anima(caminho, nos, arquivo; so_asa, trim, titulo, n_quadros = 24)
    idx = unique(round.(Int, range(1, length(caminho), length = n_quadros)))
    dados = map(caminho[idx]) do s
        tw = torcoes_spline(nos, s)
        av = avalia(tw; so_asa = so_asa, trim = trim, faixas = true)
        fx = faixas_asa(av.saida)
        (tw = tw, no = valores_no(s), eta = fx.eta,
         carga = carga_norm(fx, av.CL), CDff = av.CDff, e = av.e)
    end
    lo = minimum(minimum(d.tw) for d in dados) - 0.5
    hi = maximum(maximum(d.tw) for d in dados) + 0.5
    cd0 = dados[1].CDff
    fino = range(0, 1; length = 201)
    anim = @animate for (k, d) in enumerate(dados)
        p1 = plot(; xlabel = "η = 2y/b", ylabel = "torção [graus]",
                  title = "torção (spline monotônica)", titlefontsize = 10,
                  titlelocation = :left, legend = false, xlims = (0, 1),
                  ylims = (lo, hi), ESTILO...)
        hline!(p1, [0.0]; color = "#c3c2b7", linewidth = 0.8)
        plot!(p1, fino, pchip(nos, d.no, fino); color = PAL[1], linewidth = 2.6)
        scatter!(p1, nos, d.no; color = PAL[1], markersize = 6,
                 markerstrokecolor = INK, markerstrokewidth = 0.8)
        scatter!(p1, ETAS, d.tw; color = :white, markersize = 3.5,
                 markerstrokecolor = PAL[1], markerstrokewidth = 1.2)
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
base_full = avalia(torcoes(zeros(NV)); trim = true, faixas = true)
@printf("  aeronave compensada  CDff %.6f   Oswald %.4f   profundor %.2f°\n",
        base_full.CDff, base_full.e, base_full.de)
@printf("  CDvis (parasita, constante sob torção) %.6f\n", base_full.CDvis)

m2 = modelo_arrasto()

# A etapa de validação com a asa isolada foi feita e ficou registrada no
# README: a spline monotônica chegou à carga elíptica com Oswald 1,0035 e
# arrasto 0,33 count abaixo do piso plano. Ela não é repetida aqui.

# ====================================================================
# ETAPA 2 -- AERONAVE COMPLETA, IRRESTRITO
# ====================================================================
println("\n", "="^78)
println("ETAPA 2 -- AERONAVE COMPLETA, compensada, sem restrição")
println("="^78)

s2, cam2 = otimiza(m2, NOS_PROJ)
tw2 = torcoes_spline(NOS_PROJ, s2)
x2  = livres(tw2)
ot2 = avalia(tw2; trim = true, faixas = true)
@printf("  decrementos [graus]: %s\n",
        join([@sprintf("%5.2f", v) for v in s2], " "))
@printf("  torções [graus]: %s\n",
        join([@sprintf("%6.2f", v) for v in tw2], " "))
@printf("  CDff   %.6f -> %.6f   (%+.2f counts)\n",
        base_full.CDff, ot2.CDff, 1e4*(ot2.CDff - base_full.CDff))
@printf("  CDtot (CDvis + CDff) %.6f -> %.6f   (%+.2f counts)\n",
        cdtot_ff(base_full), cdtot_ff(ot2),
        1e4*(cdtot_ff(ot2) - cdtot_ff(base_full)))
@printf("  Oswald %.4f   -> %.4f\n", base_full.e, ot2.e)
@printf("  profundor de compensação: %.2f -> %.2f graus\n", base_full.de, ot2.de)

fx2 = faixas_asa(ot2.saida)
dv2 = maximum(abs.(carga_norm(fx2, ot2.CL) .- eliptica(fx2.eta)))
@printf("\n  desvio da carga da asa em relação à elíptica: %.3f\n", dv2)
println("  A aeronave completa não vai para a asa elíptica, e não deveria:")
println("  a empenagem carrega para compensar a arfagem e opera no downwash")
println("  da asa, de modo que o mínimo é do conjunto, não de cada parte.")

# ====================================================================
# ONDE PÔR OS NÓS DA SPLINE
# ====================================================================
# Para arrasto puro a posição dos nós quase não importa: qualquer
# distribuição razoável chega perto do mesmo ótimo, porque o arrasto
# induzido responde a poucos modos de baixa frequência da torção. Quem
# decide a colocação é a RESTRIÇÃO DE ESTOL, que precisa que a curva vire
# depressa junto da raiz do aileron: se a torção vazar para dentro de
# η = 0,56 ela alivia justamente a região que deveria estolar primeiro, e
# paga arrasto sem comprar margem. Por isso a tabela mede as duas coisas,
# e a escolha é feita pela coluna com restrição.

println("\n", "="^78)
println("ONDE PÔR OS NÓS DA SPLINE")
println("="^78)
me = modelo_estol()
const CANDIDATOS_NOS = [
    [0.0, 0.5, 1.0],
    [0.0, 0.33, 0.67, 1.0],
    [0.0, 0.25, 0.5, 0.75, 1.0],
    [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    [0.0, 0.22, 0.40, 0.56, 0.75, 1.0],
    [0.0, 0.48, 0.56, 0.75, 1.0],
    [0.0, 0.30, 0.48, 0.56, 0.70, 0.85, 1.0],
    [0.0, 0.25, 0.45, 0.56, 0.68, 0.84, 1.0],
    [0.0, 0.16, 0.33, 0.5, 0.66, 0.83, 1.0]]

@printf("  %-34s %19s   %19s\n", "", "sem restrição", "estol antes de $(ETA_ESTOL_PROJ)")
@printf("  %-34s %9s %9s   %9s %9s\n", "nós", "counts", "torç tot",
        "counts", "torç tot")
function estuda_nos(m, me, candidatos)
    reg = Tuple{Vector{Float64},Float64,Float64}[]
    for nos in candidatos
        ns = length(nos) - 1
        sk, _ = otimiza(m, nos)
        xk = livres(torcoes_spline(nos, sk))
        c_livre = 1e4*(cd_mod(m, xk) - base_full.CDff)
        sc = otimiza_com_estol(m, me, nos, ETA_ESTOL_PROJ;
                               partidas = partidas_de(ns))
        rot = join([@sprintf("%.2f", v) for v in nos], " ")
        if sc === nothing
            @printf("  %-34s %+9.2f %8.1f°   %19s\n", rot, c_livre, sum(sk),
                    "inviável")
        else
            c_rest = 1e4*(cd_mod(m, sc.x) - base_full.CDff)
            @printf("  %-34s %+9.2f %8.1f°   %+9.2f %8.1f°\n", rot, c_livre,
                    sum(sk), c_rest, sum(sc.s))
            push!(reg, (nos, c_rest, sum(sc.s)))
        end
    end
    return reg
end
reg_nos = estuda_nos(m2, me, CANDIDATOS_NOS)
if !isempty(reg_nos)
    melhor_nos = reg_nos[argmin([r[2] for r in reg_nos])]
    @printf("\n  melhor com restrição: %s, %+.2f counts\n",
            join([@sprintf("%.2f", v) for v in melhor_nos[1]], " "),
            melhor_nos[2])
end
@printf("  adotado no estudo: %s\n",
        join([@sprintf("%.2f", v) for v in NOS_PROJ], " "))

# ====================================================================
# ETAPA 3 -- PROJETO: COM A EXIGÊNCIA DE ESTOL DA FAR 25.203
# ====================================================================
println("\n", "="^78)
println("ETAPA 3 -- PROJETO: estol começando para dentro de η_lim (FAR 25.203)")
println("="^78)

@printf("  aileron: de η = %.2f a 0,90; projeto exige estol antes de η = %.2f\n",
        ETA_AILERON, ETA_ESTOL_PROJ)
@printf("  sem torção:        estol começa em η = %.3f\n", eta_critico(me, zeros(NV)))
@printf("  ótimo da etapa 2:  estol começa em η = %.3f\n", eta_critico(me, x2))
println("  (ambos começam dentro do aileron, o que a FAR 25.203 não admite)\n")

partidas = partidas_de(length(NOS_PROJ) - 1, [s2])

# A varredura vai da posição mais EXIGENTE (η_lim menor) para a mais
# folgada e leva a solução obtida como semente da próxima. Os conjuntos
# viáveis são encaixados (quem estola antes de 0,40 estola antes de 0,50),
# de modo que carregar a solução garante que a tabela não possa sair
# incoerente por falha do otimizador.
const ETAS_LIM = (0.35, 0.40, 0.45, 0.50, 0.56)

function varre_posicoes(m, me, nos, etas_lim, partidas)
    sols = Dict{Float64,Any}()
    anterior = Vector{Float64}[]
    for mg in sort(collect(etas_lim))
        s = otimiza_com_estol(m, me, nos, mg;
                              partidas = vcat(partidas, anterior))
        s === nothing && continue
        av = avalia(torcoes_spline(nos, s.s); trim = true, faixas = true)
        sols[mg] = (s = s, av = av)
        anterior = [copy(s.s)]
    end
    return sols
end

sols = varre_posicoes(m2, me, NOS_PROJ, ETAS_LIM, partidas)

println("  custo de exigir o início do estol para dentro de η_lim:")
@printf("  %-8s %10s %9s %8s %8s %8s %8s %8s\n", "η_lim", "CDff", "counts",
        "Oswald", "estol η", "torç tot", "salto", "°/m")
for mg in ETAS_LIM
    if !haskey(sols, mg)
        @printf("  %-8.2f %10s\n", mg, "inviável"); continue
    end
    r = sols[mg]
    tw = torcoes_spline(NOS_PROJ, r.s.s)
    @printf("  %-8.2f %10.6f %+9.2f %8.4f %8.3f %7.1f° %7.1f° %8.2f\n", mg,
            r.av.CDff, 1e4*(r.av.CDff - base_full.CDff), r.av.e,
            eta_critico(me, r.s.x), sum(r.s.s), salto_max(tw), taxa_max(tw))
end
println("\n  (salto = maior variação de torção entre estações vizinhas da asa;")
println("   °/m = maior taxa de torção ao longo da envergadura. São as medidas")
println("   que dizem se a asa é fabricável: um transporte fica perto de")
println("   0,3 °/m, de modo que valores bem acima disso pedem atenção)")

# Coerência: exigir o estol mais para dentro não pode sair mais barato, e
# uma posição mais folgada não pode ser inviável se uma mais exigente é.
let viaveis = sort(collect(keys(sols)))
    ruim = false
    for i in 1:length(viaveis)-1
        if sols[viaveis[i]].av.CDff < sols[viaveis[i+1]].av.CDff - 1e-6
            @printf("  ATENÇÃO: η_lim %.2f custa menos que %.2f, otimização incoerente\n",
                    viaveis[i], viaveis[i+1])
            ruim = true
        end
    end
    if !isempty(viaveis)
        for mg in ETAS_LIM
            mg > minimum(viaveis) && !haskey(sols, mg) || continue
            @printf("  ATENÇÃO: η_lim %.2f saiu inviável com %.2f viável, incoerente\n",
                    mg, minimum(viaveis))
            ruim = true
        end
    end
    ruim || println("\n  coerência da varredura conferida: custo cresce ao exigir o estol mais para dentro")
end

haskey(sols, ETA_ESTOL_PROJ) ||
    error("a posição de projeto η_lim = $(ETA_ESTOL_PROJ) saiu inviável")
esc = sols[ETA_ESTOL_PROJ]
s3, cam3, ot3 = esc.s.s, esc.s.cam, esc.av
tw3, x3 = torcoes_spline(NOS_PROJ, s3), esc.s.x
@printf("\n  projeto: estol tem de começar para dentro de η = %.2f\n",
        ETA_ESTOL_PROJ)
@printf("  decrementos [graus]: %s\n",
        join([@sprintf("%5.2f", v) for v in s3], " "))
@printf("  torções [graus]: %s\n",
        join([@sprintf("%6.2f", v) for v in tw3], " "))
@printf("  CDff  %.6f  (%+.2f counts contra a asa sem torção)\n",
        ot3.CDff, 1e4*(ot3.CDff - base_full.CDff))
@printf("  CDtot (CDvis + CDff) %.6f  (%+.2f counts)\n", cdtot_ff(ot3),
        1e4*(cdtot_ff(ot3) - cdtot_ff(base_full)))
@printf("  estol começa em η = %.3f, folga em relação a η = %.2f: %+.3f graus\n",
        eta_critico(me, x3), ETA_ESTOL_PROJ, folga_estol(me, x3, ETA_ESTOL_PROJ))
@printf("  torção total raiz-ponta: %.1f graus\n", maximum(tw3) - minimum(tw3))
@printf("  custo da exigência de estol: %+.2f counts sobre a etapa 2\n",
        1e4*(ot3.CDff - ot2.CDff))
@printf("\n  maior salto entre estações vizinhas: %.2f graus, ou %.2f °/m\n",
        salto_max(tw3), taxa_max(tw3))
if taxa_max(tw3) > 1.0
    println("  RESSALVA: essa taxa é bem maior que a de um transporte típico,")
    println("  perto de 0,3 °/m. A solução é monotônica e lisa, mas concentra")
    println("  a mudança de incidência num trecho curto junto à raiz do")
    println("  aileron, que é onde a restrição age. Suavizar mais custa")
    println("  arrasto, e o preço está na tabela de posições acima.")
end

# O caminho da solução acima parte de uma semente da busca global, que já
# cai praticamente em cima do ótimo: animá-lo mostraria só o polimento
# final, convergindo em duas ou três iterações, e nada da otimização. O
# caminho que conta a história é o que parte da asa SEM torção, sem
# semente nenhuma. Ele é refeito aqui e conferido contra o ótimo: tem de
# chegar ao mesmo ponto, e isso é mais um teste de que o ótimo é global.
println("\n  refazendo o caminho a partir da asa sem torção, sem semente...")
s_zero, cam3 = otimiza(m2, NOS_PROJ; pen = penalidade_estol(me, ETA_ESTOL_PROJ),
                       s0 = zeros(length(NOS_PROJ) - 1))
x_zero = livres(torcoes_spline(NOS_PROJ, s_zero))
@printf("  %d pontos no caminho; chega a CDff %.6f contra %.6f do ótimo (%+.2f counts)\n",
        length(cam3), cd_mod(m2, x_zero), cd_mod(m2, x3),
        1e4*(cd_mod(m2, x_zero) - cd_mod(m2, x3)))
@printf("  folga de estol no fim do caminho: %+.3f graus\n",
        folga_estol(me, x_zero, ETA_ESTOL_PROJ))

println("\n  gerando a animação do caminho da etapa 3...")
anima(cam3, NOS_PROJ, "evolucao_3_com_estol.gif"; so_asa = false, trim = true,
      titulo = "etapa 3: aeronave completa com estol restrito")

# --------------------------------------------------------------------
# MAPA DA RESTRIÇÃO NO PLANO DAS VARIÁVEIS LIVRES
#
# No ótimo os decrementos s1, s2, s5 e s6 estão no batente zero, de modo
# que só s3 e s4 ficam livres. Fixando os outros em zero, o problema cabe
# num plano e dá para ver a geometria dele: curvas de nível do arrasto,
# região inviável (estol começando dentro do aileron) e o ótimo. Se o
# ótimo é um vértice definido por duas restrições ativas, isso aparece
# como um bico na fronteira da região viável.
println("\n  mapa do problema no plano (s3, s4)...")
let g3 = range(0, 7; length = 141), g4 = range(0, 5; length = 101)
    plano(a, b) = livres(torcoes_spline(NOS_PROJ, [0.0, 0.0, a, b, 0.0, 0.0]))
    Z = [1e4*cd_mod(m2, plano(a, b)) for b in g4, a in g3]
    Mg = [folga_estol(me, plano(a, b), ETA_ESTOL_PROJ) for b in g4, a in g3]
    q = contour(g3, g4, Z; levels = 25, color = :viridis,
                colorbar_title = "CDff [counts]",
                xlabel = "s3: washout entre η 0,48 e 0,56 [graus]",
                ylabel = "s4: washout entre η 0,56 e 0,70 [graus]",
                title = "arrasto e restrição de estol (s1 = s2 = s5 = s6 = 0)",
                titlefontsize = 10, titlelocation = :left, ESTILO...)
    contourf!(q, g3, g4, Mg .< DESEMPATE; levels = [0.5, 1.5], color = [:gray],
              alpha = 0.18, colorbar_entry = false)
    contour!(q, g3, g4, Mg; levels = [DESEMPATE], color = PAL[2], linewidth = 2.6,
             colorbar_entry = false)
    plot!(q, [c[3] for c in cam3], [c[4] for c in cam3]; color = INK,
          linewidth = 1.2, linestyle = :dot,
          label = "caminho a partir de s = 0 (projeção)")
    scatter!(q, [s3[3]], [s3[4]]; color = PAL[2], marker = :star5, markersize = 11,
             markerstrokecolor = INK, label = "ótimo")
    annotate!(q, 1.2, 4.4, text("região inviável:\nestol começa fora de η_lim", 8,
                                INK2, :left))
    savefig(plot(q; size = (900, 650), dpi = 200),
            joinpath(SAIDA, "mapa_restricao_estol.png"))
end
println("  figura: $SAIDA/mapa_restricao_estol.png")

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

println("\n", "="^78)
println("O QUE O PERFIL DA RAIZ RESOLVERIA")
println("="^78)
al2 = alfa_estol(me, x2)
dentro = me.eta .< ETA_ESTOL_PROJ
alfa_ail = minimum(al2[.!dentro])
@printf("  com a torção da etapa 2, a região externa a η = %.2f estola em α = %.2f°\n",
        ETA_ESTOL_PROJ, alfa_ail)
println("  para o estol começar para dentro dela, o clmax interno teria")
println("  de cair para:")
@printf("  %-10s %12s %12s %10s\n", "η", "clmax atual", "clmax alvo", "redução")
for k in findall(dentro)
    me.eta[k] > 0.45 && continue
    alvo = me.a[k] + me.C[:, k]'x2 + me.b[k]*(alfa_ail - DESEMPATE)
    @printf("  %-10.3f %12.3f %12.3f %9.1f%%\n", me.eta[k], me.lim[k], alvo,
            100*(alvo/me.lim[k] - 1))
end
@printf("\n  ganho potencial: a torção voltaria a ser a da etapa 2 e o arrasto\n")
@printf("  cairia de %.6f para %.6f, ou seja %.1f counts.\n",
        ot3.CDff, ot2.CDff, 1e4*(ot3.CDff - ot2.CDff))

# ====================================================================
# FIGURA DE SÍNTESE
# ====================================================================
etapas = [("2: aeronave completa", s2, cam2, m2, ot2, PAL[2]),
          ("3: com estol restrito", s3, cam3, m2, ot3, PAL[3])]
fino = range(0, 1; length = 201)

# O caminho da etapa 3 passa por cinco estágios de penalidade e é muito
# mais longo que o das outras, então o eixo é a fração do caminho e não a
# iteração bruta: assim as três convergências ficam comparáveis. Os saltos
# da etapa 3 são reais, e não ruído: a cada aumento do peso a solução é
# empurrada para dentro da região viável e paga arrasto por isso.
p1 = plot(; xlabel = "fração do caminho percorrido", ylabel = "CDff em counts",
          title = "(a) caminho da otimização", titlefontsize = 10,
          titlelocation = :left, legend = :topright, ESTILO...)
for (nome, _, cam, m, _, cor) in etapas
    y = [1e4*cd_mod(m, livres(torcoes_spline(NOS_PROJ, s))) for s in cam]
    plot!(p1, range(0, 1; length = length(y)), y; color = cor, linewidth = 2,
          label = nome)
end

p2 = plot(; xlabel = "η = 2y/b", ylabel = "torção [graus]",
          title = "(b) torção ótima de cada etapa", titlefontsize = 10,
          titlelocation = :left, legend = :bottomleft, xlims = (0, 1),
          ESTILO...)
hline!(p2, [0.0]; color = "#c3c2b7", linewidth = 0.8, label = "")
vspan!(p2, [ETA_AILERON, 0.90]; color = "#eceae0", alpha = 0.7, linewidth = 0,
       label = "aileron")
for (nome, s, _, _, _, cor) in etapas
    plot!(p2, fino, pchip(NOS_PROJ, valores_no(s), fino); color = cor,
          linewidth = 2.2, label = nome)
    scatter!(p2, NOS_PROJ, valores_no(s); color = cor, markersize = 5,
             markerstrokecolor = INK, markerstrokewidth = 0.7, label = "")
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

p4 = plot(; xlabel = "o estol tem de começar para dentro de η_lim",
          ylabel = "CDff em counts", legend = :topleft,
          title = "(d) custo da posição exigida para o início do estol", titlefontsize = 10,
          titlelocation = :left, ESTILO...)
mg = sort(collect(keys(sols)))
plot!(p4, mg, [1e4*sols[k].av.CDff for k in mg]; color = PAL[3],
      linewidth = 2.4, marker = :circle, markersize = 5,
      markerstrokecolor = PAL[3], label = "ótimo com restrição")
hline!(p4, [1e4*base_full.CDff]; color = CINZA, linewidth = 2,
       linestyle = :dash, label = "asa sem torção")
hline!(p4, [1e4*ot2.CDff]; color = PAL[2], linewidth = 2, linestyle = :dot,
       label = "ótimo sem restrição (etapa 2)")
scatter!(p4, [ETA_ESTOL_PROJ], [1e4*ot3.CDff]; color = PAL[3], markersize = 10,
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
@printf("%-24s %10.6f %+10.2f %9.4f %9.3f\n", "2: completa irrestrito",
        ot2.CDff, 1e4*(ot2.CDff - base_full.CDff), ot2.e, eta_critico(me, x2))
@printf("%-24s %10.6f %+10.2f %9.4f %9.3f\n", "3: projeto final", ot3.CDff,
        1e4*(ot3.CDff - base_full.CDff), ot3.e, eta_critico(me, x3))
println("\ntodas as torções são splines PCHIP monotônicas de $(length(NOS_PROJ)) nós")

open(joinpath(SAIDA, "torcao_otimizada.json"), "w") do io
    JSON.print(io, Dict(
        "etas" => ETAS,
        "parametrizacao" => Dict(
            "tipo" => "spline PCHIP monotônica sobre decrementos não negativos",
            "nos" => NOS_PROJ,
            "torcao_total_max" => TW_TOTAL),
        "criterio_estol" => Dict(
            "norma" => "FAR 25.203(a): rolamento eficaz até e durante o estol",
            "eta_aileron" => ETA_AILERON,
            "eta_limite_adotado" => ETA_ESTOL_PROJ,
            "desempate_graus" => DESEMPATE),
        "sem_torcao" => Dict("CDff" => base_full.CDff, "e" => base_full.e,
                             "eta_estol" => eta_critico(me, zeros(NV))),
        "geometria" => Dict("versao_designTool" => "v3",
                            "winglet" => true, "Lc_h" => 4.6,
                            "zr_w" => -2.5, "zr_h" => 3.85, "zr_v" => 3.85),
        "etapa_2_completa" => Dict("decrementos" => s2, "torcoes" => tw2,
                                   "CDff" => ot2.CDff, "e" => ot2.e,
                                   "eta_estol" => eta_critico(me, x2)),
        "etapa_3_projeto" => Dict("decrementos" => s3, "torcoes" => tw3,
                                  "CDff" => ot3.CDff,
                                  "CDtot_ff" => cdtot_ff(ot3),
                                  "e" => ot3.e,
                                  "eta_estol" => eta_critico(me, x3),
                                  "folga_estol" => folga_estol(me, x3, ETA_ESTOL_PROJ)),
        "custo_da_posicao" => Dict(string(k) => sols[k].av.CDff for k in mg)), 2)
end
println("\ngravado: $SAIDA/torcao_otimizada.json")
