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
#   durante o estol. A tradução geométrica disso não é uma estação
#   arbitrária: é que o estol NÃO comece na região do aileron, que nesta
#   asa vai de η = 0,56 a η = 0,90. A margem é dada em ÂNGULO DE ATAQUE:
#   quantos graus a região do aileron ainda aguenta depois que a região
#   interna estola.
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
            push!(ys, parse(Float64, split(t)[2])); espera = false
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
const ETA_AILERON = 0.56          # raiz do aileron: o estol tem de vir antes
const MARGEM_PROJ = 0.5           # [graus] de margem adotada no projeto
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

function modelos_arrasto()
    d = cache_valido()
    if d !== nothing
        println("modelos de arrasto lidos do cache ($CACHE)")
        mat(v) = reduce(hcat, [Float64.(r) for r in v])'
        le(k) = (f0 = d[k]["f0"], g = Float64.(d[k]["g"]),
                 H = mat(d[k]["H"]), cdvis = d[k]["cdvis"])
        return le("asa"), le("completa")
    end
    n = 1 + 2NV + NV*(NV-1)÷2
    println("construindo os modelos de arrasto ($(2n) rodadas de AVL)...")
    ma = constroi_arrasto(so_asa = true,  trim = false)
    mc = constroi_arrasto(so_asa = false, trim = true)
    open(CACHE, "w") do io
        emp(m) = Dict("f0" => m.f0, "g" => m.g, "cdvis" => m.cdvis,
                      "H" => [m.H[i, :] for i in 1:size(m.H, 1)])
        JSON.print(io, Dict("nv" => NV, "asa" => emp(ma), "completa" => emp(mc)))
    end
    return ma, mc
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
Restrição de estol como penalidade direta sobre a margem do aileron.

A exigência é uma diferença de mínimos, portanto não é convexa nem suave.
Quando as torções eram livres por estação isso obrigava a reescrevê-la de
forma linear, elegendo uma faixa interna de referência e varrendo a
escolha. Com a spline o problema já é não linear e de poucas variáveis, de
modo que a penalidade direta é ao mesmo tempo mais simples e exata: não há
reformulação a validar, e a viabilidade é conferida no fim pelo próprio
critério.
"""
penalidade_estol(me, margem) = x -> max(margem - margem_aileron(me, x), 0.0)^2

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

"Um candidato é viável se o estol começa fora do aileron com a margem pedida."
function viavel(me, nos, s, margem)
    x = livres(torcoes_spline(nos, s))
    return eta_critico(me, x) < ETA_AILERON &&
           margem_aileron(me, x) >= margem - 0.05 &&
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
function amostra_global(m, me, nos, margem; n = 60000, semente = 20240917)
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
        viavel(me, nos, s, margem) || continue
        push!(achados, (cd_mod(m, livres(torcoes_spline(nos, s))), s))
    end
    sort!(achados; by = first)
    return [a[2] for a in achados]
end

"Etapa 3: busca global seguida de polimento local, melhor ótimo viável."
function otimiza_com_estol(m, me, nos, margem; partidas, n_polir = 8)
    pen = penalidade_estol(me, margem)
    brutos = amostra_global(m, me, nos, margem)
    sementes = vcat(brutos[1:min(end, n_polir)], partidas)
    melhor = nothing
    for s0 in sementes
        s, cam = otimiza(m, nos; pen = pen, s0 = s0)
        viavel(me, nos, s, margem) || continue
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
base_asa  = avalia(torcoes(zeros(NV)); so_asa = true, faixas = true)
base_full = avalia(torcoes(zeros(NV)); trim = true, faixas = true)
piso = CL_PROJ^2/(pi*AR)
@printf("  asa isolada          CDff %.6f   Oswald %.4f\n",
        base_asa.CDff, base_asa.e)
@printf("  aeronave trimada     CDff %.6f   Oswald %.4f   profundor %.2f°\n",
        base_full.CDff, base_full.e, base_full.de)
@printf("  CDvis (parasita, constante sob torção) %.6f\n", base_full.CDvis)
@printf("  piso analítico CL²/(π AR) = %.6f  (carga elíptica)\n", piso)

m1, m2 = modelos_arrasto()

# ====================================================================
# ETAPA 1 -- VALIDAÇÃO: ASA ISOLADA, IRRESTRITO
# ====================================================================
println("\n", "="^78)
println("ETAPA 1 -- VALIDAÇÃO: asa isolada, sem restrição")
println("="^78)
println("  A resposta é conhecida: carga elíptica, Oswald 1, CDff no piso.")
println("  A torção é uma spline monotônica de $(length(NOS_PROJ)) nós.")

s1, cam1 = otimiza(m1, NOS_PROJ)
tw1 = torcoes_spline(NOS_PROJ, s1)
x1  = livres(tw1)
ot1 = avalia(tw1; so_asa = true, faixas = true)

# Referência: o mesmo problema com as torções livres por estação, que é o
# menor arrasto que o modelo admite. A diferença mede o que a exigência de
# monotonicidade e suavidade custa nesta etapa.
function otimo_livre(m)
    lo, hi = fill(TW_MIN, NV), fill(TW_MAX, NV)
    r = optimize(z -> 1e4*cd_mod(m, z), lo, hi, zeros(NV), Fminbox(LBFGS()),
                 Optim.Options(iterations = 400, g_tol = 1e-12))
    return Optim.minimizer(r)
end
xl1 = otimo_livre(m1)

@printf("\n  decrementos [graus]: %s\n",
        join([@sprintf("%5.2f", v) for v in s1], " "))
@printf("  torções [graus]: %s\n",
        join([@sprintf("%6.2f", v) for v in tw1], " "))
@printf("  CDff   %.6f -> %.6f   (%+.2f counts)\n",
        base_asa.CDff, ot1.CDff, 1e4*(ot1.CDff - base_asa.CDff))
@printf("  Oswald %.4f   -> %.4f\n", base_asa.e, ot1.e)
@printf("  distância ao piso analítico: %+.2f counts\n", 1e4*(ot1.CDff - piso))
@printf("  o modelo previu %.6f e o AVL deu %.6f (erro %.2f counts)\n",
        cd_mod(m1, x1), ot1.CDff, 1e4*abs(cd_mod(m1, x1) - ot1.CDff))
@printf("  torção livre por estação daria %.6f, ou seja a spline monotônica\n",
        cd_mod(m1, xl1))
@printf("  custa %+.2f counts e entrega uma asa construível\n",
        1e4*(cd_mod(m1, x1) - cd_mod(m1, xl1)))

# A asa tem 6 graus de diedro, de modo que NÃO é plana. O piso CL²/(π AR)
# vale para asa plana; para uma asa não plana o mínimo de Munk fica um
# pouco abaixo dele, e o fator de Oswald pode passar de 1.
ok1 = ot1.e > 0.99 && ot1.CDff <= piso + 5e-6
if ok1
    println("\n  >> VALIDADO: Oswald ~ 1 e arrasto no piso analítico.")
    ot1.e > 1.0 && @printf("     (Oswald %.4f passa de 1 porque a asa tem diedro:\n      com asa não plana o mínimo de Munk fica abaixo do elíptico plano)\n",
                           ot1.e)
else
    @printf("\n  >> Oswald %.4f, %+.2f counts do piso plano.\n",
            ot1.e, 1e4*(ot1.CDff - piso))
    println("     A spline monotônica não alcança exatamente a elíptica porque")
    println("     a carga elíptica desta asa exigiria torção não monótona perto")
    println("     da raiz. O desvio mede esse preço.")
end
println("\n  gerando a animação da convergência à elíptica...")
anima(cam1, NOS_PROJ, "evolucao_1_asa_isolada.gif"; so_asa = true, trim = false,
      titulo = "etapa 1: asa isolada, sem restrição")

# ====================================================================
# ETAPA 2 -- AERONAVE COMPLETA, IRRESTRITO
# ====================================================================
println("\n", "="^78)
println("ETAPA 2 -- AERONAVE COMPLETA, trimada, sem restrição")
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

@printf("  %-34s %19s   %19s\n", "", "sem restrição", "com estol $(MARGEM_PROJ)°")
@printf("  %-34s %9s %9s   %9s %9s\n", "nós", "counts", "torç tot",
        "counts", "torç tot")
function estuda_nos(m, me, candidatos)
    reg = Tuple{Vector{Float64},Float64,Float64}[]
    for nos in candidatos
        ns = length(nos) - 1
        sk, _ = otimiza(m, nos)
        xk = livres(torcoes_spline(nos, sk))
        c_livre = 1e4*(cd_mod(m, xk) - base_full.CDff)
        sc = otimiza_com_estol(m, me, nos, MARGEM_PROJ;
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
println("ETAPA 3 -- PROJETO: estol antes do aileron (FAR 25.203)")
println("="^78)

@printf("  aileron: de η = %.2f a 0,90\n", ETA_AILERON)
@printf("  sem torção:        estol em η = %.3f, margem do aileron %+.2f°\n",
        eta_critico(me, zeros(NV)), margem_aileron(me, zeros(NV)))
@printf("  ótimo da etapa 2:  estol em η = %.3f, margem do aileron %+.2f°\n",
        eta_critico(me, x2), margem_aileron(me, x2))
println("  (margem negativa significa que o estol começa DENTRO do aileron,")
println("   o que a FAR 25.203 não admite: o rolamento perde eficácia)\n")

partidas = partidas_de(length(NOS_PROJ) - 1, [s2])

# A varredura vai da margem MAIOR para a menor e leva a solução obtida
# como semente da próxima. Os conjuntos viáveis são encaixados (o que
# atende 2 graus atende 1,5), de modo que carregar a solução garante que
# a tabela não possa sair incoerente por falha do otimizador.
const MARGENS = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0)

function varre_margens(m, me, nos, margens, partidas)
    sols = Dict{Float64,Any}()
    anterior = Vector{Float64}[]
    for mg in reverse(margens)
        s = otimiza_com_estol(m, me, nos, mg;
                              partidas = vcat(partidas, anterior))
        s === nothing && continue
        av = avalia(torcoes_spline(nos, s.s); trim = true, faixas = true)
        sols[mg] = (s = s, av = av)
        anterior = [copy(s.s)]
    end
    return sols
end

sols = varre_margens(m2, me, NOS_PROJ, MARGENS, partidas)

println("  custo da margem de estol exigida:")
@printf("  %-8s %10s %9s %8s %8s %8s %8s %8s\n", "margem", "CDff", "counts",
        "Oswald", "estol η", "torç tot", "salto", "°/m")
for mg in MARGENS
    if !haskey(sols, mg)
        @printf("  %-8.1f %10s\n", mg, "inviável"); continue
    end
    r = sols[mg]
    tw = torcoes_spline(NOS_PROJ, r.s.s)
    @printf("  %-8.1f %10.6f %+9.2f %8.4f %8.3f %7.1f° %7.1f° %8.2f\n", mg,
            r.av.CDff, 1e4*(r.av.CDff - base_full.CDff), r.av.e,
            eta_critico(me, r.s.x), sum(r.s.s), salto_max(tw), taxa_max(tw))
end
println("\n  (salto = maior variação de torção entre estações vizinhas da asa;")
println("   °/m = maior taxa de torção ao longo da envergadura. São as medidas")
println("   que dizem se a asa é fabricável: um transporte fica perto de")
println("   0,3 °/m, de modo que valores bem acima disso pedem atenção)")

# Coerência: exigir mais margem não pode sair mais barato, e uma margem
# menor não pode ser inviável se uma maior é viável.
let viaveis = sort(collect(keys(sols)))
    ruim = false
    for i in 1:length(viaveis)-1
        if sols[viaveis[i]].av.CDff > sols[viaveis[i+1]].av.CDff + 1e-6
            @printf("  ATENÇÃO: margem %.1f custa mais que %.1f, otimização incoerente\n",
                    viaveis[i], viaveis[i+1])
            ruim = true
        end
    end
    if !isempty(viaveis)
        for mg in MARGENS
            mg < maximum(viaveis) && !haskey(sols, mg) || continue
            @printf("  ATENÇÃO: margem %.1f saiu inviável com %.1f viável, incoerente\n",
                    mg, maximum(viaveis))
            ruim = true
        end
    end
    ruim || println("\n  coerência da varredura conferida: custo cresce com a margem")
end

haskey(sols, MARGEM_PROJ) ||
    error("a margem de projeto $(MARGEM_PROJ)° saiu inviável")
esc = sols[MARGEM_PROJ]
s3, cam3, ot3 = esc.s.s, esc.s.cam, esc.av
tw3, x3 = torcoes_spline(NOS_PROJ, s3), esc.s.x
@printf("\n  margem adotada no projeto: %.1f grau de ângulo de ataque\n",
        MARGEM_PROJ)
@printf("  decrementos [graus]: %s\n",
        join([@sprintf("%5.2f", v) for v in s3], " "))
@printf("  torções [graus]: %s\n",
        join([@sprintf("%6.2f", v) for v in tw3], " "))
@printf("  CDff  %.6f  (%+.2f counts contra a asa sem torção)\n",
        ot3.CDff, 1e4*(ot3.CDff - base_full.CDff))
@printf("  CDtot (CDvis + CDff) %.6f  (%+.2f counts)\n", cdtot_ff(ot3),
        1e4*(cdtot_ff(ot3) - cdtot_ff(base_full)))
@printf("  estol em η = %.3f com margem de %+.2f graus\n",
        eta_critico(me, x3), margem_aileron(me, x3))
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
    println("  arrasto, e o preço está na tabela de margens acima.")
end

println("\n  gerando a animação do caminho da etapa 3...")
anima(cam3, NOS_PROJ, "evolucao_3_com_estol.gif"; so_asa = false, trim = true,
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
etapas = [("1: asa isolada", s1, cam1, m1, ot1, PAL[1]),
          ("2: aeronave completa", s2, cam2, m2, ot2, PAL[2]),
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
            "margem_adotada_graus" => MARGEM_PROJ),
        "sem_torcao" => Dict("CDff" => base_full.CDff, "e" => base_full.e,
                             "eta_estol" => eta_critico(me, zeros(NV))),
        "etapa_1_asa_isolada" => Dict("decrementos" => s1, "torcoes" => tw1,
                                      "CDff" => ot1.CDff, "e" => ot1.e,
                                      "piso_analitico" => piso),
        "etapa_2_completa" => Dict("decrementos" => s2, "torcoes" => tw2,
                                   "CDff" => ot2.CDff, "e" => ot2.e,
                                   "eta_estol" => eta_critico(me, x2)),
        "etapa_3_projeto" => Dict("decrementos" => s3, "torcoes" => tw3,
                                  "CDff" => ot3.CDff,
                                  "CDtot_ff" => cdtot_ff(ot3),
                                  "e" => ot3.e,
                                  "eta_estol" => eta_critico(me, x3),
                                  "margem_aileron" => margem_aileron(me, x3)),
        "custo_da_margem" => Dict(string(k) => sols[k].av.CDff for k in mg)), 2)
end
println("\ngravado: $SAIDA/torcao_otimizada.json")
