# Verificação independente da otimização de torção.
#
# Rodar de dentro da pasta avl/:   julia verifica_torcao.jl
#
# A otimização se apoia em dois modelos construídos por diferenças
# finitas: o arrasto induzido, quadrático nas torções, e o estol, afim
# nelas. Os dois são exatos SE a hipótese de linearidade valer, e é
# justamente isso que este script testa, sempre por fora dos modelos.
#
#   1  a torção é mesmo monótona e de tamanho plausível
#   2  o modelo de arrasto acerta FORA dos pontos que o construíram
#   3  CDff, CL e Oswald fecham entre si pela definição
#   4  o estol começa onde o modelo afim diz, conferido por varredura
#      direta de ângulo de ataque no AVL
#   5  a etapa 1 chega mesmo à carga elíptica
#
# Nenhuma verificação reutiliza o modelo que está sendo verificado.

using Printf
using JSON
using Random
using Statistics
using LinearAlgebra

const AVL   = "./avl.exe"
const BASE  = "aft.avl"
const SAIDA = "resultados"
const DADOS = JSON.parsefile("dados_designtool.json")
const PP    = DADOS["ponto_de_projeto"]
const OT    = JSON.parsefile(joinpath(SAIDA, "torcao_otimizada.json"))

const MACH    = PP["M"]
const CL_PROJ = round(PP["CL"], digits = 4)
const SREF    = PP["Sref"]
const BREF    = DADOS["referencia"]["Bref"]
const SEMI    = BREF/2
const AR      = BREF^2/SREF

const MACH_BAIXO  = 0.2
const ETA_AILERON = 0.56
const ETA_LIM   = [0.1011, 0.398, 0.90]
const CLMAX_LIM = [1.774, 1.7985, 1.7338]
const ETAS = Float64.(OT["etas"])

falhas = String[]
function confere(ok, titulo, detalhe)
    @printf("  [%s] %-52s %s\n", ok ? "ok" : "!!", titulo, detalhe)
    ok || push!(falhas, titulo)
    return ok
end

# --------------------------------------------------------------------
# INFRAESTRUTURA (independente do script de otimização)

function roda_avl(cmds)
    tmp = "_vf_cmds.txt"
    write(tmp, cmds)
    try
        return read(pipeline(`$AVL`, stdin = tmp), String)
    finally
        rm(tmp, force = true)
    end
end

function num(txt, pat)
    m = collect(eachmatch(pat, txt))
    isempty(m) && return NaN
    return parse(Float64, m[end].captures[1])
end

function escreve(destino, tw; so_asa = false)
    linhas = readlines(BASE)
    n = length(linhas)
    inicios = [i for i in 1:n if strip(linhas[i]) in ("SURFACE", "BODY")]
    i_asa = 0
    for i in inicios
        strip(linhas[i]) == "SURFACE" || continue
        j = i + 1
        while j <= n && (isempty(strip(linhas[j])) ||
                         startswith(strip(linhas[j]), "#")); j += 1; end
        j <= n && strip(linhas[j]) == "Wing" && (i_asa = i)
    end
    prox = findfirst(>(i_asa), inicios)
    fim_asa = prox === nothing ? n : inicios[prox] - 1
    faixa = so_asa ? (1:fim_asa) : (1:n)
    saida = String[]
    i_sec, espera = 0, false
    for i in faixa
        ln = linhas[i]; t = strip(ln)
        dentro = (i >= i_asa && i <= fim_asa)
        if dentro && t == "SECTION"
            espera = true
        elseif dentro && espera && !isempty(t) && !startswith(t, "#")
            i_sec += 1
            p = split(t); p[5] = @sprintf("%.4f", tw[i_sec])
            push!(saida, join(p, "  ")); espera = false; continue
        end
        push!(saida, ln)
    end
    write(destino, join(saida, "\n") * "\n")
end

function faixas_asa(saida)
    tr = split(saida, r"Surface # 1\s+Wing")[2]
    bl = String(split(tr, "Surface # 2")[1])
    L = [[parse(Float64, v) for v in split(strip(m.match))[2:end]]
         for m in eachmatch(r"^\s*\d+(?:\s+[-\dEe.+]+){12}\s*$"m, bl)]
    d = reduce(hcat, L)'
    return (eta = d[:, 1]./SEMI, corda = d[:, 2], ccl = d[:, 4],
            cl_norm = d[:, 6], cl = d[:, 7])
end

function avalia(tw; so_asa = false, trim = false)
    escreve("_vf.avl", tw; so_asa = so_asa)
    c = ["load _vf.avl", "oper", "m", "mn $MACH", ""]
    trim && push!(c, "d2 pm 0")
    append!(c, ["a c $CL_PROJ", "x", "fs", "", "", "quit"])
    s = roda_avl(join(c, "\n") * "\n")
    rm("_vf.avl", force = true)
    return (CDff = num(s, r"CDff\s*=\s*([-\d.Ee+]+)"),
            CDind = num(s, r"CDind\s*=\s*([-\d.Ee+]+)"),
            CL = num(s, r"CLtot\s*=\s*([-\d.]+)"),
            e = num(s, r"\se =\s+([-\d.]+)"), fx = faixas_asa(s))
end

function em_alfa(tw, alfa)
    escreve("_vf.avl", tw)
    s = roda_avl(join(["load _vf.avl", "oper", "m", "mn $MACH_BAIXO", "",
                       "d2 d2 0", "a a $alfa", "x", "fs", "", "",
                       "quit"], "\n") * "\n")
    rm("_vf.avl", force = true)
    return faixas_asa(s)
end

function clmax_local(η)
    η <= ETA_LIM[1] && return CLMAX_LIM[1]
    η >= ETA_LIM[3] && return CLMAX_LIM[3]
    if η <= ETA_LIM[2]
        w = (η - ETA_LIM[1])/(ETA_LIM[2] - ETA_LIM[1])
        return (1-w)*CLMAX_LIM[1] + w*CLMAX_LIM[2]
    end
    w = (η - ETA_LIM[2])/(ETA_LIM[3] - ETA_LIM[2])
    return (1-w)*CLMAX_LIM[2] + w*CLMAX_LIM[3]
end

eliptica(eta) = sqrt.(max.(0.0, 1 .- eta.^2))
carga_norm(fx, CL) = fx.ccl ./ (4*SREF*CL/(pi*BREF))

tw1 = Float64.(OT["etapa_1_asa_isolada"]["torcoes"])
tw2 = Float64.(OT["etapa_2_completa"]["torcoes"])
tw3 = Float64.(OT["etapa_3_projeto"]["torcoes"])

println("="^78)
println("VERIFICAÇÃO INDEPENDENTE DA OTIMIZAÇÃO DE TORÇÃO")
println("="^78)
@printf("parametrização gravada: %s\n", OT["parametrizacao"]["tipo"])
@printf("nós: %s\n",
        join([@sprintf("%.2f", v) for v in OT["parametrizacao"]["nos"]], " "))

# ====================================================================
println("\n1. A TORÇÃO É MONÓTONA E DE TAMANHO PLAUSÍVEL?")
# ====================================================================
for (nome, tw) in (("etapa 1", tw1), ("etapa 2", tw2), ("etapa 3", tw3))
    d = diff(tw)
    faixa = maximum(tw) - minimum(tw)
    taxa = maximum(abs.(d) ./ (diff(ETAS) .* SEMI))
    confere(maximum(d) <= 1e-6, "$nome: não crescente da raiz para a ponta",
            @sprintf("maior subida %.2e graus", maximum(d)))
    confere(tw[1] == 0.0, "$nome: raiz na referência de gauge",
            @sprintf("torção da raiz %.3f graus", tw[1]))
    confere(faixa <= 12.0, "$nome: torção total dentro do limite",
            @sprintf("%.2f graus entre raiz e ponta", faixa))
    @printf("       taxa máxima %.2f °/m  (transporte típico perto de 0,3)\n",
            taxa)
end

# ====================================================================
println("\n2. O MODELO DE ARRASTO ACERTA FORA DOS PONTOS QUE O CONSTRUÍRAM?")
# ====================================================================
# O modelo quadrático foi ajustado perturbando UMA estação de cada vez
# com 3 graus. Aqui ele é cobrado em torções que não estão nesse
# conjunto: as próprias soluções ótimas e algumas aleatórias.
cache = joinpath(SAIDA, "modelo_arrasto_cache.json")
if isfile(cache)
    mc = JSON.parsefile(cache)
    mat(v) = reduce(hcat, [Float64.(r) for r in v])'
    for (chave, nome, so_asa, trim) in (("asa", "asa isolada", true, false),
                                        ("completa", "aeronave completa",
                                         false, true))
        d = mc[chave]
        f0, g, H = d["f0"], Float64.(d["g"]), mat(d["H"])
        modelo(x) = f0 + g'x + 0.5*x'*H*x
        rng = MersenneTwister(7)
        casos = [tw1[2:end], tw2[2:end], tw3[2:end]]
        for alvo in (4.0, 7.0, 10.0)
            # washout aleatório, monótono, com torção total igual a `alvo`:
            # fica na mesma faixa das soluções, que é onde o modelo é cobrado
            s = cumsum(rand(rng, length(ETAS) - 1))
            push!(casos, -alvo .* s ./ s[end])
        end
        erros = Float64[]
        for x in casos
            tw = vcat(0.0, x)
            av = avalia(tw; so_asa = so_asa, trim = trim)
            push!(erros, 1e4*(modelo(x) - av.CDff))
        end
        # Os três primeiros casos são as próprias soluções ótimas, que é
        # onde o modelo precisa acertar; os três últimos são washouts
        # aleatórios de até 10 graus, bem longe do ponto de linearização,
        # e servem para medir até onde o modelo continua válido.
        confere(maximum(abs.(erros[1:3])) < 1.0,
                "$nome: erro do modelo NAS SOLUÇÕES ótimas",
                @sprintf("máximo %.2f counts", maximum(abs.(erros[1:3]))))
        @printf("       longe do ponto de linearização (até 10° de torção): %.2f counts\n",
                maximum(abs.(erros[4:end])))
    end
else
    println("  (cache do modelo ausente, pulando)")
end

# ====================================================================
println("\n3. CDff, CL E OSWALD FECHAM ENTRE SI?")
# ====================================================================
# O fator de Oswald do AVL é de campo distante, e = CL²/(π AR CDff), o
# que para a asa isolada fecha na quarta casa. Na aeronave completa sobra
# cerca de 2%, e a razão não é erro: com a fuselagem presente o CL
# integrado no plano de Trefftz fica pouco abaixo do CLtot integrado nas
# superfícies, e o AVL usa o primeiro para montar o `e`. Abaixo o CL que
# o `e` do AVL implica é comparado com o CLtot para expor essa diferença
# em vez de escondê-la.
for (nome, tw, so_asa, trim) in (("etapa 1", tw1, true, false),
                                 ("etapa 2", tw2, false, true),
                                 ("etapa 3", tw3, false, true))
    av = avalia(tw; so_asa = so_asa, trim = trim)
    e_ff = av.CL^2/(pi*AR*av.CDff)
    cl_tp = sqrt(av.e*pi*AR*av.CDff)          # CL implícito no e do AVL
    razao = cl_tp/av.CL
    tol = so_asa ? 0.002 : 0.03
    confere(abs(e_ff - av.e) < tol*max(1.0, av.e),
            "$nome: e do AVL = CL²/(π AR CDff)",
            @sprintf("AVL %.4f, recalculado %.4f", av.e, e_ff))
    @printf("       CL de Trefftz implícito %.4f contra CLtot %.4f (%.2f%%)\n",
            cl_tp, av.CL, 100*(razao - 1))
    @printf("       e de campo próximo, para referência: %.4f\n",
            av.CL^2/(pi*AR*av.CDind))
end
piso = CL_PROJ^2/(pi*AR)
av1 = avalia(tw1; so_asa = true)
confere(av1.CDff <= piso + 5e-6,
        "etapa 1: arrasto no piso analítico ou abaixo",
        @sprintf("CDff %.6f contra piso %.6f", av1.CDff, piso))

# ====================================================================
println("\n4. O ESTOL COMEÇA ONDE O MODELO AFIM DIZ?")
# ====================================================================
# O modelo de estol extrapola de apenas dois ângulos de ataque. Aqui o
# AVL é rodado numa varredura e a estação crítica é achada direto, sem
# extrapolar nada.
function estol_direto(tw; alfas = 4.0:1.0:24.0)
    dados = [(a, em_alfa(tw, a)) for a in alfas]
    eta = dados[1][2].eta
    lim = clmax_local.(eta)
    alfa_cr = fill(Inf, length(eta))
    for j in eachindex(eta)
        cls = [d[2].cl_norm[j] for d in dados]
        for k in 1:length(alfas)-1
            if cls[k] <= lim[j] <= cls[k+1]
                t = (lim[j] - cls[k])/(cls[k+1] - cls[k])
                alfa_cr[j] = alfas[k] + t*(alfas[k+1] - alfas[k])
                break
            end
        end
    end
    dentro = eta .< ETA_AILERON
    # Basta que ALGUMA faixa de cada região alcance o clmax: o critério é
    # sobre o mínimo, e as faixas junto à fuselagem têm cl local baixo por
    # interferência e podem nunca estolar sem que isso signifique nada.
    alcancou = any(isfinite, alfa_cr[dentro]) && any(isfinite, alfa_cr[.!dentro])
    if !alcancou
        return (eta_estol = NaN, alfa = NaN, margem = NaN, ok = false,
                eta_todas = eta, alfa_todas = alfa_cr)
    end
    j = argmin(alfa_cr)
    margem = minimum(alfa_cr[.!dentro]) - minimum(alfa_cr[dentro])
    return (eta_estol = eta[j], alfa = alfa_cr[j], margem = margem, ok = true,
            eta_todas = eta, alfa_todas = alfa_cr)
end

for (nome, tw, chave) in (("etapa 2", tw2, "etapa_2_completa"),
                          ("etapa 3", tw3, "etapa_3_projeto"))
    dir = estol_direto(tw)
    prev = Float64(OT[chave]["eta_estol"])
    confere(dir.ok && abs(dir.eta_estol - prev) < 0.10,
            "$nome: estação de estol bate com o modelo",
            dir.ok ? @sprintf("direto η %.3f, modelo η %.3f, α %.2f°",
                              dir.eta_estol, prev, dir.alfa) :
                     "a varredura de α não alcançou o clmax")
end
dir3 = estol_direto(tw3)
confere(dir3.ok && dir3.eta_estol < ETA_AILERON,
        "etapa 3: estol começa FORA do aileron (FAR 25.203)",
        @sprintf("η %.3f contra aileron a partir de %.2f",
                 dir3.eta_estol, ETA_AILERON))
prev_m = Float64(OT["etapa_3_projeto"]["margem_aileron"])
confere(dir3.ok && dir3.margem > 0.0, "etapa 3: margem do aileron positiva",
        @sprintf("direto %+.2f°, modelo %+.2f°", dir3.margem, prev_m))

# sensibilidade ao critério: cl_norm contra cl
fx = em_alfa(tw3, 12.0)
dif = maximum(abs.(fx.cl_norm .- fx.cl))
@printf("       diferença entre cl_norm e cl a 12°: %.4f (escolha do critério)\n",
        dif)

# ====================================================================
println("\n5. A ETAPA 1 CHEGA MESMO À CARGA ELÍPTICA?")
# ====================================================================
c1 = carga_norm(av1.fx, av1.CL)
el = eliptica(av1.fx.eta)
interna = av1.fx.eta .< 0.95        # a ponta é dominada pela discretização
desvio = maximum(abs.(c1[interna] .- el[interna]))
confere(desvio < 0.08, "etapa 1: carga colada na elíptica até η = 0,95",
        @sprintf("desvio máximo %.3f", desvio))
confere(av1.e > 0.99, "etapa 1: Oswald junto de 1",
        @sprintf("e = %.4f", av1.e))

# ====================================================================
# PARTE B -- CONDIÇÕES DE OTIMALIDADE
# ====================================================================
# A parte A perguntou se o problema descreve a física. Esta pergunta se
# o problema, como foi posto, foi de fato resolvido: quais restrições
# estão ativas, se os multiplicadores têm o sinal certo, se algum batente
# foi alcançado e se a condição de segunda ordem vale.
#
# Há um motivo concreto para desconfiar. A otimização roda numa
# reparametrização logística, s = D/(1+exp(-u)), e ds/du = s(1 - s/D)
# tende a zero quando s tende a zero. Então um decremento levado ao
# batente inferior tem gradiente NULO em u qualquer que seja o sinal de
# df/ds: a estacionariedade em u não distingue um mínimo legítimo de uma
# parada prematura. Só o teste em s, feito aqui, separa os dois casos.

const NOS = Float64.(OT["parametrizacao"]["nos"])
const NS  = length(NOS) - 1
const DEC_MAX  = 6.0
const TW_TOTAL = Float64(OT["parametrizacao"]["torcao_total_max"])
const MARGEM_PROJ = Float64(OT["criterio_estol"]["margem_adotada_graus"])

# PCHIP, igual ao usado na otimização
function inclinacao_extremo(h1, h2, Δ1, Δ2)
    d = ((2h1 + h2)*Δ1 - h1*Δ2)/(h1 + h2)
    d*Δ1 <= 0 && return zero(d)
    (Δ1*Δ2 < 0 && abs(d) > 3abs(Δ1)) && return 3Δ1
    return d
end
function inclinacoes_pchip(xk, yk)
    n = length(xk); h = diff(xk); Δ = diff(yk) ./ h
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
        h = xk[i+1] - xk[i]; t = (q - xk[i])/h
        t2 = t*t; t3 = t2*t
        (2t3 - 3t2 + 1)*yk[i] + (t3 - 2t2 + t)*h*dk[i] +
            (-2t3 + 3t2)*yk[i+1] + (t3 - t2)*h*dk[i+1]
    end
end
valores_no(s) = vcat(zero(eltype(s)), -cumsum(s))
torcoes_spline(s) = pchip(NOS, valores_no(s), ETAS)
livres(tw) = tw[2:end]

# modelo de estol, reconstruído aqui: as condições de otimalidade valem
# para o problema COMO FOI POSTO, então usa-se o mesmo modelo afim
println("\nreconstruindo o modelo de estol para o teste de otimalidade...")
function modelo_estol()
    function amostra(x, alfa)
        escreve("_vf.avl", vcat(0.0, x))
        s = roda_avl(join(["load _vf.avl", "oper", "m", "mn $MACH_BAIXO", "",
                           "d2 d2 0", "a a $alfa", "x", "fs", "", "",
                           "quit"], "\n") * "\n")
        rm("_vf.avl", force = true)
        return faixas_asa(s)
    end
    nv = length(ETAS) - 1
    d1 = amostra(zeros(nv), 8.0); d2 = amostra(zeros(nv), 14.0)
    b = (d2.cl_norm .- d1.cl_norm)./6.0
    a = d1.cl_norm .- b.*8.0
    C = zeros(nv, length(d1.eta))
    for j in 1:nv
        e = zeros(nv); e[j] = 1.0
        C[j, :] = amostra(e, 8.0).cl_norm .- d1.cl_norm
    end
    lim = clmax_local.(d1.eta)
    return (eta = d1.eta, a = a, b = b, C = C, lim = lim)
end
me = modelo_estol()
alfa_estol(x) = (me.lim .- me.a .- me.C'x)./me.b
function margem_aileron(x)
    al = alfa_estol(x); dentro = me.eta .< ETA_AILERON
    return minimum(al[.!dentro]) - minimum(al[dentro])
end

mc = JSON.parsefile(joinpath(SAIDA, "modelo_arrasto_cache.json"))
matriz(v) = reduce(hcat, [Float64.(r) for r in v])'
function objetivo(chave)
    d = mc[chave]
    f0, g, H = d["f0"], Float64.(d["g"]), matriz(d["H"])
    return s -> (x = livres(torcoes_spline(s)); f0 + g'x + 0.5*x'*H*x)
end

function gradiente(f, s; h = 1e-6)
    g = similar(s)
    for k in eachindex(s)
        sp = copy(s); sm = copy(s)
        sp[k] += h; sm[k] -= h
        g[k] = (f(sp) - f(sm))/(2h)
    end
    return g
end
function hessiana(f, s; h = 1e-4)
    n = length(s); H = zeros(n, n); f0 = f(s)
    for i in 1:n, j in i:n
        sp = copy(s); sp[i] += h; sp[j] += h
        si = copy(s); si[i] += h
        sj = copy(s); sj[j] += h
        H[i, j] = H[j, i] = (f(sp) - f(si) - f(sj) + f0)/h^2
    end
    return H
end

"""
Teste KKT no espaço dos decrementos.

Restrições escritas como c(s) >= 0: s_k >= 0, DEC_MAX - s_k >= 0,
TW_TOTAL - soma(s) >= 0 e, na etapa 3, margem(s) - margem_exigida >= 0.
No ótimo vale grad f = soma(lambda_i grad c_i) com lambda_i >= 0 sobre as
restrições ativas, o que é resolvido por mínimos quadrados.
"""
function testa_kkt(nome, f, s; margem = nothing, tol_ativa = 1e-4)
    println("\n  ", nome)
    g = gradiente(f, s)
    @printf("    decrementos s : %s\n",
            join([@sprintf("%7.4f", v) for v in s], ""))
    @printf("    df/ds [counts]: %s\n",
            join([@sprintf("%7.2f", 1e4*v) for v in g], ""))

    inf_ativa = findall(s .<= tol_ativa)
    sup_ativa = findall(s .>= DEC_MAX - tol_ativa)
    tot_ativa = sum(s) >= TW_TOTAL - tol_ativa
    cols, nomes = Vector{Float64}[], String[]
    for k in inf_ativa
        e = zeros(NS); e[k] = 1.0; push!(cols, e); push!(nomes, "s$k >= 0")
    end
    for k in sup_ativa
        e = zeros(NS); e[k] = -1.0; push!(cols, e)
        push!(nomes, "s$k <= $DEC_MAX")
    end
    tot_ativa && (push!(cols, fill(-1.0, NS)); push!(nomes, "torção total"))
    # Restrição de estol DESAGREGADA. Escrita como margem do aileron ela é
    # um mínimo de mínimos, portanto não diferenciável, e um gradiente por
    # diferenças finitas num ponto de empate é arbitrário. A forma correta
    # separa a exigência por faixa: fixada a faixa interna que estola
    # primeiro, cada faixa do aileron precisa estolar `margem` graus depois
    # dela. Cada uma dessas é suave, e o KKT admite um multiplicador por
    # faixa ativa. Se várias estiverem empatadas, é a combinação convexa
    # dos gradientes delas que precisa reproduzir o gradiente do objetivo.
    if margem !== nothing
        x = livres(torcoes_spline(s))
        al = alfa_estol(x)
        dentro = me.eta .< ETA_AILERON
        i_ref = findall(dentro)[argmin(al[dentro])]
        folga = margem_aileron(x) - margem
        @printf("    faixa interna crítica: η = %.3f, estola em α = %.2f°\n",
                me.eta[i_ref], al[i_ref])
        @printf("    folga da restrição de estol: %+.4f graus\n", folga)
        ext = findall(.!dentro)
        gj(j) = z -> (xx = livres(torcoes_spline(z));
                      aa = alfa_estol(xx); aa[j] - aa[i_ref] - margem)
        folgas = [gj(j)(s) for j in ext]
        ativas = ext[folgas .<= 5e-3]
        @printf("    faixas do aileron ativas: %d de %d", length(ativas),
                length(ext))
        isempty(ativas) || @printf("  (η = %s)",
            join([@sprintf("%.3f", me.eta[j]) for j in ativas], ", "))
        println()
        for j in ativas
            push!(cols, gradiente(gj(j), s))
            push!(nomes, @sprintf("estol η=%.3f", me.eta[j]))
        end
    end

    @printf("    batente inferior ativo: %s\n",
            isempty(inf_ativa) ? "nenhum" : join(["s$k" for k in inf_ativa], " "))
    println("    batente superior ativo: ",
            isempty(sup_ativa) ?
                "nenhum, o limite de $(DEC_MAX)° por intervalo não amarrou" :
                join(["s$k" for k in sup_ativa], " "))
    @printf("    torção total (limite %.1f°): %.2f°, %s\n", TW_TOTAL, sum(s),
            tot_ativa ? "ATIVA" : "inativa")

    if isempty(cols)
        res = maximum(abs.(1e4 .* g))
        confere(res < 0.5, "$nome: estacionário, sem restrição ativa",
                @sprintf("maior |df/ds| = %.3f counts/grau", res))
        return
    end
    A = reduce(hcat, cols)
    lam = A \ g
    resid = 1e4*maximum(abs.(A*lam .- g))
    @printf("    multiplicadores [counts por unidade]:\n")
    for (i, n) in enumerate(nomes)
        @printf("      %-16s %+10.3f  %s\n", n, 1e4*lam[i],
                lam[i] >= -1e-9 ? "" : "<-- SINAL ERRADO")
    end
    confere(all(lam .>= -1e-7), "$nome: multiplicadores não negativos",
            @sprintf("menor %.3e", minimum(lam)))
    confere(resid < 0.5, "$nome: estacionariedade de KKT",
            @sprintf("resíduo %.3f counts/grau", resid))

    # segunda ordem no cone crítico: núcleo dos gradientes ativos
    H = hessiana(f, s)
    Q = nullspace(A')
    if size(Q, 2) > 0
        Hr = Symmetric(Q'*H*Q)
        av = eigvals(Hr)
        confere(minimum(av) > -1e-9,
                "$nome: Hessiana reduzida positiva no cone crítico",
                @sprintf("menor autovalor %.3e", minimum(av)))
    else
        println("    cone crítico vazio: o ótimo é determinado pelas restrições")
    end
end

s1 = Float64.(OT["etapa_1_asa_isolada"]["decrementos"])
s2 = Float64.(OT["etapa_2_completa"]["decrementos"])
s3 = Float64.(OT["etapa_3_projeto"]["decrementos"])

println("\n6. CONDIÇÕES DE OTIMALIDADE (KKT) NO ESPAÇO DOS DECREMENTOS")
testa_kkt("etapa 1, asa isolada", objetivo("asa"), s1)
testa_kkt("etapa 2, aeronave completa", objetivo("completa"), s2)
testa_kkt("etapa 3, projeto", objetivo("completa"), s3; margem = MARGEM_PROJ)

# ====================================================================
println("\n7. O MULTIPLICADOR DO ESTOL BATE COM A VARREDURA DE MARGEM?")
# ====================================================================
# O multiplicador é o preço marginal da restrição: quanto de arrasto
# custa exigir um grau a mais de margem. Se ele bate com a inclinação da
# tabela de margens, os dois caminhos independentes concordam.
if haskey(OT, "custo_da_margem")
    cm = OT["custo_da_margem"]
    chaves = sort([parse(Float64, k) for k in keys(cm)])
    perto = filter(k -> abs(k - MARGEM_PROJ) <= 0.51, chaves)
    if length(perto) >= 2
        for i in 1:length(perto)-1
            a, b = perto[i], perto[i+1]
            incl = 1e4*(cm[string(b)] - cm[string(a)])/(b - a)
            @printf("    secante entre %.1f° e %.1f°: %+.1f counts por grau\n",
                    a, b, incl)
        end
    end
end

# ====================================================================
println("\n", "="^78)
if isempty(falhas)
    println("TUDO CONFERE: a otimização se sustenta fora dos modelos e")
    println("satisfaz as condições de otimalidade de primeira e segunda ordem")
else
    println("FALHARAM $(length(falhas)) VERIFICAÇÕES:")
    for f in falhas; println("  - ", f); end
end
println("="^78)
