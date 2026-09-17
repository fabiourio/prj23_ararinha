# Convergência de malha do modelo AVL -- estudo próprio, separado das
# análises do Lab 04. Define a malha que o lab04.jl usa depois.
#
# Rodar de dentro da pasta avl/:   julia convergencia_malha.jl
#
# MÉTODO
#
# 1. A grandeza de interesse é o ARRASTO, que é o que a otimização da asa
#    vai minimizar. O estudo é conduzido sobre ele; os demais coeficientes
#    são verificados ao final, para garantir que a malha escolhida também
#    os resolve.
#
# 2. O arrasto induzido é lido no plano de Trefftz (CDff) e não pela
#    integração de campo próximo (CDind), que oscila com o refino sem
#    tendência clara.
#
# 3. Corda e envergadura resolvem coisas diferentes. O arrasto induzido
#    depende da distribuição de sustentação ao longo da ENVERGADURA, então
#    é o refino em envergadura que o governa. O refino em CORDA resolve a
#    distribuição de pressão ao longo dela, que importa para momento,
#    ponto neutro e derivadas de controle. As varreduras são separadas.
#
# 4. O critério não é comparar com a malha mais fina. Além de um certo
#    refino o resultado sai do platô e volta a oscilar, porque a alocação
#    de painéis entre as seções da asa muda de forma descontínua. A
#    tolerância útil nunca é menor que esse ruído: adota-se a menor malha
#    a partir da qual o resultado permanece dentro da banda de ruído do
#    platô.

using Printf
using Plots
using Statistics

gr()

# --------------------------------------------------------------------
# CONFIGURAÇÃO

const AVL      = "./avl.exe"
const BASE     = "aft.avl"
const VARIANTE = "_malha.avl"

const MACH    = 0.85
const CL_PROJ = 0.5053
const MAC     = 6.8995

const N_PLATO = 4       # níveis finais usados para estimar a banda de ruído
const TOL_CD  = 0.005   # tolerância no CD (0,5% ~ 0,5 count)

# Tolerâncias dos demais coeficientes, usadas para escolher o refino na
# CORDA (que o arrasto não define, por ser insensível a ele).
const TOL_ALPHA = 0.05          # [graus]
const TOL_XNP   = 0.005*6.8995  # [m], 0,5 %MAC
const TOL_CLA   = 0.01          # relativo
const TOL_CMDE  = 0.02          # relativo

const PAL   = ["#2a78d6", "#eb6834", "#1baf7a"]
const INK   = "#0b0b0b"
const INK2  = "#52514e"
const ESTILO = (framestyle = :axes, gridcolor = "#e1e0d9", gridalpha = 1.0,
                gridlinewidth = 0.7, foreground_color_axis = INK2,
                foreground_color_border = "#c3c2b7",
                foreground_color_text = INK2, tickfontsize = 9,
                guidefontsize = 10, legendfontsize = 9,
                background_color = :white)

const W, H, V = "Wing", "Horizontal tail", "Vertical tail"

# --------------------------------------------------------------------
# INFRAESTRUTURA

function roda_avl(cmds::AbstractString)
    tmp = "_cmds.txt"
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

"Variante do arquivo base com outras malhas. A nacele fica intacta."
function escreve_variante(malhas::Dict{String,Tuple{Int,Int}})
    saida, sup = String[], ""
    espera_nome = espera_malha = false
    for ln in readlines(BASE)
        t = strip(ln)
        if t == "SURFACE"
            espera_nome = true; push!(saida, ln)
        elseif espera_nome && !isempty(t) && !startswith(t, "#")
            sup, espera_nome = t, false
            espera_malha = haskey(malhas, sup); push!(saida, ln)
        elseif espera_malha && !isempty(t) && !startswith(t, "#")
            nc, ns = malhas[sup]
            push!(saida, "$nc 1.0 $ns 1.0"); espera_malha = false
        else
            push!(saida, ln)
        end
    end
    write(VARIANTE, join(saida, "\n") * "\n")
end

function mede(malhas::Dict{String,Tuple{Int,Int}})
    escreve_variante(malhas)
    s = roda_avl(join(["load $VARIANTE", "oper", "m", "mn $MACH", "",
                       "a c $CL_PROJ", "x", "st", "", "", "quit"], "\n") * "\n")
    rm(VARIANTE, force = true)
    # o bloco st termina antes do ponto neutro: a razão de estabilidade
    # espiral logo abaixo também contém "Cnb" e contaminaria a leitura
    st = split(split(s, "Stability-axis derivatives")[end], "Neutral point")[1]
    return (vortices = num(s, r"(\d+)\s+Vortices"),
            alpha = num(s, r"Alpha\s*=\s*([-\d.]+)"),
            CDff  = num(s, r"CDff\s*=\s*([-\d.Ee+]+)"),
            CDind = num(s, r"CDind\s*=\s*([-\d.Ee+]+)"),
            e     = num(s, r"\se =\s+([-\d.]+)"),
            CLa   = num(st, r"CLa\s*=\s*([-\d.]+)"),
            CMde  = num(st, r"Cmd2\s*=\s*([-\d.]+)"),
            xnp   = num(s, r"Xnp\s*=\s*([-\d.]+)"))
end

malha_asa(nc, ns) = Dict(W => (nc, ns),
                         H => (nc, max(4, ns ÷ 2)),
                         V => (nc, max(3, ns ÷ 3)))

# --------------------------------------------------------------------
# PLATÔ E BANDA DE RUÍDO

"""
Estima o platô pelos N_PLATO níveis mais finos: a referência é a média
deles e a banda de ruído é o maior desvio dentro desse grupo. A tolerância
efetiva nunca é menor que essa banda, porque exigir mais que o ruído do
próprio método não faz sentido.
"""
function plato(vals)
    fim = vals[end-N_PLATO+1:end]
    ref = mean(fim)
    banda = maximum(abs.(fim .- ref))/abs(ref)
    tol = max(TOL_CD, banda)
    i_rec = length(vals)
    for i in 1:length(vals)
        if all(abs(vals[j] - ref)/abs(ref) <= tol for j in i:length(vals))
            i_rec = i; break
        end
    end
    return ref, banda, tol, i_rec
end

function tabela(titulo, rotulos, res)
    vals = [x.CDff for x in res]
    ref, banda, tol, i_rec = plato(vals)
    println("\n", "="^84)
    println(titulo)
    println("="^84)
    @printf("%-12s %8s %11s %10s %10s %11s %8s\n",
            "malha", "vórtices", "CDff", "desvio", "d.sucess.", "CDind", "e")
    for (i, (r, x)) in enumerate(zip(rotulos, res))
        d  = 100*(x.CDff - ref)/ref
        ds = i == 1 ? "" : @sprintf("%9.3f%%", 100*abs(x.CDff - res[i-1].CDff)/res[i-1].CDff)
        marca = i == i_rec ? "  <== platô" : ""
        @printf("%-12s %8.0f %11.6f %9.2f%% %10s %11.6f %8.4f%s\n",
                r, x.vortices, x.CDff, d, ds, x.CDind, x.e, marca)
    end
    @printf("\nplatô: CDff = %.6f   banda de ruído = %.2f%%   tolerância efetiva = %.2f%%\n",
            ref, 100*banda, 100*tol)
    println(">> malha recomendada: $(rotulos[i_rec]) ($(Int(res[i_rec].vortices)) vórtices)")
    return ref, banda, i_rec
end

"""
Figura em dois painéis:

  (a) CDff contra o refino, em escala linear, com a banda de tolerância
      sombreada em torno do valor de platô. Mostra a subida, a entrada na
      banda e a permanência nela.

  (b) Prova de suficiência: para cada malha, a MAIOR variação que qualquer
      refino posterior ainda provoca. Quando essa curva cruza a tolerância
      e não volta a subir, está demonstrado que refinar mais não muda o
      resultado além do tolerado. Esse gráfico é preferível à extrapolação
      de Richardson neste caso porque a alocação de painéis entre as
      seções da asa é discreta: o erro não cai de forma suave com 1/N e
      não há faixa assintótica limpa para extrapolar.
"""
function figura(titulo, xlab, xs, res, ref, banda, i_rec, arquivo)
    cd = [x.CDff for x in res]
    tol = max(TOL_CD, banda)
    n = length(xs)

    pa = plot(; xlabel = xlab, ylabel = "CDff (plano de Trefftz)",
              title = "(a) " * titulo, titlefontsize = 10,
              titlelocation = :left, legend = :bottomright, ESTILO...)
    plot!(pa, [xs[1], xs[end]], fill(ref*(1 + tol), 2);
          fillrange = fill(ref*(1 - tol), 2), fillalpha = 0.16,
          fillcolor = PAL[3], linealpha = 0,
          label = @sprintf("tolerância de %.1f%%", 100*tol))
    hline!(pa, [ref]; color = PAL[3], linewidth = 1.5, linestyle = :dash,
           label = @sprintf("platô  CDff = %.5f", ref))
    plot!(pa, xs, cd; color = PAL[1], linewidth = 2.4, marker = :circle,
          markersize = 5, markerstrokecolor = PAL[1], label = "CDff")
    scatter!(pa, [xs[i_rec]], [cd[i_rec]]; color = PAL[2], markersize = 10,
             marker = :star5, markerstrokecolor = INK, label = "malha adotada")

    # ---- painel (b): maior variação causada por qualquer refino posterior
    resto = [maximum(abs(cd[j] - cd[i])/abs(cd[i]) for j in i:n) for i in 1:n-1]
    pb = plot(; xlabel = xlab,
              ylabel = "maior variação de qualquer refino posterior",
              title = "(b) suficiência: o que ainda muda se refinar mais",
              titlefontsize = 10, titlelocation = :left,
              legend = :topright, ESTILO...)
    hline!(pb, [tol]; color = PAL[3], linewidth = 1.5, linestyle = :dash,
           label = @sprintf("tolerância de %.1f%%", 100*tol))
    plot!(pb, xs[1:n-1], resto; color = PAL[1], linewidth = 2.4,
          marker = :circle, markersize = 5, markerstrokecolor = PAL[1],
          label = "variação remanescente")
    if i_rec <= n-1
        scatter!(pb, [xs[i_rec]], [resto[i_rec]]; color = PAL[2],
                 markersize = 10, marker = :star5, markerstrokecolor = INK,
                 label = @sprintf("adotada: %.2f%% remanescente",
                                  100*resto[i_rec]))
    end
    ylims!(pb, 0, max(tol*3, minimum([maximum(resto), tol*6])))

    savefig(plot(pa, pb; layout = (1, 2), size = (1240, 470), dpi = 200,
                 left_margin = 6Plots.mm, bottom_margin = 6Plots.mm), arquivo)
    @printf("  variação remanescente na malha adotada: %.2f%%\n",
            100*resto[min(i_rec, n-1)])
    println("figura: ", arquivo)
end

# ====================================================================
# ETAPA 1 -- ENVERGADURA: é ela que governa o arrasto induzido
# ====================================================================
# A asa tem seis seções (cinco segmentos), então o AVL precisa de pelo
# menos um punhado de painéis para distribuir entre elas: níveis muito
# grosseiros simplesmente não montam e são descartados.

NS_TENTA = [8, 10, 12, 14, 16, 20, 24, 28, 32, 40]
res_bruto = [(ns, mede(malha_asa(8, ns))) for ns in NS_TENTA]
validos = [(ns, r) for (ns, r) in res_bruto if isfinite(r.CDff)]
descartados = [ns for (ns, r) in res_bruto if !isfinite(r.CDff)]
isempty(descartados) ||
    println("\nníveis descartados (o AVL não montou a malha): ",
            join(descartados, ", "))

NS = [ns for (ns, _) in validos]
res_s = [r for (_, r) in validos]
rot_s = ["$ns painéis" for ns in NS]
ref_s, banda_s, i_s = tabela(
    "ETAPA 1 -- refino na ENVERGADURA (corda fixa em 8)", rot_s, res_s)
NS_ESC = NS[i_s]
figura("arrasto contra o refino na envergadura", "painéis na semi-envergadura",
       NS, res_s, ref_s, banda_s, i_s, "resultados/conv_envergadura.png")

# ====================================================================
# ETAPA 2 -- CORDA: o arrasto é insensível a ela
# ====================================================================
NC_TENTA = [2, 3, 4, 6, 8, 12, 16, 24]
res_bruto_c = [(nc, mede(malha_asa(nc, NS_ESC))) for nc in NC_TENTA]
validos_c = [(nc, r) for (nc, r) in res_bruto_c if isfinite(r.CDff)]
NC = [nc for (nc, _) in validos_c]
res_c = [r for (_, r) in validos_c]
rot_c = ["$nc painéis" for nc in NC]
ref_c, banda_c, i_c = tabela(
    "ETAPA 2 -- refino na CORDA (envergadura fixa em $NS_ESC)", rot_c, res_c)
figura("arrasto contra o refino na corda", "painéis na corda",
       NC, res_c, ref_c, banda_c, i_c, "resultados/conv_corda.png")

uteis = [x.CDff for (nc, x) in zip(NC, res_c) if nc >= 4]
faixa_cd = 100*(maximum(uteis) - minimum(uteis))/ref_c
@printf("\nO arrasto varia apenas %.2f%% em toda a faixa de refino na corda,\n",
        faixa_cd)
println("de modo que NÃO é o arrasto que define esse refino. A corda é")
println("escolhida na etapa seguinte, pelos coeficientes que dependem dela.")

# ====================================================================
# ETAPA 3 -- CORDA pelos coeficientes de momento e controle
# ====================================================================
# O AVL encaixa a charneira do profundor na borda de painel mais próxima.
# Com espaçamento cosseno essa borda muda de lugar conforme o refino, de
# modo que parte da variação do CMδe não é falta de resolução, e sim a
# charneira efetiva andando. O bom Nchordwise é o que põe a borda perto do
# Xhinge pedido, e não simplesmente o maior.
charneira(nc, xh = 0.70) = (e = [(1 - cos(pi*i/nc))/2 for i in 0:nc];
                            e[argmin(abs.(e .- xh))])

println("\n", "="^84)
println("ETAPA 3 -- refino na CORDA pelos demais coeficientes")
println("="^84)
println("charneira efetiva do profundor (Xhinge pedido = 0,70):")
for nc in (4, 6, 8, 10, 12, 16, 20, 24)
    bom = abs(charneira(nc) - 0.70) < 0.02 ? "  <== bem posicionada" : ""
    @printf("  Nc = %2d  ->  x/c = %.4f  (%+.1f%% da corda)%s\n",
            nc, charneira(nc), 100*(charneira(nc) - 0.70), bom)
end

fina = mede(Dict(W => (24, 32), H => (24, 16), V => (16, 10)))
@printf("\nreferência fina: %d vórtices | α %.4f°  CDff %.6f  xnp %.4f  CLα %.4f  CMδe %.5f\n\n",
        fina.vortices, fina.alpha, fina.CDff, fina.xnp, fina.CLa, fina.CMde)

# Só refinos que põem a charneira perto de 0,70 entram como candidatos.
ns_h, ns_v = max(4, NS_ESC ÷ 2), max(3, NS_ESC ÷ 3)
candidatos = [("asa 8 / EH 8",  Dict(W => (8, NS_ESC), H => (8, ns_h),  V => (8, ns_v))),
              ("asa 8 / EH 16", Dict(W => (8, NS_ESC), H => (16, ns_h), V => (8, ns_v))),
              ("asa 8 / EH 24", Dict(W => (8, NS_ESC), H => (24, ns_h), V => (8, ns_v))),
              ("asa 16 / EH 16", Dict(W => (16, NS_ESC), H => (16, ns_h), V => (8, ns_v)))]

@printf("%-16s %8s | %8s %8s %9s %8s %8s   %s\n", "candidato", "vórtices",
        "d.α [°]", "d.CDff", "d.xnp [m]", "d.CLα", "d.CMδe", "atende")
println("-"^88)
aprovados = []
for (nome, m) in candidatos
    r = mede(m)
    d = (al = abs(r.alpha - fina.alpha),
         cd = abs(r.CDff - fina.CDff)/fina.CDff,
         xn = abs(r.xnp - fina.xnp),
         cl = abs(r.CLa - fina.CLa)/abs(fina.CLa),
         cm = abs(r.CMde - fina.CMde)/abs(fina.CMde))
    ok = d.al <= TOL_ALPHA && d.cd <= TOL_CD && d.xn <= TOL_XNP &&
         d.cl <= TOL_CLA && d.cm <= TOL_CMDE
    @printf("%-16s %8.0f | %8.3f %7.2f%% %9.3f %7.2f%% %7.2f%%   %s\n",
            nome, r.vortices, d.al, 100*d.cd, d.xn, 100*d.cl, 100*d.cm,
            ok ? "sim" : "não")
    ok && push!(aprovados, (nome, m, r))
end

if isempty(aprovados)
    error("nenhum candidato atendeu todas as tolerâncias; reveja as faixas")
end
nome, m, r = aprovados[argmin([x[3].vortices for x in aprovados])]

# ====================================================================
# MALHA ADOTADA
# ====================================================================
println("\n", "="^84)
println("MALHA ADOTADA -- $nome")
println("="^84)
for s in (W, H, V)
    @printf("  %-18s %2d (corda) x %2d (envergadura)\n", s, m[s][1], m[s][2])
end
@printf("  %-18s %d vórtices\n", "total", r.vortices)
@printf("  charneira efetiva do profundor: x/c = %.4f\n", charneira(m[H][1]))
println("\ncritérios:")
@printf("  envergadura pelo ARRASTO      -> %d painéis (platô, banda de ruído %.2f%%)\n",
        NS_ESC, 100*banda_s)
@printf("  corda pelos MOMENTOS/CONTROLE -> asa %d, EH %d\n", m[W][1], m[H][1])

open("resultados/malha_adotada.txt", "w") do io
    for s in (W, H, V)
        println(io, "$s $(m[s][1]) $(m[s][2])")
    end
    println(io, "vortices $(Int(r.vortices))")
end
println("\ngravado: malha_adotada.txt")
