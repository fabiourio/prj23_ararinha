# Vamos rodar tudo que o LAB 04 exige dentro desse .jl e gerar os gráficos para o relatório.
#
# Rodar de dentro da pasta avl/:   julia lab04.jl
#
# ETAPA 1 -- CONVERGÊNCIA DE MALHA
#
# Em malha de vórtices, corda e envergadura convergem grandezas diferentes:
# o número de painéis na envergadura governa a distribuição de sustentação
# (CL, arrasto induzido), e o número na corda governa a distribuição de
# pressão ao longo dela (Cm, ponto neutro, derivadas de controle). Por isso
# as varreduras são separadas: primeiro a corda com a envergadura fixa,
# depois a envergadura com a corda já escolhida.
#
# O arrasto induzido é julgado pelo CDff do plano de Trefftz, e não pela
# integração de campo próximo (CDind), que oscila com o refino sem
# tendência clara.

using Printf
using Plots

gr()

# --------------------------------------------------------------------
# CONFIGURAÇÃO

const AVL      = "./avl.exe"
const BASE     = "aft.avl"          # geometria de referência do estudo
const VARIANTE = "_malha.avl"       # arquivo temporário de cada nível

const MACH    = 0.85                # ponto de projeto
const CL_PROJ = 0.5053
const MAC     = 6.8995              # corda de referência [m]

# Tolerâncias para considerar a malha convergida (desvio contra a malha
# mais fina da varredura).
const TOL_ALPHA = 0.05              # [graus]
const TOL_CDFF  = 0.01              # relativo (1% ~ 1 count)
const TOL_XNP   = 0.005*MAC         # [m]  (0,5 %MAC)
const TOL_CLA   = 0.01              # relativo
const TOL_CMDE  = 0.01              # relativo

# Estilo das figuras (mesmo dos relatórios anteriores da equipe)
const PAL   = ["#2a78d6", "#eb6834", "#1baf7a"]
const INK2  = "#52514e"
const GRIDC = "#e1e0d9"
const ESTILO = (framestyle = :axes, gridcolor = GRIDC, gridalpha = 1.0,
                gridlinewidth = 0.7, foreground_color_axis = INK2,
                foreground_color_border = "#c3c2b7",
                foreground_color_text = INK2, tickfontsize = 9,
                guidefontsize = 10, legendfontsize = 9,
                background_color = :white)

# --------------------------------------------------------------------
# INFRAESTRUTURA: rodar o AVL e ler números da saída

function roda_avl(cmds::AbstractString)
    tmp = "_cmds.txt"
    write(tmp, cmds)
    saida = try
        read(pipeline(`$AVL`, stdin = tmp), String)
    finally
        rm(tmp, force = true)
    end
    return saida
end

"Último número que casa com o padrão, ou NaN."
function num(txt::AbstractString, pat::Regex)
    m = collect(eachmatch(pat, txt))
    isempty(m) && return NaN
    return parse(Float64, m[end].captures[1])
end

"""
Escreve uma variante do arquivo base com outras malhas.
`malhas` mapeia o nome da superfície para (Nchordwise, Nspanwise).
Superfícies fora do dicionário (nacele) ficam intactas.
"""
function escreve_variante(base, destino, malhas::Dict{String,Tuple{Int,Int}})
    saida = String[]
    sup = ""
    espera_nome = false
    espera_malha = false
    for ln in readlines(base)
        t = strip(ln)
        if t == "SURFACE"
            espera_nome = true
            push!(saida, ln)
        elseif espera_nome && !isempty(t) && !startswith(t, "#")
            sup = t
            espera_nome = false
            espera_malha = haskey(malhas, sup)
            push!(saida, ln)
        elseif espera_malha && !isempty(t) && !startswith(t, "#")
            nc, ns = malhas[sup]
            push!(saida, "$nc 1.0 $ns 1.0")
            espera_malha = false
        else
            push!(saida, ln)
        end
    end
    write(destino, join(saida, "\n") * "\n")
end

"Roda um nível de malha e devolve as métricas de convergência."
function mede(malhas::Dict{String,Tuple{Int,Int}})
    escreve_variante(BASE, VARIANTE, malhas)
    cmds = join(["load $VARIANTE", "oper", "m", "mn $MACH", "",
                 "a c $CL_PROJ", "x", "st", "", "", "quit"], "\n") * "\n"
    s = roda_avl(cmds)
    rm(VARIANTE, force = true)

    # O bloco st termina antes da linha do ponto neutro: a razão de
    # estabilidade espiral logo abaixo também contém "Cnb" e contaminaria
    # a leitura.
    st = split(split(s, "Stability-axis derivatives")[end], "Neutral point")[1]

    return (vortices = num(s, r"(\d+)\s+Vortices"),
            alpha    = num(s, r"Alpha\s*=\s*([-\d.]+)"),
            CL       = num(s, r"CLtot\s*=\s*([-\d.]+)"),
            CDff     = num(s, r"CDff\s*=\s*([-\d.Ee+]+)"),
            CDind    = num(s, r"CDind\s*=\s*([-\d.Ee+]+)"),
            CLa      = num(st, r"CLa\s*=\s*([-\d.]+)"),
            CMa      = num(st, r"Cma\s*=\s*([-\d.]+)"),
            CMde     = num(st, r"Cmd2\s*=\s*([-\d.]+)"),
            xnp      = num(s,  r"Xnp\s*=\s*([-\d.]+)"))
end

# --------------------------------------------------------------------
# RELATÓRIO DE UMA VARREDURA

const METRICAS = [(:alpha, "α [graus]",  :abs, TOL_ALPHA),
                  (:CDff,  "CDff",       :rel, TOL_CDFF),
                  (:xnp,   "xnp [m]",    :abs, TOL_XNP),
                  (:CLa,   "CLα [1/rad]", :rel, TOL_CLA),
                  (:CMde,  "CMδe",       :rel, TOL_CMDE)]

desvio(v, ref, tipo) = tipo === :rel ? abs(v - ref)/abs(ref) : abs(v - ref)

"Imprime a tabela da varredura e devolve o índice do nível recomendado."
function analisa(titulo, rotulos, res)
    fino = res[end]
    println("\n", "="^78)
    println(titulo)
    println("="^78)
    @printf("%-12s %8s %9s %10s %10s %10s %10s\n",
            "malha", "vórtices", "α [°]", "CDff", "xnp [m]", "CLα", "CMδe")
    for (r, x) in zip(rotulos, res)
        @printf("%-12s %8.0f %9.4f %10.6f %10.4f %10.4f %10.5f\n",
                r, x.vortices, x.alpha, x.CDff, x.xnp, x.CLa, x.CMde)
    end

    # Desvio contra a malha mais fina
    println("\ndesvio contra a malha mais fina ($(rotulos[end])):")
    @printf("%-12s %10s %10s %10s %10s %10s   %s\n",
            "malha", "α [°]", "CDff", "xnp [m]", "CLα", "CMδe", "dentro da tol.")
    ok_idx = Int[]
    for (i, (r, x)) in enumerate(zip(rotulos, res))
        ds = [desvio(getfield(x, k), getfield(fino, k), t) for (k, _, t, _) in METRICAS]
        tols = [tol for (_, _, _, tol) in METRICAS]
        ok = all(ds .<= tols)
        ok && push!(ok_idx, i)
        @printf("%-12s %10.4f %10.2e %10.4f %10.2e %10.2e   %s\n",
                r, ds[1], ds[2], ds[3], ds[4], ds[5], ok ? "sim" : "não")
    end

    # Sanidade: as variações entre níveis sucessivos devem estar caindo.
    println("\nvariação entre níveis sucessivos (CDff relativo, xnp [m]):")
    for i in 2:length(res)
        dc = desvio(res[i].CDff, res[i-1].CDff, :rel)
        dx = desvio(res[i].xnp,  res[i-1].xnp,  :abs)
        @printf("  %-12s -> %-12s   CDff %8.2e   xnp %7.4f\n",
                rotulos[i-1], rotulos[i], dc, dx)
    end

    escolhido = isempty(ok_idx) ? length(res) : minimum(ok_idx)
    println("\n>> malha recomendada: $(rotulos[escolhido]) " *
            "($(Int(res[escolhido].vortices)) vórtices)")
    return escolhido
end

function figura(titulo, rotulos, res, arquivo)
    fino = res[end]
    nn = [x.vortices for x in res[1:end-1]]
    p = plot(; xlabel = "número de vórtices",
             ylabel = "desvio da malha mais fina",
             title = titulo, titlefontsize = 11, titlelocation = :left,
             yscale = :log10, legend = :topright, ESTILO...)
    series = [(:CDff, "CDff", :rel, PAL[1], :circle),
              (:xnp, "xnp", :abs, PAL[2], :square),
              (:CLa, "CLα", :rel, PAL[3], :utriangle)]
    for (k, rot, tipo, cor, mk) in series
        d = [max(desvio(getfield(x, k), getfield(fino, k), tipo), 1e-6)
             for x in res[1:end-1]]
        plot!(p, nn, d; color = cor, linewidth = 2, marker = mk,
              markersize = 5, markerstrokecolor = cor, label = rot)
    end
    savefig(plot(p; size = (800, 480), dpi = 200), arquivo)
    println("figura: ", arquivo)
end

# ====================================================================
# VARREDURA A -- CORDA (envergadura fixa)
# ====================================================================
# O refino na corda é o que resolve a distribuição de pressão: Cm, ponto
# neutro e, sobretudo, as derivadas de controle (a charneira do profundor
# precisa cair entre painéis bem resolvidos).

NC = [4, 6, 8, 10, 12, 16]
res_c = [mede(Dict("Wing" => (nc, 24),
                   "Horizontal tail" => (nc, 12),
                   "Vertical tail" => (nc, 8))) for nc in NC]
rot_c = ["$(nc) x corda" for nc in NC]
i_c = analisa("VARREDURA A -- refino na CORDA (envergadura fixa: 24/12/8)",
              rot_c, res_c)
NC_ESCOLHIDO = NC[i_c]
figura("convergência no refino da corda", rot_c, res_c,
       "convergencia_corda.png")

# ====================================================================
# VARREDURA B -- ENVERGADURA (corda já escolhida)
# ====================================================================
# O refino na envergadura resolve a distribuição de sustentação, que é o
# que governa o arrasto induzido e o CLmax pelo método da seção crítica.

NS = [8, 12, 16, 24, 32, 40]
res_s = [mede(Dict("Wing" => (NC_ESCOLHIDO, ns),
                   "Horizontal tail" => (NC_ESCOLHIDO, max(4, ns ÷ 2)),
                   "Vertical tail" => (NC_ESCOLHIDO, max(4, ns ÷ 3))))
         for ns in NS]
rot_s = ["$(ns) x enverg." for ns in NS]
i_s = analisa("VARREDURA B -- refino na ENVERGADURA (corda: $NC_ESCOLHIDO)",
              rot_s, res_s)
NS_ESCOLHIDO = NS[i_s]
figura("convergência no refino da envergadura", rot_s, res_s,
       "convergencia_envergadura.png")

# ====================================================================
# RECOMENDAÇÃO
# ====================================================================
println("\n", "="^78)
println("MALHA RECOMENDADA PARA AS ANÁLISES DO LAB 04")
println("="^78)
@printf("  asa                %2d x %2d\n", NC_ESCOLHIDO, NS_ESCOLHIDO)
@printf("  empenagem horiz.   %2d x %2d\n", NC_ESCOLHIDO, max(4, NS_ESCOLHIDO ÷ 2))
@printf("  empenagem vert.    %2d x %2d\n", NC_ESCOLHIDO, max(4, NS_ESCOLHIDO ÷ 3))
@printf("  total              %d vórtices\n", res_s[i_s].vortices)
println("\nTolerâncias usadas (contra a malha mais fina de cada varredura):")
@printf("  α    %.2f graus\n", TOL_ALPHA)
@printf("  CDff %.1f%%\n", 100*TOL_CDFF)
@printf("  xnp  %.1f%%MAC (%.3f m)\n", 100*TOL_XNP/MAC, TOL_XNP)
@printf("  CLα e CMδe  %.1f%%\n", 100*TOL_CLA)
