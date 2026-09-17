# Vamos rodar tudo que o LAB 04 exige dentro desse .jl e gerar os gráficos
# para o relatório. A malha já vem pronta e verificada (ver
# convergencia_malha.jl) e os dados do designTool vêm de
# dados_designtool.json (ver dados_designtool.py).
#
# Rodar de dentro da pasta avl/:   julia lab04.jl
#
# ETAPA 1 -- PONTO DE PROJETO
# Estabelece a condição de referência e levanta a distribuição de
# sustentação da asa sem torção, que é a linha de base da otimização de
# torção.

using Printf
using Plots
using JSON

gr()

# --------------------------------------------------------------------
# CONFIGURAÇÃO

const AVL = "./avl.exe"
const DADOS = JSON.parsefile("dados_designtool.json")
const PP = DADOS["ponto_de_projeto"]

const MACH    = PP["M"]
const CL_PROJ = round(PP["CL"], digits = 4)
const SREF    = PP["Sref"]
const BREF    = DADOS["referencia"]["Bref"]
const CREF    = DADOS["referencia"]["Cref"]
const SEMI    = BREF/2

const PAL   = ["#2a78d6", "#eb6834", "#1baf7a"]
const INK   = "#0b0b0b"
const INK2  = "#52514e"
const ESTILO = (framestyle = :axes, gridcolor = "#e1e0d9", gridalpha = 1.0,
                gridlinewidth = 0.7, foreground_color_axis = INK2,
                foreground_color_border = "#c3c2b7",
                foreground_color_text = INK2, tickfontsize = 9,
                guidefontsize = 10, legendfontsize = 9,
                background_color = :white)

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

"""
Roda um caso no ponto de projeto. `restricao` é o comando de arfagem
("d2 d2 0" para profundor fixo, "d2 pm 0" para arfagem compensada) e `it`
é a incidência da empenagem.
"""
function caso(arquivo; restricao = "d2 d2 0", it = 0.0, faixas = false)
    cmds = ["load $arquivo", "oper", "m", "mn $MACH", "",
            "de", "1 $it", "", restricao, "a c $CL_PROJ", "x"]
    faixas && append!(cmds, ["fs", ""])
    append!(cmds, ["st", "", "", "quit"])
    s = roda_avl(join(cmds, "\n") * "\n")
    st = split(split(s, "Stability-axis derivatives")[end], "Neutral point")[1]
    return (saida = s,
            alpha = num(s, r"Alpha\s*=\s*([-\d.]+)"),
            CL    = num(s, r"CLtot\s*=\s*([-\d.]+)"),
            CD    = num(s, r"CDtot\s*=\s*([-\d.]+)"),
            CDff  = num(s, r"CDff\s*=\s*([-\d.Ee+]+)"),
            CDind = num(s, r"CDind\s*=\s*([-\d.Ee+]+)"),
            CDvis = num(s, r"CDvis\s*=\s*([-\d.Ee+]+)"),
            e     = num(s, r"\se =\s+([-\d.]+)"),
            Cm    = num(s, r"Cmtot\s*=\s*([-\d.]+)"),
            de    = num(s, r"elevator\s*=\s*([-\d.]+)"),
            CLa   = num(st, r"CLa\s*=\s*([-\d.]+)"),
            CMa   = num(st, r"Cma\s*=\s*([-\d.]+)"),
            xnp   = num(s, r"Xnp\s*=\s*([-\d.]+)"))
end

"""
Distribuição por faixas da asa direita: η, corda, c·cl e cl local.
As colunas do bloco `fs` do AVL são, depois do índice da faixa:
Yle, Chord, Area, c·cl, ai, cl_norm, cl, cd, cdv, cm_c/4, cm_LE, C.P.x/c.
A saída do AVL usa terminador de linha do Windows, então a classe final
precisa aceitar o retorno de carro (por isso \\s e não apenas espaço).
"""
function faixas_asa(saida)
    trecho = split(saida, r"Surface # 1\s+Wing")[2]
    bloco = String(split(trecho, "Surface # 2")[1])
    linhas = [[parse(Float64, v) for v in split(strip(m.match))[2:end]]
              for m in eachmatch(r"^\s*\d+(?:\s+[-\dEe.+]+){12}\s*$"m, bloco)]
    isempty(linhas) && error("não achei as faixas da asa na saída do AVL")
    d = reduce(hcat, linhas)'
    return (eta = d[:, 1]./SEMI, corda = d[:, 2], ccl = d[:, 4], cl = d[:, 7])
end

# ====================================================================
# ETAPA 1 -- PONTO DE PROJETO
# ====================================================================

println("="^78)
println("ETAPA 1 -- PONTO DE PROJETO")
println("="^78)
@printf("  W0  = %10.1f N   (%.1f kgf)\n", PP["W0_N"], PP["W0_N"]/9.81)
@printf("  W   = %10.1f N   (%.1f kgf, %.0f%% de combustível, 100%% de carga)\n",
        PP["W_N"], PP["W_kgf"], 100*PP["fuel_frac"])
@printf("  h   = %10.1f m\n", PP["h_m"])
@printf("  rho = %10.5f kg/m3      a = %.2f m/s\n", PP["rho"], PP["a_inf"])
@printf("  M   = %10.2f            V = %.1f m/s\n", PP["M"], PP["V"])
@printf("  CL  = %10.4f            Sref = %.3f m2\n", PP["CL"], SREF)

# --------------------------------------------------------------------
# Estado da aeronave sem torção nos dois CGs

println("\n", "-"^78)
println("ESTADO NO PONTO DE PROJETO (asa sem torção, profundor neutro)")
println("-"^78)
@printf("%-16s %8s %9s %9s %9s %8s %9s %9s\n",
        "arquivo", "α [°]", "CD", "CDff", "CDind", "e", "Cm", "xnp [m]")
casos = Dict{String,Any}()
for f in ("fwd.avl", "aft.avl")
    r = caso(f; faixas = true)
    casos[f] = r
    @printf("%-16s %8.3f %9.5f %9.6f %9.6f %8.4f %9.5f %9.4f\n",
            f, r.alpha, r.CD, r.CDff, r.CDind, r.e, r.Cm, r.xnp)
end

ref_dt = DADOS["polar_designtool"]
@printf("\n  comparação com o designTool no mesmo ponto:\n")
@printf("    CD0 (parasita)   AVL usa %.5f (CDp do cabeçalho) | designTool %.5f\n",
        casos["aft.avl"].CDvis, ref_dt["CD0"])
@printf("    CD induzido      AVL %.5f (Trefftz)              | designTool %.5f\n",
        casos["aft.avl"].CDff, ref_dt["CDind"])
@printf("    CD total         AVL %.5f                        | designTool %.5f  (%+.1f%%)\n",
        casos["aft.avl"].CD, ref_dt["CD"],
        100*(casos["aft.avl"].CD/ref_dt["CD"] - 1))
@printf("    fator de Oswald  AVL %.4f                         | designTool %.4f\n",
        casos["aft.avl"].e, ref_dt["e"])
# A diferença de arrasto induzido tem uma causa identificada: o designTool
# credita ao winglet um alongamento efetivo maior que o geométrico, e o
# modelo do AVL não tem winglet.
@printf("\n    o designTool usa alongamento efetivo %.3f contra %.3f geométrico\n",
        ref_dt["AR_eff"], ref_dt["AR_geom"])
@printf("    (winglet, %+.0f%%), que o modelo do AVL não tem. Conferindo:\n",
        100*(ref_dt["AR_eff"]/ref_dt["AR_geom"] - 1))
@printf("      CL²/(π AR_eff e) = %.5f  reproduz o CDind do designTool\n",
        CL_PROJ^2/(pi*ref_dt["AR_eff"]*ref_dt["e"]))
@printf("      com o AR geométrico daria %.5f, próximo do que o AVL calcula\n",
        CL_PROJ^2/(pi*ref_dt["AR_geom"]*ref_dt["e"]))

# --------------------------------------------------------------------
# Distribuição de sustentação contra a elíptica
#
# Para sustentação total fixa, a distribuição que minimiza o arrasto
# induzido de uma asa isolada é a elíptica. Comparar a distribuição atual
# com ela mostra onde a asa está sobrecarregada, que é exatamente a
# informação que a otimização de torção vai usar.
#
# Carga por envergadura: L'(y) = q c(y) cl(y). Na elíptica,
# L'(η) = L'(0) sqrt(1 - η²) com L'(0) = 4 q Sref CL / (π b), de modo que
# normalizando c·cl por 4 Sref CL / (π b) a elíptica vira sqrt(1 - η²).

r = casos["aft.avl"]
fx = faixas_asa(r.saida)
norma = 4*SREF*r.CL/(pi*BREF)
carga = fx.ccl ./ norma
eliptica = sqrt.(max.(0.0, 1 .- fx.eta.^2))

excesso = carga .- eliptica
i_max = argmax(excesso)
println("\n", "-"^78)
println("DISTRIBUIÇÃO DE SUSTENTAÇÃO CONTRA A ELÍPTICA (linha de base)")
println("-"^78)
@printf("  maior excesso sobre a elíptica: %+.3f em η = %.3f\n",
        excesso[i_max], fx.eta[i_max])
@printf("  carga na ponta (η > 0,8): %+.3f em média sobre a elíptica\n",
        sum(excesso[fx.eta .> 0.8])/count(fx.eta .> 0.8))
@printf("  fator de Oswald atual: e = %.4f  (elíptica daria e = 1)\n", r.e)
@printf("  arrasto induzido atual: %.6f\n", r.CDff)
@printf("  mínimo teórico (e = 1): %.6f  ->  margem de %.1f counts\n",
        r.CL^2/(pi*BREF^2/SREF), 1e4*(r.CDff - r.CL^2/(pi*BREF^2/SREF)))
println("""
  A margem acima é um LIMITE SUPERIOR otimista, não uma meta. Ela supõe
  carga elíptica na asa isolada, e a aeronave real não chega lá por três
  motivos: a empenagem carrega para compensar a arfagem, e o ótimo do
  conjunto asa mais empenagem não é a asa elíptica; a fuselagem e as
  naceles perturbam o carregamento; e a restrição de estol vai impedir
  parte do alívio de ponta. O valor serve para dimensionar a expectativa
  da otimização de torção, que atua justamente sobre o excesso de carga
  na ponta mostrado acima.""")

p = plot(; xlabel = "η = 2y/b", ylabel = "carga normalizada  c·cl / (4 S CL / π b)",
         title = "distribuição de sustentação da asa sem torção",
         titlefontsize = 11, titlelocation = :left, legend = :bottomleft,
         xlims = (0, 1), ESTILO...)
plot!(p, fx.eta, eliptica; color = INK, linewidth = 2, linestyle = :dash,
      label = "elíptica (arrasto induzido mínimo)")
plot!(p, fx.eta, carga; color = PAL[1], linewidth = 2.4, marker = :circle,
      markersize = 3, markerstrokecolor = PAL[1], label = "asa sem torção")
plot!(p, fx.eta, excesso; color = PAL[2], linewidth = 2,
      label = "excesso sobre a elíptica")
hline!(p, [0.0]; color = "#c3c2b7", linewidth = 0.8, label = "")
savefig(plot(p; size = (860, 520), dpi = 200), "resultados/ponto_projeto_carga.png")
println("\nfigura: ponto_projeto_carga.png")

# --------------------------------------------------------------------
# Guarda o estado de referência para as etapas seguintes
open("resultados/ponto_de_projeto.txt", "w") do io
    println(io, "# estado no ponto de projeto, asa sem torção")
    println(io, "M $MACH")
    println(io, "CL $CL_PROJ")
    for f in ("fwd.avl", "aft.avl")
        c = casos[f]
        println(io, "$f alpha $(c.alpha) CD $(c.CD) CDff $(c.CDff) ",
                "e $(c.e) Cm $(c.Cm) xnp $(c.xnp)")
    end
end
println("gravado: ponto_de_projeto.txt")
