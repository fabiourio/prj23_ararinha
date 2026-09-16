# Lab 04 -- figuras dos itens 4 a 7 (secao critica e polares).
#
# Le os CSVs de relatorio_lab04/tables/ gerados pelos scripts Python e
# desenha, no estilo dos relatorios da equipe:
#   - cl_norm x eta na condicao de estol, com a curva limite (item 4);
#   - polar CD x CL dos quatro casos (item 5);
#   - CL x alpha (item 6) e CL x delta_e (item 7).
# O ponto de projeto (CL = 0,5053) e marcado em todas as curvas.
#
# Rodar da raiz do repo:  julia lab04_figuras.jl

using CSV
using DataFrames
using Plots

gr()

PAL = ["#2a78d6", "#eb6834", "#1baf7a"]
INK = "#0b0b0b"
INK2 = "#52514e"
GRIDC = "#e1e0d9"
CL_PROJ = 0.5053

CASOS = [("fwd_livre", "CG diant. sem deflexão", PAL[1], :solid),
         ("fwd_trim", "CG diant. compensado", PAL[1], :dash),
         ("aft_livre", "CG tras. sem deflexão", PAL[2], :solid),
         ("aft_trim", "CG tras. compensado", PAL[2], :dash)]

estilo = (framestyle = :axes, gridcolor = GRIDC, gridalpha = 1.0,
          gridlinewidth = 0.7, foreground_color_axis = INK2,
          foreground_color_border = "#c3c2b7",
          foreground_color_text = INK2, tickfontsize = 10,
          guidefontsize = 11, legendfontsize = 9,
          background_color = :white)

interpola(x, y, x0) = begin
    k = findlast(v -> v <= x0, x)
    k === nothing && return y[1]
    k >= length(x) && return y[end]
    y[k] + (y[k+1] - y[k])*(x0 - x[k])/(x[k+1] - x[k])
end

os_dir = "relatorio_lab04/images"
mkpath("$os_dir/05_secao_critica")
mkpath("$os_dir/06_polares")

# ---------------- item 4: cl x eta na condicao de estol ----------------
# As distribuicoes da asa coincidem entre os casos (mesma asa, mesmo
# alpha de estol), entao cada caso ganha um painel proprio com seus
# alpha_max e CLmax no titulo.
import JSON
critica = JSON.parsefile("avl/saidas/secao_critica.json")

paineis4 = []
for (idx, (caso, rotulo, cor, traco)) in enumerate(CASOS)
    df = CSV.read("relatorio_lab04/tables/clxy_$caso.csv", DataFrame)
    r = critica[caso]
    amax = replace(string(round(r["alpha_max_deg"], digits = 1)),
                   "." => ",")
    clm = replace(string(round(r["CLmax"], digits = 3)), "." => ",")
    p = plot(; title = "$rotulo   (αmax = $(amax)°, CLmax = $clm)",
             titlefontsize = 10, titlelocation = :left,
             legend = (idx == 1 ? :bottomleft : false),
             ylims = (0.0, 2.1), xlims = (0, 1),
             xlabel = (idx > 2 ? "η = 2y/b" : ""),
             ylabel = (idx % 2 == 1 ? "cl local (plano normal)" : ""),
             estilo...)
    plot!(p, df.eta, df.clmax_lim; color = INK, linewidth = 2.0,
          linestyle = :dot, label = "limite clmax (Lab 03)")
    plot!(p, df.eta, df.cl_norm; color = cor, linewidth = 2.2,
          linestyle = traco, label = "distribuição no estol")
    scatter!(p, [r["eta_critica"]],
             [interpola(df.eta, df.clmax_lim, r["eta_critica"])];
             color = INK, markersize = 7, marker = :star5,
             label = "seção crítica")
    push!(paineis4, p)
end
savefig(plot(paineis4...; layout = (2, 2), size = (1300, 820), dpi = 200,
             left_margin = 5Plots.mm, bottom_margin = 5Plots.mm),
        "$os_dir/05_secao_critica/clxy_estol.png")
println("figura: clxy_estol.png")

# ---------------- otimizacao de torcao: antes e depois ----------------
mkpath("$os_dir/06_torcao")
torcao = JSON.parsefile("avl/saidas/torcao.json")
etas_t = [0.0, 0.1011, 0.398, 0.56, 0.90, 1.0]
tors = [torcao["torcoes"][string(e)] for e in etas_t]

pt1 = plot(; xlabel = "η = 2y/b", ylabel = "torção  [graus]",
           title = "distribuição de torção otimizada",
           titlefontsize = 11, titlelocation = :left,
           legend = false, xlims = (0, 1), estilo...)
hline!(pt1, [0.0]; color = "#c3c2b7", linewidth = 0.8)
plot!(pt1, etas_t, tors; color = PAL[1], linewidth = 2.4,
      marker = :circle, markersize = 6)

antes = CSV.read("relatorio_lab04/tables/clxy_fwd_livre_semtorcao.csv",
                 DataFrame)
depois = CSV.read("relatorio_lab04/tables/clxy_fwd_livre.csv", DataFrame)
pt2 = plot(; xlabel = "η = 2y/b", ylabel = "cl local (plano normal)",
           title = "distribuição no estol, antes e depois",
           titlefontsize = 11, titlelocation = :left,
           legend = :bottomleft, ylims = (0, 2.1), xlims = (0, 1),
           estilo...)
plot!(pt2, depois.eta, depois.clmax_lim; color = INK, linewidth = 2.0,
      linestyle = :dot, label = "limite clmax (Lab 03)")
plot!(pt2, antes.eta, antes.cl_norm; color = "#9a99944f",
      linewidth = 2.4, label = "sem torção (estol em η = 0,89)")
plot!(pt2, depois.eta, depois.cl_norm; color = PAL[1], linewidth = 2.4,
      label = "com torção (estol em η = 0,45)")
scatter!(pt2, [0.888], [interpola(antes.eta, antes.clmax_lim, 0.888)];
         color = "#52514e", markersize = 7, marker = :star5,
         label = false)
scatter!(pt2, [0.446], [interpola(depois.eta, depois.clmax_lim, 0.446)];
         color = PAL[1], markersize = 8, marker = :star5,
         markerstrokecolor = INK, label = false)
savefig(plot(pt1, pt2; layout = (1, 2), size = (1300, 480), dpi = 200,
             left_margin = 5Plots.mm, bottom_margin = 6Plots.mm),
        "$os_dir/06_torcao/torcao_antes_depois.png")
println("figura: torcao_antes_depois.png")

# ---------------- itens 5 a 7: polares ----------------
function painel(xcol, ycol, xlab, ylab, pos)
    p = plot(; xlabel = xlab, ylabel = ylab, legend = pos, estilo...)
    for (caso, rotulo, cor, traco) in CASOS
        df = CSV.read("relatorio_lab04/tables/polar_$caso.csv", DataFrame)
        plot!(p, df[!, xcol], df[!, ycol]; color = cor, linewidth = 2.2,
              linestyle = traco, label = rotulo)
        x0 = interpola(df.CL, df[!, xcol], CL_PROJ)
        y0 = interpola(df.CL, df[!, ycol], CL_PROJ)
        scatter!(p, [x0], [y0]; color = cor, markersize = 7,
                 marker = :diamond, markerstrokecolor = INK,
                 label = false)
    end
    p
end

p5 = painel(:CD, :CL, "CD", "CL", :bottomright)
savefig(plot(p5; size = (900, 560), dpi = 200),
        "$os_dir/06_polares/polar_cd_cl.png")
println("figura: polar_cd_cl.png")

p6 = painel(:alpha_deg, :CL, "α  [graus]", "CL", :bottomright)
savefig(plot(p6; size = (900, 560), dpi = 200),
        "$os_dir/06_polares/cl_alpha.png")
println("figura: cl_alpha.png")

p7 = painel(:delta_e_deg, :CL, "δe  [graus]", "CL", :bottomright)
savefig(plot(p7; size = (900, 560), dpi = 200),
        "$os_dir/06_polares/cl_de.png")
println("figura: cl_de.png")
