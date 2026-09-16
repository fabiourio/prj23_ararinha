# Lab 04 -- figura do estudo de convergencia de malha do AVL.
#
# Le relatorio_lab04/tables/convergencia_malha.csv (gerado por
# lab04_convergencia.py) e desenha o desvio percentual de cada metrica em
# relacao a malha mais fina, em tres paineis com a mesma escala log.
#
# Rodar da raiz do repo:  julia lab04_convergencia_fig.jl

using CSV
using DataFrames
using Plots
using Measures

gr()

df = CSV.read("relatorio_lab04/tables/convergencia_malha.csv", DataFrame)
fino = df[end, :]
i_adot = findfirst(==("12x40"), df.malha_asa)

# Paleta e tintas do estilo dos relatorios (estilo.py do Lab 03)
PAL = ["#2a78d6", "#eb6834", "#1baf7a"]
INK2 = "#52514e"
GRIDC = "#e1e0d9"

NN = df.vortices[1:end-1]
x_adot = df.vortices[i_adot]

erro(col) = max.(100 .* abs.(df[1:end-1, col] .- fino[col]) ./
                 abs(fino[col]), 1e-3)

paineis = [
    (:alpha_deg, "α para CL = 0,5053", PAL[1]),
    (:CDff, "CD,ind (plano de Trefftz)", PAL[2]),
    (:xnp_m, "ponto neutro", PAL[3]),
]

graficos = []
for (idx, (col, rotulo, cor)) in enumerate(paineis)
    p = plot(NN, erro(col);
             yscale = :log10,
             xlims = (-150, 4400),
             ylims = (0.002, 4.0),
             yticks = ([0.01, 0.1, 1.0], ["0,01", "0,1", "1"]),
             xticks = ([0, 1000, 2000, 3000, 4000],
                       ["0", "1000", "2000", "3000", "4000"]),
             marker = :circle, markersize = 6,
             markerstrokecolor = cor, linewidth = 2.5, color = cor,
             title = rotulo, titlefontsize = 12, titlelocation = :left,
             tickfontsize = 10, guidefontsize = 11,
             legend = false, framestyle = :axes,
             gridcolor = GRIDC, gridalpha = 1.0, gridlinewidth = 0.7,
             foreground_color_axis = INK2,
             foreground_color_border = "#c3c2b7",
             foreground_color_text = INK2,
             xlabel = "número de vórtices")
    vline!(p, [x_adot]; color = INK2, linewidth = 1.0, linestyle = :dash)
    if idx == 1
        plot!(p; ylabel = "desvio da malha 24×80  [%]")
        annotate!(p, x_adot + 150, 2.6,
                  text("malha adotada", 9, INK2, :left))
    end
    push!(graficos, p)
end

fig = plot(graficos...; layout = (1, 3), size = (1500, 460), dpi = 200,
           left_margin = 6mm, bottom_margin = 8mm, top_margin = 3mm,
           background_color = :white)

caminho = "relatorio_lab04/images/03_convergencia/convergencia_malha.png"
savefig(fig, caminho)
println("figura: ", caminho)
