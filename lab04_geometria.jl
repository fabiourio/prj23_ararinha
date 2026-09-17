# Lab 04 -- vistas da geometria da aeronave a partir dos arquivos do AVL.
#
# Le avl/fwd.avl (secoes de cada superficie, com Xle Yle Zle Chord Ainc) e
# desenha, no estilo dos relatorios da equipe:
#   - planta da aeronave (vista de cima), com asa, EH e EV;
#   - vista lateral da asa comparando as cordas sem torcao e com a torcao
#     otimizada, que mostra o washout crescente da raiz para a ponta.
#
# Rodar da raiz do repo:  julia lab04_geometria.jl

using Plots

gr()

PAL = ["#2a78d6", "#eb6834", "#1baf7a"]
INK = "#0b0b0b"
INK2 = "#52514e"
GRIDC = "#e1e0d9"

estilo = (framestyle = :axes, gridcolor = GRIDC, gridalpha = 1.0,
          gridlinewidth = 0.7, foreground_color_axis = INK2,
          foreground_color_border = "#c3c2b7",
          foreground_color_text = INK2, tickfontsize = 10,
          guidefontsize = 11, legendfontsize = 9,
          background_color = :white)

# ---- parser das SECTIONs por superficie ----
function le_superficies(caminho)
    superficies = Dict{String, Vector{NTuple{5, Float64}}}()
    nome = ""
    esperando_nome = false
    lin = readlines(caminho)
    i = 1
    while i <= length(lin)
        t = strip(lin[i])
        if t == "SURFACE"
            esperando_nome = true
        elseif esperando_nome && !isempty(t) && !startswith(t, "#")
            nome = t
            superficies[nome] = NTuple{5, Float64}[]
            esperando_nome = false
        elseif t == "SECTION"
            j = i + 1
            while startswith(strip(lin[j]), "#"); j += 1; end
            vals = parse.(Float64, split(strip(lin[j])))
            push!(superficies[nome], (vals[1], vals[2], vals[3],
                                      vals[4], vals[5]))
            i = j
        end
        i += 1
    end
    return superficies
end

sup = le_superficies("avl/fwd.avl")

# contorno de uma superficie (LE ida + TE volta), espelhado em Y
function contorno(secs; espelha = true)
    xle = [s[1] for s in secs]; yle = [s[2] for s in secs]
    chord = [s[4] for s in secs]
    xs = vcat(xle, reverse(xle .+ chord))
    ys = vcat(yle, reverse(yle))
    if espelha
        xs = vcat(xs, reverse(xs)); ys = vcat(ys, -reverse(ys))
    end
    return xs, ys
end

# ---------------- planta (vista de cima) ----------------
p1 = plot(; xlabel = "x  [m]", ylabel = "y  [m]",
          title = "planta da aeronave", titlefontsize = 11,
          titlelocation = :left, legend = false, aspect_ratio = :equal,
          yflip = false, estilo...)
for (nome, cor) in (("Wing", PAL[1]), ("Horizontal tail", PAL[2]))
    xs, ys = contorno(sup[nome])
    plot!(p1, xs, ys; seriestype = :shape, fillcolor = cor,
          fillalpha = 0.25, linecolor = cor, linewidth = 1.8)
end
# EV: vista de cima e so a linha de centro (envergadura vertical)
vt = sup["Vertical tail"]
plot!(p1, [vt[1][1], vt[1][1] + vt[1][4], vt[end][1] + vt[end][4],
           vt[end][1]], [0, 0, 0, 0]; linecolor = PAL[3], linewidth = 2.5)

# ---------------- vista lateral: cordas com e sem torcao ----------------
asa = sup["Wing"]
p2 = plot(; xlabel = "x  [m]", ylabel = "z  [m]  (para cima)",
          title = "vista lateral da asa: cordas sem torção e com torção",
          titlefontsize = 11, titlelocation = :left, legend = :topright,
          aspect_ratio = :equal, estilo...)
for (k, s) in enumerate(asa)
    xle, _, zle, c, ainc = s
    # sem torcao: corda horizontal
    plot!(p2, [xle, xle + c], [zle, zle]; color = "#b9b8b3",
          linewidth = 2.0, label = (k == 1 ? "sem torção" : ""))
    # com torcao: corda girada de ainc em torno do bordo de ataque
    a = deg2rad(ainc)
    plot!(p2, [xle, xle + c*cos(a)], [zle, zle + c*sin(a)];
          color = PAL[1], linewidth = 2.2,
          label = (k == 1 ? "com torção otimizada" : ""))
end

fig = plot(p1, p2; layout = (2, 1), size = (1000, 900), dpi = 200,
           left_margin = 6Plots.mm, bottom_margin = 5Plots.mm)
mkpath("relatorio_lab04/images/01_geometria")
savefig(fig, "relatorio_lab04/images/01_geometria/geometria_avl.png")
println("figura: geometria_avl.png")
