# Distribuição de sustentação da aeronave INTEIRA e de onde vem o fator
# de Oswald abaixo de 1 na etapa 2.
#
# Rodar de dentro da pasta avl/:   julia carga_completa.jl
#
# A etapa 2 minimiza o arrasto induzido da aeronave completa e compensada e
# ainda assim devolve Oswald 0,86, bem abaixo do 1,00 da asa isolada. A
# pergunta é se isso é esperado.
#
# A explicação de manual seria arrasto de compensação: a empenagem carrega
# para baixo, a asa precisa carregar mais que o peso e o induzido vai com
# o quadrado, de modo que o e cai mesmo com carga elíptica. ESSA
# EXPLICAÇÃO NÃO SE APLICA AQUI, e os números deste script mostram por
# quê: no ótimo da etapa 2 a empenagem está praticamente descarregada e a
# asa carrega MENOS que o CLff, de modo que o fator de compensação
# empurraria o e para CIMA de 1.
#
# A causa real é outra, e aparece na decomposição por peça: as naceles.
# Elas estão modeladas como SUPERFÍCIE sustentadora em forma de anel, com
# cerca de 91 m² de malha cada, quase metade da área da asa. No VLM um
# anel em ângulo de ataque gera circulação e esteira, e isso vira arrasto
# induzido: sozinhas elas respondem por 9,6 dos 11,6 counts que separam a
# asa isolada da aeronave completa, e ainda tiram 11% da carga local da
# asa na estação onde ficam.
#
# O script mede as duas coisas: a carga por superfície e a decomposição
# do arrasto montando a aeronave peça por peça.

using Printf
using JSON
using Plots
using Statistics

gr()

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

const PAL   = ["#2a78d6", "#eb6834", "#1baf7a"]
const INK, INK2 = "#0b0b0b", "#52514e"
const CINZA = "#b9b8b3"
const ESTILO = (framestyle = :axes, gridcolor = "#e1e0d9", gridalpha = 1.0,
                gridlinewidth = 0.7, foreground_color_axis = INK2,
                foreground_color_border = "#c3c2b7",
                foreground_color_text = INK2, tickfontsize = 9,
                guidefontsize = 10, legendfontsize = 8,
                background_color = :white)

function roda_avl(cmds)
    tmp = "_cc_cmds.txt"
    write(tmp, cmds)
    try
        return read(pipeline(`$AVL`, stdin = tmp), String)
    finally
        rm(tmp, force = true)
    end
end

num(t, p) = (m = collect(eachmatch(p, t)); isempty(m) ? NaN :
             parse(Float64, m[end].captures[1]))

function escreve(destino, tw)
    linhas = readlines(BASE); n = length(linhas)
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
    saida = String[]; i_sec, espera = 0, false
    for i in 1:n
        ln = linhas[i]; t = strip(ln)
        dentro = i_asa <= i <= fim
        if dentro && t == "SECTION"
            espera = true
        elseif dentro && espera && !isempty(t) && !startswith(t, "#")
            i_sec += 1; p = split(t); p[5] = @sprintf("%.4f", tw[min(i_sec, length(tw))])
            push!(saida, join(p, "  ")); espera = false; continue
        end
        push!(saida, ln)
    end
    write(destino, join(saida, "\n") * "\n")
end

"""
Faixas de TODAS as superfícies, e não só da asa.

Cada bloco de `Surface #` traz a tabela de faixas daquela superfície. Com
YDUPLICATE o AVL repete a superfície espelhada, então basta o primeiro
bloco de cada nome: a carga do outro lado é a mesma.
"""
function faixas_todas(saida)
    blocos = split(saida, "Surface #")
    out = Dict{String,NamedTuple}()
    for b in blocos[2:end]
        linhas = split(b, "\n")
        cab = strip(linhas[1])
        nome = strip(join(split(cab)[2:end], " "))
        occursin("YDUP", uppercase(nome)) && continue
        haskey(out, nome) && continue
        L = [[parse(Float64, v) for v in split(strip(m.match))[2:end]]
             for m in eachmatch(r"^\s*\d+(?:\s+[-\dEe.+]+){12}\s*$"m, String(b))]
        isempty(L) && continue
        d = reduce(hcat, L)'
        d = d[d[:, 1] .< SEMI - 1e-3, :]   # sem as faixas do winglet, todas em y = b/2
        out[nome] = (y = d[:, 1], corda = d[:, 2], area = d[:, 3],
                     ccl = d[:, 4], cl_norm = d[:, 6], cl = d[:, 7])
    end
    return out
end

function roda(tw; trim = true)
    escreve("_cc.avl", tw)
    c = ["load _cc.avl", "oper", "m", "mn $MACH", ""]
    trim && push!(c, "d2 pm 0")
    append!(c, ["a c $CL_PROJ", "x", "fs", "", "fn", "", "", "quit"])
    s = roda_avl(join(c, "\n") * "\n")
    rm("_cc.avl", force = true)
    return (CDff = num(s, r"CDff\s*=\s*([-\d.Ee+]+)"),
            CDind = num(s, r"CDind\s*=\s*([-\d.Ee+]+)"),
            CL = num(s, r"CLtot\s*=\s*([-\d.]+)"),
            CLff = num(s, r"CLff\s*=\s*([-\d.Ee+]+)"),
            e = num(s, r"\se =\s+([-\d.]+)"),
            de = num(s, r"elevator\s*=\s*([-\d.]+)"),
            alpha = num(s, r"Alpha\s*=\s*([-\d.]+)"),
            sup = forcas_superficie(s), fx = faixas_todas(s))
end

"""
Forças por superfície, lidas da tabela que o próprio AVL monta.

Integrar c cl dy à mão só valeria para uma superfície plana e horizontal.
A nacele deste modelo é um anel, com as seções dispostas num círculo, de
modo que as faixas dela nem estão num plano: somar c cl dy ali dá um
número sem significado, e grande. A tabela do AVL resolve as direções
corretamente, então é ela que vale.
"""
function forcas_superficie(saida)
    bloco = split(saida, "Surface Forces (referred to Sref")
    length(bloco) < 2 && return Dict{String,NamedTuple}()
    out = Dict{String,NamedTuple}()
    for ln in split(bloco[2], "\n")
        m = match(r"^\s*\d+\s+((?:[-\d.Ee+]+\s+){9})(\S.*?)\s*$", ln)
        m === nothing && continue
        v = [parse(Float64, x) for x in split(m.captures[1])]
        nome = replace(m.captures[2], "(YDUP)" => "") |> strip |> String
        a = get(out, nome, (area = 0.0, CL = 0.0, CDi = 0.0))
        out[nome] = (area = a.area + v[1], CL = a.CL + v[2],
                     CDi = a.CDi + v[8])
    end
    return out
end

eliptica(η) = sqrt.(max.(0.0, 1 .- η.^2))

tw0 = zeros(length(OT["etas"]))
tw2 = Float64.(OT["etapa_2_completa"]["torcoes"])
tw3 = Float64.(OT["etapa_3_projeto"]["torcoes"])

println("="^78)
println("DISTRIBUIÇÃO DE SUSTENTAÇÃO DA AERONAVE COMPLETA")
println("="^78)

casos = [("sem torção", tw0), ("etapa 2, irrestrito", tw2),
         ("etapa 3, projeto", tw3)]
res = Dict{String,Any}()
for (nome, tw) in casos
    r = roda(tw)
    res[nome] = r
    println("\n", nome)
    @printf("  CLtot %.4f   CLff %.4f   CDff %.6f   Oswald %.4f   α %.2f°   profundor %.2f°\n",
            r.CL, r.CLff, r.CDff, r.e, r.alpha, r.de)
    soma = 0.0
    for sup in sort(collect(keys(r.sup)))
        f = r.sup[sup]
        soma += f.CL
        @printf("    %-20s CL %+8.4f   CDi %+9.6f\n", sup, f.CL, f.CDi)
    end
    @printf("    %-20s CL %+8.4f   (o que falta para CLtot %+.4f é a fuselagem)\n",
            "soma das superfícies", soma, r.CL)
end

# ====================================================================
println("\n", "="^78)
println("POR QUE O OSWALD FICA ABAIXO DE 1")
println("="^78)
println("""
  O e que o AVL imprime é e = CLff²/(π AR_asa CDff): usa o CL do plano de
  Trefftz no numerador e o alongamento da ASA no denominador, mas o CDff
  é o da AERONAVE INTEIRA. Ou seja, ele cobra da asa o arrasto induzido
  que o conjunto faz. Três coisas empurram esse número para baixo sem que
  a asa esteja mal otimizada:

    1. a empenagem carrega para fechar a arfagem, de modo que a asa tem
       de carregar mais do que o peso, e o induzido vai com o quadrado;
    2. a empenagem e a nacele fazem induzido próprio, que entra no CDff;
    3. asa e empenagem estão em alturas diferentes, o sistema é NÃO
       PLANAR no plano de Trefftz e as esteiras interferem.

  O primeiro efeito é quantificável: se a asa carrega k vezes o CLff, o
  induzido dela vai a k² e o e cai para 1/k² mesmo com carga elíptica.""")

@printf("\n  %-22s %8s %8s %8s %8s %9s %9s\n", "caso", "CL asa", "CL EH",
        "CL nac", "k asa", "1/k²", "e do AVL")
for (nome, _) in casos
    r = res[nome]
    pega(chave) = (k = findfirst(n -> occursin(chave, n), sort(collect(keys(r.sup))));
                   k === nothing ? 0.0 : r.sup[sort(collect(keys(r.sup)))[k]].CL)
    cl_asa, cl_eh, cl_nac = pega("Wing"), pega("Horizontal"), pega("Nacelle")
    k = cl_asa/r.CLff
    @printf("  %-22s %+8.4f %+8.4f %+8.4f %8.4f %9.4f %9.4f\n", nome,
            cl_asa, cl_eh, cl_nac, k, 1/k^2, r.e)
end

let r = res["etapa 2, irrestrito"]
    cdi_nao_asa = sum(v.CDi for (k, v) in r.sup if !occursin("Wing", k))
    @printf("  induzido de campo próximo que NÃO é da asa: %.6f (%.1f counts)\n",
            cdi_nao_asa, 1e4*cdi_nao_asa)
end

# ====================================================================
println("\n", "="^78)
println("DECOMPOSIÇÃO: O QUE CADA PEÇA ACRESCENTA AO INDUZIDO")
println("="^78)
# A conta acima mostra que NÃO é arrasto de compensação: na etapa 2 a
# empenagem está praticamente descarregada e a asa carrega menos que o
# CLff, de modo que o fator 1/k² empurraria o e para cima de 1. Para
# achar a causa, o jeito direto é montar a aeronave por partes e medir o
# e a cada peça acrescentada, sempre com a MESMA torção e o mesmo CL.

"Escreve o modelo mantendo apenas as superfícies e corpos pedidos."
function escreve_parcial(destino, tw, manter)
    linhas = readlines(BASE); n = length(linhas)
    ini = [i for i in 1:n if strip(linhas[i]) in ("SURFACE", "BODY")]
    cabecalho = 1:(isempty(ini) ? n : ini[1] - 1)
    saida = String[linhas[i] for i in cabecalho]
    i_sec = 0
    for (b, i0) in enumerate(ini)
        i1 = b < length(ini) ? ini[b+1] - 1 : n
        j = i0 + 1
        while j <= i1 && (isempty(strip(linhas[j])) ||
                          startswith(strip(linhas[j]), "#")); j += 1; end
        nome = j <= i1 ? strip(linhas[j]) : ""
        e_asa = nome == "Wing"
        (nome in manter) || continue
        espera = false
        for i in i0:i1
            ln = linhas[i]; t = strip(ln)
            if e_asa && t == "SECTION"
                espera = true
            elseif e_asa && espera && !isempty(t) && !startswith(t, "#")
                i_sec += 1; p = split(t); p[5] = @sprintf("%.4f", tw[min(i_sec, length(tw))])
                push!(saida, join(p, "  ")); espera = false; continue
            end
            push!(saida, ln)
        end
    end
    write(destino, join(saida, "\n") * "\n")
end

function roda_parcial(tw, manter; trim = false)
    escreve_parcial("_cp.avl", tw, manter)
    c = ["load _cp.avl", "oper", "m", "mn $MACH", ""]
    trim && push!(c, "d2 pm 0")
    append!(c, ["a c $CL_PROJ", "x", "", "quit"])
    s = roda_avl(join(c, "\n") * "\n")
    rm("_cp.avl", force = true)
    return (CDff = num(s, r"CDff\s*=\s*([-\d.Ee+]+)"),
            CLff = num(s, r"CLff\s*=\s*([-\d.Ee+]+)"),
            e = num(s, r"\se =\s+([-\d.]+)"))
end

# A lista de peças é conferida contra as superfícies do próprio modelo, de
# modo que uma superfície nova no aft.avl não fica de fora sem aviso.
const ASA = ["Wing", "Winglet"]
const VARIANTES = [
    ("só a asa, sem winglet",          ["Wing"], false),
    ("asa com winglet",                ASA, false),
    ("asa + fuselagem",                [ASA; "Fuselage"], false),
    ("asa + fuselagem + naceles",      [ASA; "Fuselage"; "Nacelle"], false),
    ("tudo, profundor em zero",        [ASA; "Fuselage"; "Nacelle";
                                        "Horizontal tail"; "Vertical tail"], false),
    ("tudo, compensado (caso de projeto)", [ASA; "Fuselage"; "Nacelle";
                                         "Horizontal tail"; "Vertical tail"], true)]
let nomes = String[]
    L = readlines(BASE)
    for (i, ln) in enumerate(L)
        strip(ln) in ("SURFACE", "BODY") || continue
        j = i + 1
        while isempty(strip(L[j])) || startswith(strip(L[j]), "#"); j += 1; end
        push!(nomes, strip(L[j]))
    end
    faltam = setdiff(nomes, VARIANTES[end][2])
    isempty(faltam) || error("a decomposição não inclui: $(join(faltam, ", "))")
end

@printf("\n  torção da etapa 2, CL fixo em %.4f\n\n", CL_PROJ)
@printf("  %-34s %10s %9s %9s\n", "configuração", "CDff", "e", "Δ counts")
let anterior = nothing
    for (nome, manter, trim) in VARIANTES
        r = roda_parcial(tw2, manter; trim = trim)
        d = anterior === nothing ? 0.0 : 1e4*(r.CDff - anterior)
        @printf("  %-34s %10.6f %9.4f %+9.2f\n", nome, r.CDff, r.e, d)
        anterior = r.CDff
    end
end
println("""
  A coluna Δ é o que cada peça acrescenta ao arrasto induzido da
  aeronave inteira, com a torção congelada na da etapa 2. É a resposta
  direta para de onde vem o Oswald abaixo de 1.

  A nacele deste modelo é uma SUPERFÍCIE sustentadora em forma de anel,
  com cerca de 91 m² de malha cada, quase metade da área da asa. No VLM
  um anel em ângulo de ataque gera circulação, esteira e portanto
  arrasto induzido. Uma nacele real tem escoamento passante e não se
  comporta assim. A modelagem foi mantida porque é a do 737.avl da
  disciplina, inclusive o COMPONENT 1 junto com a asa, e com ela o 737
  de referência também fica com e perto de 0,75.""")

# A pergunta que importa para o projeto: a nacele distorce a CARGA da
# asa? Se distorcer, a torção ótima foi ajustada em cima de um artefato.
println("\n  a nacele distorce a carga da asa? (mesma torção, mesmo CL)")
function carga_asa(manter)
    escreve_parcial("_cp.avl", tw2, manter)
    c = ["load _cp.avl", "oper", "m", "mn $MACH", "", "a c $CL_PROJ", "x",
         "fs", "", "", "quit"]
    s = roda_avl(join(c, "\n") * "\n")
    rm("_cp.avl", force = true)
    f = faixas_todas(s)
    k = findfirst(n -> occursin("Wing", n), sort(collect(keys(f))))
    return f[sort(collect(keys(f)))[k]]
end
com = carga_asa([ASA; "Fuselage"; "Nacelle"])
sem = carga_asa([ASA; "Fuselage"])
if length(com.y) == length(sem.y)
    rel = (com.ccl .- sem.ccl) ./ maximum(sem.ccl)
    j = argmax(abs.(rel))
    @printf("    maior diferença de carga local: %+.1f%% em η = %.3f\n",
            100*rel[j], com.y[j]/SEMI)
    @printf("    diferença média em módulo: %.1f%%\n", 100*mean(abs.(rel)))
    println(maximum(abs.(rel)) > 0.05 ?
        "    >> a nacele MEXE na carga da asa, então ela influencia a torção ótima" :
        "    >> a nacele quase não mexe na carga da asa: a torção ótima fica de pé")
end

# ====================================================================
# FIGURA
# ====================================================================
# As duas superfícies num mesmo eixo de envergadura absoluta, que é o que
# mostra a empenagem ocupando a mesma faixa de y da parte interna da asa.

p1 = plot(; xlabel = "y [m]", ylabel = "c·cl local [m]",
          title = "(a) carga por envergadura, asa e empenagem",
          titlefontsize = 10, titlelocation = :left, legend = :topright,
          ESTILO...)
hline!(p1, [0.0]; color = "#c3c2b7", linewidth = 0.8, label = "")
for (i, (nome, _)) in enumerate(casos)
    r = res[nome]
    for (sup, f) in r.fx
        occursin("Wing", sup) || occursin("Horizontal", sup) || continue
        plot!(p1, f.y, f.ccl; color = PAL[i], linewidth = 2.2,
              linestyle = occursin("Wing", sup) ? :solid : :dash,
              label = occursin("Wing", sup) ? nome : "")
    end
end
annotate!(p1, 6.0, -0.35,
          text("tracejado: empenagem horizontal", 8, INK2, :left))

# carga da asa normalizada, contra a elíptica
p2 = plot(; xlabel = "η = 2y/b", ylabel = "carga normalizada",
          title = "(b) só a asa, contra a elíptica", titlefontsize = 10,
          titlelocation = :left, legend = :bottomleft, xlims = (0, 1),
          ESTILO...)
for (i, (nome, _)) in enumerate(casos)
    r = res[nome]
    for (sup, f) in r.fx
        occursin("Wing", sup) || continue
        η = f.y ./ SEMI
        plot!(p2, η, f.ccl ./ (4*SREF*r.CL/(pi*BREF)); color = PAL[i],
              linewidth = 2.2, label = nome)
        i == 1 && plot!(p2, η, eliptica(η); color = INK, linewidth = 2,
                        linestyle = :dash, label = "elíptica")
    end
end

savefig(plot(p1, p2; layout = (1, 2), size = (1250, 480), dpi = 200,
             left_margin = 6Plots.mm, bottom_margin = 6Plots.mm),
        joinpath(SAIDA, "carga_completa.png"))
println("\nfigura: $SAIDA/carga_completa.png")
