# O resultado da otimização faz sentido no contexto do projeto?
#
# Rodar de dentro da pasta avl/:   julia contexto_resultado.jl
#
# Três perguntas que não se respondem olhando só para a otimização:
#
#   1. O critério de estol está calibrado? O designTool estima CLmax da
#      aeronave limpa; se o nosso critério fizer a aeronave estolar num
#      CL bem diferente disso, ele está apertado ou frouxo demais, e a
#      restrição do Lab 04 herda esse erro.
#
#   2. O fator de Oswald que o AVL devolve é compatível com o de um
#      transporte real? Cuidado: existem dois "e" diferentes na
#      literatura e compará-los direto é erro de categoria.
#
#   3. O ganho de arrasto significa alguma coisa no total da aeronave?

using Printf
using JSON
using Statistics

const AVL   = "./avl.exe"
const BASE  = "aft.avl"
const SAIDA = "resultados"
const DADOS = JSON.parsefile("dados_designtool.json")
const PP    = DADOS["ponto_de_projeto"]
const POL   = DADOS["polar_designtool"]
const OT    = JSON.parsefile(joinpath(SAIDA, "torcao_otimizada.json"))

const SREF = PP["Sref"]
const BREF = DADOS["referencia"]["Bref"]
const SEMI = BREF/2
const AR   = BREF^2/SREF
const CL_PROJ = round(PP["CL"], digits = 4)
const MACH_BAIXO = 0.2
const ETA_AILERON = 0.56
const ETA_LIM   = [0.1011, 0.398, 0.90]
const CLMAX_LIM = [1.774, 1.7985, 1.7338]

function roda_avl(cmds)
    tmp = "_ct.txt"; write(tmp, cmds)
    try; return read(pipeline(`$AVL`, stdin = tmp), String)
    finally; rm(tmp, force = true); end
end
num(t, p) = (m = collect(eachmatch(p, t)); isempty(m) ? NaN :
             parse(Float64, m[end].captures[1]))

function escreve(destino, tw)
    L = readlines(BASE); n = length(L)
    ini = [i for i in 1:n if strip(L[i]) in ("SURFACE", "BODY")]
    i_asa = 0
    for i in ini
        strip(L[i]) == "SURFACE" || continue
        j = i + 1
        while j <= n && (isempty(strip(L[j])) || startswith(strip(L[j]), "#"))
            j += 1
        end
        j <= n && strip(L[j]) == "Wing" && (i_asa = i)
    end
    prox = findfirst(>(i_asa), ini)
    fim = prox === nothing ? n : ini[prox] - 1
    saida = String[]; i_sec, espera = 0, false
    for i in 1:n
        ln = L[i]; t = strip(ln); dentro = i_asa <= i <= fim
        if dentro && t == "SECTION"
            espera = true
        elseif dentro && espera && !isempty(t) && !startswith(t, "#")
            i_sec += 1; p = split(t); p[5] = @sprintf("%.4f", tw[i_sec])
            push!(saida, join(p, "  ")); espera = false; continue
        end
        push!(saida, ln)
    end
    write(destino, join(saida, "\n") * "\n")
end

function faixas_asa(s)
    tr = split(s, r"Surface # 1\s+Wing")[2]
    bl = String(split(tr, "Surface # 2")[1])
    L = [[parse(Float64, v) for v in split(strip(m.match))[2:end]]
         for m in eachmatch(r"^\s*\d+(?:\s+[-\dEe.+]+){12}\s*$"m, bl)]
    d = reduce(hcat, L)'
    return (eta = d[:, 1]./SEMI, cl_norm = d[:, 6])
end

function em_alfa(tw, alfa)
    escreve("_ct.avl", tw)
    s = roda_avl(join(["load _ct.avl", "oper", "m", "mn $MACH_BAIXO", "",
                       "d2 d2 0", "a a $alfa", "x", "fs", "", "",
                       "quit"], "\n") * "\n")
    rm("_ct.avl", force = true)
    return (CL = num(s, r"CLtot\s*=\s*([-\d.]+)"), fx = faixas_asa(s))
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

"Primeiro α em que alguma faixa alcança o clmax do seu perfil, e o CL ali."
function estol(tw; alfas = 6.0:0.5:26.0)
    dados = [(a, em_alfa(tw, a)) for a in alfas]
    eta = dados[1][2].fx.eta
    lim = clmax_local.(eta)
    melhor = (alfa = Inf, CL = NaN, eta = NaN)
    for j in eachindex(eta)
        cls = [d[2].fx.cl_norm[j] for d in dados]
        for k in 1:length(alfas)-1
            if cls[k] <= lim[j] <= cls[k+1]
                t = (lim[j] - cls[k])/(cls[k+1] - cls[k])
                a = alfas[k] + t*(alfas[k+1] - alfas[k])
                cl = dados[k][2].CL +
                     t*(dados[k+1][2].CL - dados[k][2].CL)
                a < melhor.alfa && (melhor = (alfa = a, CL = cl, eta = eta[j]))
                break
            end
        end
    end
    return melhor
end

tw0 = zeros(length(OT["etas"]))
tw2 = Float64.(OT["etapa_2_completa"]["torcoes"])
tw3 = Float64.(OT["etapa_3_projeto"]["torcoes"])

println("="^78)
println("1. O CRITÉRIO DE ESTOL ESTÁ CALIBRADO?")
println("="^78)
println("""
  O critério do Lab 04 diz que a asa estola quando a primeira faixa
  alcança o clmax do perfil daquela estação. Isso é uma aproximação
  clássica de projeto conceitual, mas ignora alívio tridimensional e
  trata o VLM invíscido como se resolvesse separação. O jeito de saber
  se ela está no lugar certo é comparar o CLmax que ela prevê com o que
  o designTool estima para a aeronave limpa.""")

@printf("\n  %-22s %9s %9s %9s\n", "caso", "α estol", "CL no estol", "estol η")
for (nome, tw) in (("sem torção", tw0), ("etapa 2", tw2), ("etapa 3", tw3))
    e = estol(tw)
    @printf("  %-22s %8.2f° %9.3f %9.3f\n", nome, e.alfa, e.CL, e.eta)
end
@printf("\n  designTool, CLmax da aeronave limpa: %.3f\n", POL["CLmax_limpa"])
println("""
  O designTool usa correlação de asa enflechada com os perfis; o nosso
  critério usa faixa a faixa. Se os dois caírem perto, a restrição do
  Lab 04 está pedindo o que a aeronave de fato entrega. Se o nosso CL de
  estol ficar bem ACIMA, o critério é frouxo; bem ABAIXO, é conservador
  e a restrição está cobrando caro por uma margem que não existe.""")

println("\n", "="^78)
println("2. O OSWALD É COMPATÍVEL COM UM TRANSPORTE REAL?")
println("="^78)
println("""
  Existem dois fatores de eficiência com o mesmo nome:

    e de ENVERGADURA (span efficiency), só do induzido invíscido. Para
    asa de transporte bem projetada fica entre 0,95 e 1,00, e pode
    passar de 1 com winglet ou diedro. É o que o AVL devolve.

    e de OSWALD da polar completa, que absorve também a parcela viscosa
    que cresce com a sustentação. Para jato de transporte fica tipico
    entre 0,75 e 0,85. É o que o designTool devolve.

  Comparar um com o outro direto é erro de categoria: são coisas
  diferentes e o segundo é sempre menor.""")

cdind_dt = POL["CDind"]
e_dt = POL["e"]
ar_eff = POL["AR_eff"]
e_dt_geom = CL_PROJ^2/(pi*AR*cdind_dt)
@printf("\n  designTool: e = %.4f referido a AR_eff = %.3f (com winglet)\n",
        e_dt, ar_eff)
@printf("              CDind = %.6f\n", cdind_dt)
@printf("              o mesmo CDind referido ao AR geométrico %.3f dá e = %.4f\n",
        AR, e_dt_geom)
@printf("\n  AVL, asa + fuselagem, sem as naceles: e = 0.9713  (ver carga_completa.jl)\n")
@printf("  AVL, aeronave completa compensada:    e = %.4f\n",
        OT["etapa_2_completa"]["e"])
println("""
  O modelo do AVL NÃO tem winglet, e o designTool conta com ele através
  do AR efetivo. Por isso o CDind do designTool é menor: a comparação
  honesta é contra a nossa configuração sem winglet.""")

println("\n", "="^78)
println("3. O GANHO SIGNIFICA ALGUMA COISA NO TOTAL?")
println("="^78)
cd0 = POL["CD0"]
cdwave = POL["CDwave"]
base_ff = 0.010982
for (nome, chave) in (("etapa 2, sem restrição de estol", "etapa_2_completa"),
                      ("etapa 3, projeto", "etapa_3_projeto"))
    cdff = OT[chave]["CDff"]
    d = 1e4*(cdff - base_ff)
    total = cd0 + cdwave + base_ff
    @printf("  %-34s %+7.2f counts   %+.2f%% do arrasto total\n", nome, d,
            100*(cdff - base_ff)/total)
end
@printf("\n  arrasto total da aeronave no ponto de projeto: %.5f", cd0 + cdwave + base_ff)
@printf("  (%.0f counts)\n", 1e4*(cd0 + cdwave + base_ff))
@printf("  parasita %.0f, onda %.0f, induzido sem torção %.0f counts\n",
        1e4*cd0, 1e4*cdwave, 1e4*base_ff)
println("""
  Regra de bolso de transporte: 1% de arrasto de cruzeiro vale perto de
  1% de combustível na etapa. A leitura honesta é essa, e não o número
  em counts isolado.""")

println("\n", "="^78)
println("4. A RESTRIÇÃO DE ESTOL ESTÁ PESADA?")
println("="^78)
@printf("  sem torção, não conforme:        CDff %.6f  (estol em η 0,842)\n",
        base_ff)
@printf("  etapa 2, ignora o estol:         CDff %.6f  (%+.2f counts)\n",
        OT["etapa_2_completa"]["CDff"],
        1e4*(OT["etapa_2_completa"]["CDff"] - base_ff))
if haskey(OT, "custo_da_margem")
    cm = OT["custo_da_margem"]
    for k in sort([parse(Float64, x) for x in keys(cm)])
        @printf("  conforme com margem de %.1f°:      CDff %.6f  (%+.2f counts)\n",
                k, cm[string(k)], 1e4*(cm[string(k)] - base_ff))
    end
end
println("""
  A leitura: a asa sem torção JÁ É não conforme, então a comparação
  relevante não é contra o ótimo sem restrição e sim contra a linha de
  base. Com margem nula a torção entrega conformidade E arrasto menor;
  é só a partir de meio grau que ela começa a cobrar.""")
