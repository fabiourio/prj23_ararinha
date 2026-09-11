# Lab 03 — Projeto de perfil transônico | Equipe Ararinha

Campanha de otimização de perfil para a aeronave **B** (joelho da frente de
Pareto do Lab 02). O material do professor está intacto em `../eulerblock/`,
`../airfoil_mod/` e `../xfoil/`; tudo aqui é trabalho da equipe.

A definição completa do problema, com as justificativas de cada escolha, está
em `docs/superpowers/specs/2026-09-10-otimizacao-aerofolio-design.md`.

## Como rodar

```bash
pip install -r ../../requirements.txt
cd otimizacao_aerofolio/airfoil_mod/neysecco-airfoil_mod-*/ && pip install . && cd -

python test_xfoil_runner.py        # 9 testes de regressão  (~1 min)
python test_otimiza_secao.py       # 7 testes do otimizador (~10 s)

python verificacao_malha.py        # convergência de malha   (~30 min)
python verificacao_adjunto.py      # adjunto x dif. finitas  (~10 min)

python doe_clmax_xfoil.py          # DOE de cl_max           (~1 h)
python analise_doe_clmax.py        # escolhe o preditor
python calibra_limiar.py           # limiar por estação
python figuras_doe_clmax.py

python otimiza_secao.py meio       # ~20 min por estação
python otimiza_secao.py meio --sem-bluntez
python otimiza_secao.py raiz
python otimiza_secao.py ponta

python analisa_otimos.py           # comparação + XFoil + malha (~40 min)
python figuras_otimizacao.py       # itens 4, 5 e 6 do roteiro
python polares_finais.py           # itens 7 e 8
python compara_korn.py             # item 9 / realimentação
python tabelas_relatorio.py        # tabelas 1 e 2
```

As quatro otimizações são independentes e usam um núcleo cada — rodar em
paralelo custa o mesmo que rodar uma.

## Ponto de projeto

O perfil **não** vê M = 0,85. A asa tem 32,158° de enflechamento a meia corda,
e a mesma transformação que a equação de Korn do `designTool` já usa
internamente dá **M_n = M·cosΛ = 0,7196**. Peso médio de cruzeiro
(CL = 0,5053) e distribuição elíptica de sustentação — a elíptica é o alvo da
futura otimização de torção, então projetar para ela é projetar para a asa
que pretendemos ter.

| estação | η | (t/c)_n | cl_n | Re_n |
|---|---|---|---|---|
| raiz (junção asa-fuselagem) | 0,101 | 0,2179 | 0,600 | 4,34×10⁷ |
| meio (MAC — a seção do roteiro) | 0,398 | 0,1772 | 0,732 | 3,28×10⁷ |
| ponta (fim do tanque) | 0,900 | 0,1082 | 0,768 | 1,49×10⁷ |

## Resultados

| | raiz | meio | ponta |
|---|---|---|---|
| c_d partida | 0,0716 | 0,0560 | 0,0248 |
| c_d otimizado | 0,0199 | 0,0116 | 0,0089 |
| ganho | −72,3 % | −79,2 % | −64,3 % |
| ganho em malha fina | −71,0 % | −85,9 % | −81,7 % |
| cl_max no XFoil | 1,959 | 2,036 | 2,142 |
| cl_max real (viés −0,24) | 1,72 | **1,80** | 1,90 |
| restrição de cl_max | ativa | ativa | ativa |

Na MAC, a restrição de sustentação máxima **custa +9,2 % de arrasto e compra
0,99 de cl_max**: sem ela o perfil cai para cl_max = 1,042, e a aeronave
ficaria sem margem de empuxo na decolagem.

O `cl_max real` sai do viés medido em `afere_xfoil.py`: rodando perfis NACA
com valor experimental publicado pela nossa configuração, o XFoil lê **+0,24
acima** do experimental (mediana de cinco perfis). O α de estol confirma:
19-20,5° medidos contra 14-17° reais.

## Achados que valem o relatório

**O otimizador redescobriu o perfil supercrítico sozinho.** Partindo de um
NACA de 4 dígitos, moveu a espessura máxima de x/c = 0,30 para 0,359, deixou o
arqueamento dianteiro negativo e jogou a carga para trás. Na distribuição de
Mach, trocou um choque normal forte (M = 1,51 caindo a pique) por um platô
supersônico em M ≈ 1,15–1,20 com desaceleração gradual.

**O `k_korn = 0,95` do Lab 02 está essencialmente correto** — os perfis
otimizados dão k implícito de 0,894 a 0,940. O que falha é a **agregação**:
`(t/c)_m = 0,25·tcr + 0,75·tct = 0,109` é a espessura da estação η = 0,750,
onde a corda já caiu para 43 % da raiz. Integrando os c_d das seções ao longo
da envergadura dá CDwave = 0,00403 contra 0,00051 da fórmula — **7,9×**.

**A polar do perfil otimizado tem poço de arrasto estreito** (c_d sobe 57 %
abaixo do cl de projeto), a assinatura da otimização mono-ponto. Como o CL de
cruzeiro varia de 0,612 a 0,417 ao longo da missão, isso justifica uma rodada
multiponto.

**Os ótimos estavam no batente, não no interior** (`batentes_ativos.py`).
Com os batentes sugeridos pelo roteiro havia cinco batentes ativos, todos em
coeficientes do intradorso. Soltar o último — o que constrói o *cusp* côncavo
do bordo de fuga, que o RAE2822 tem em +0,052 e o roteiro proíbe — rende
**28 a 41 % de arrasto** por estação, e os ótimos passam a ser interiores.
Verificar se o ótimo é interior deveria vir **antes** de discutir o custo das
restrições: aqui a limitação dominante valia quatro vezes mais que a
restrição de cl_max, e estava invisível.

**Mas o cusp aproximadamente dobra o momento de arfagem.** Medido no Euler:
c_m de −0,051 para −0,150 na raiz, de −0,097 para −0,165 na MAC e de −0,109
para −0,169 na ponta.

O custo disso **não foi quantificado**, e por um motivo concreto: o
`designTool` não modela o `c_m` do perfil — o `balance.py` trata CG, ponto
neutro e margem estática, e o único momento presente é a inclinação dCm/dα da
fuselagem. Estimar arrasto de compensação exigiria uma análise de equilíbrio
que não fizemos (AVL ou balanço completo com CG reposicionado); o `c_m` 2D em
torno de c/4 não é o momento da aeronave em torno do CG.

Então o item 9 fica assim: o cusp entrega 29 a 54 % menos arrasto de seção e
dobra o momento — um **trade em aberto**. O que se pode afirmar é que a função
objetivo é 2D e não paga pelo momento, então o otimizador gasta essa moeda à
vontade. Sempre que uma grandeza relevante fica fora do objetivo, o ótimo se
desloca na direção de gastá-la.

## Três defeitos no material fornecido

1. **`export_airfoil` mata o XFoil em silêncio.** Com o padrão
   `close_te=True` ele repete o ponto do bordo de fuga, que o `cstfoil` já
   fechou — o terceiro ponto coincidente vira painel de comprimento zero e o
   XFoil encerra com código 0 e sem mensagem. Afeta o `single_run.py`.
2. **`cstfoil` devolve `x_max_thickness` errado** (0,117 em vez de 0,30 no
   NACA 1411): a linha 216 filtra o vetor de espessuras e a 219 indexa o de
   abscissas completo. Afeta a coluna `x_t/c,max` da Tab. 2.
3. **`cstfoil(Au, Al, …)` e `run_cst(Al, Au, …)` têm ordem invertida.**
   Trocar gera perfil de dentro para fora; o erro que chega ao Python é um
   `FileNotFoundError` que não diz nada.

## Duas armadilhas numéricas que custaram tempo

**O `c_d` do nível 1,0 de malha erra +70 %** (extrapolação de Richardson,
ordem observada 2,43). O bom é que o *ganho* da otimização **cresce** com o
refinamento (−78,9 % → −85,9 % na MAC), porque a dissipação numérica é um
piso quase absoluto que pesa proporcionalmente mais no perfil bom. Ou seja, o
nível 1,0 subestima o benefício. Mas os c_d **absolutos** que forem para o
`designTool` precisam vir da extrapolação.

**`ftol` tem que ficar acima do piso de ruído.** Com `ftol = 1e-6`, a raiz
convergiu na avaliação 16 e gastou mais 123 sem parar: o c_d oscila 1,5×10⁻⁵
depois de convergido, quinze vezes o `ftol`. Usamos 1e-5, o valor do próprio
professor.

**Restrição redundante quebra o subproblema do SLSQP.** Acrescentamos
`mint ≥ 0` (não cruzar superfícies) ao lado do `mint ≤ 0,01` do professor. As
duas linhas da jacobiana são `−dmint` e `+dmint`, exatamente antiparalelas,
deixando a matriz deficiente em posto **por construção** — o SLSQP morreu com
"Singular matrix E in LSQ subproblem" e entregou um c_d pior. E era
desnecessária: o `mint` converge prensado contra o teto de 0,01, ou seja o
otimizador empurra na direção oposta ao cruzamento. Restrições redundantes
não são inofensivas, mesmo inativas.

**Ótimos planos existem e o SLSQP não sai deles.** A ponta convergiu na
avaliação 21 e gastou mais 99: o t/c passeia de 0,1094 a 0,1111 com o c_d
constante na sétima casa. Não é `ftol` — são projetos distintos com o mesmo
arrasto.

**Não anote número à mão.** O primeiro balanço de trimagem usou `c_m` copiados
de rodadas parciais e concluiu que o cusp valia a pena (+13 % de saldo). Com
os valores lidos do Euler nos ótimos convergidos, o saldo é **−19 %**. O
atalho quase colocou a conclusão invertida no relatório.

## Arquivos

| arquivo | o quê |
|---|---|
| `xfoil_runner.py` | interface batch com o XFoil, com detecção de estol e repetição |
| `descritores.py` | descritores geométricos analíticos do perfil |
| `doe_clmax_xfoil.py` | DOE de cl_max (corte controlado, LHS, perfis grossos) |
| `analise_doe_clmax.py` | escolhe o preditor de cl_max |
| `calibra_limiar.py` | limiar de bluntez por estação |
| `otimiza_secao.py` | a otimização (SLSQP + adjunto) |
| `verificacao_malha.py` | convergência de malha |
| `verificacao_adjunto.py` | adjunto contra diferenças finitas |
| `testa_gradiente_clmax.py` | por que o XFoil não pode entrar no laço |
| `analisa_otimos.py` | comparação, verificação em XFoil e em malha fina |
| `compara_korn.py` | realimentação para o `designTool` |
| `figuras_*.py`, `polares_finais.py`, `tabelas_relatorio.py` | entregáveis |
| `estilo.py` | paleta (a mesma do Lab 02, verificada para daltonismo) |
| `test_*.py` | 16 testes de regressão |
