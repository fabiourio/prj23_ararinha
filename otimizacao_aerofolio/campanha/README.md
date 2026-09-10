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
| c_d otimizado | 0,0199 | 0,0118 | 0,0089 |
| ganho (malha fina) | −71,0 % | −85,9 % | −81,7 % |
| cl_max verificado no XFoil | — | 2,058 | 2,128 |
| restrição de cl_max | ativa | ativa | ativa |

Na MAC, a restrição de sustentação máxima **custa +10,9 % de arrasto e compra
0,79 de cl_max**: sem ela o perfil cai para cl_max = 1,264, e a aeronave
ficaria sem margem de empuxo na decolagem.

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
