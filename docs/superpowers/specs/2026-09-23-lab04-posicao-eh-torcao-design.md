# Lab 04 — posição longitudinal da EH e torção da asa

Data: 23/09/2026 · Entrega do lab: 25/09/2026

## Contexto

O professor, olhando o AVL da equipe, apontou dois problemas: a asa sombreia
a EH em alto ângulo de ataque, e a EH pode ir mais para trás para aproveitar
o braço e reduzir o MTOW. A orientação foi recuar a EH manualmente até o limite
da fuselagem, já que a otimização do Lab 02 já rodou. O Lab 04 recomeça do
zero: só os `.avl`, os perfis e o `avl.exe` continuam na pasta.

Hoje (`Lc_h` = 4,6) a esteira da raiz da asa cruza a altura da EH por volta de
α ≈ 11°, praticamente no estol.

## Decisões já tomadas

- Scripts em **Python**, que chamam o designTool diretamente. É uma linguagem
  só dentro do zip.
- **O deslocamento em Z do AVL fica como está** (asa −1,20 m, EH +1,85 m,
  EV +0,85 m, nacele acompanhando a asa). É exigência do professor, para
  manter distância entre os painéis. Por isso o sombreamento é medido com as
  alturas **reais** do designTool, e nunca com as do AVL.
- A posição da EH é escolhida por **varredura** de `Lc_h`, não por otimização
  conjunta. A torção é otimizada depois, com a EH já fixada.
- Objetivo da torção: **mínimo CD compensado (Cm = 0) no ponto de projeto**.

## Etapa 0 — verificações antes de usar os resultados

### 0a. Convergência de malha

Roda sobre a geometria de base (`Lc_h` = 4,6, sem torção, `aft.avl`), no ponto
de projeto compensado (M = 0,85, CL de projeto, Cm = 0 pelo `it`).

- Refina uma direção por vez, mantendo as outras na malha atual: painéis ao
  longo da envergadura e da corda da asa, do winglet e da EH. O corpo
  (`Nbody`) entra só se as superfícies já tiverem convergido.
- Monitora CDff (em counts), CL, Cm, o `it` de compensação e o cl da faixa
  mais crítica no estol (M = 0,2).
- Critério de platô: a malha adotada é a menor em que dobrar o número de
  painéis muda o CDff menos de 0,1 count e o `it` menos de 0,02°.
- Saída: gráficos de convergência e a tabela da malha adotada. O
  `gera_avl.py` passa a usar essa malha.
- A malha é conferida de novo no `Lc_h` escolhido, com um único refinamento.

### 0b. Verificação do otimizador: asa limpa → distribuição elíptica

O otimizador da etapa 2 é aplicado a um caso com resposta conhecida.

- Modelo: só a asa, sem fuselagem, naceles, empenagens ou winglet, com
  diedro zero (asa plana, para valer a teoria de Prandtl/Munk). A planta, o
  enflechamento e os perfis são os do projeto.
- O mesmo código da etapa 2 (mesmas variáveis de torção, limites e SLSQP)
  minimiza o CDi no CL de projeto, sem restrição de estol e sem compensação.
- Esperado: carregamento c·cl/c_ref × y elíptico, e fator de Oswald do plano
  de Trefftz próximo de 1.
- Critérios de aprovação:
  - e ≥ 0,98;
  - desvio RMS do carregamento em relação à elipse de mesma sustentação
    ≤ 2% do carregamento na raiz;
  - duas partidas diferentes (sem torção e com washout linear de −4°) chegam
    ao mesmo CDi, com diferença ≤ 0,1 count.
- Saída: gráfico do carregamento otimizado contra a elipse, e um e antes e
  depois da otimização.
- Se o e ficar abaixo de 0,98, verifica-se primeiro se a culpa é da
  parametrização (poucas seções de torção). Para isso a mesma otimização roda
  com o dobro de seções, antes de se concluir que o otimizador falhou.

## Etapa 1 — varredura da posição da EH

Variável: `Lc_h` de 4,0 até o limite geométrico, com passo de 0,1, mais o
próprio limite. O `Cht` fica fixo, então `S_h` cai com o braço, como no
designTool. A asa fica sem torção nesta etapa.

Para cada `Lc_h`:

1. designTool: W0, S_h e a geometria da EH, CD0, CG dianteiro e traseiro, e
   dε/dα (Roskam, a mesma fórmula do `balance.py`).
2. Peso no ponto de projeto, com a mesma definição do Lab 03 (50% de
   combustível e 100% de carga paga), recalculado → CL de projeto.
   M = 0,85, h = 10.668 m.
3. `.avl` regerado (fwd e aft) com a EH nova e o mesmo deslocamento em Z.
4. AVL no cruzeiro: `it` que zera o δe (item 3 do roteiro), CD compensado,
   ponto neutro → margem estática nos dois CGs.
5. AVL em baixa velocidade (M = 0,2): α_estol e CLmax pelo método da seção
   crítica, e δe de compensação no CLmax com o CG dianteiro.
6. Sombreamento (ver abaixo), com α_estol.

### Critério de sombreamento

A "sombra" da asa é a faixa entre as retas que partem do bordo de ataque e do
bordo de fuga de uma seção da asa na direção do escoamento local. Essa direção
fica inclinada de γ = α − ε em relação ao eixo do corpo, com ε = (dε/dα)·α.
Para cada estação y da semienvergadura da EH, usa-se a seção da asa na mesma
y, com x e z reais. A EH está **fora da sombra** quando o bordo de ataque e o
extradorso dela ficam abaixo da reta do bordo de fuga. Como γ cresce com α, a
reta só sobe, então depois que a EH sai da sombra ela não volta.

Saídas: α_saída, o α a partir do qual toda a EH fica fora da sombra, e a folga
α_estol − α_saída.

### Restrições para escolher `Lc_h`

| Restrição | Valor adotado |
|---|---|
| Fim da fuselagem | bordo de fuga da raiz da EH ≤ L_f − 0,5 m |
| Sombreamento | α_saída ≤ α_estol − 2° |
| Profundor | \|δe\| ≤ 20° na compensação no CLmax com o CG dianteiro |
| Margem estática | MS ≥ 5% com o CG traseiro (ponto neutro do AVL) |
| Leme | fração do leme encoberta pelas linhas de 60°/30° do `plots.py`, só informada |

Escolha: o **maior `Lc_h` viável**. W0 e CD devem cair de forma monótona com o
braço; se não caírem, a escolha passa a ser o mínimo de CD entre os pontos
viáveis. O resultado é um gráfico de W0, CD compensado, δe, MS e folga de
sombreamento × `Lc_h`, com a região viável marcada.

## Etapa 2 — torção na posição escolhida

- Variáveis: `Ainc` das seções da asa (a raiz fica em 0°), limites de −8° a
  +3°. O winglet continua como superfície própria, no COMPONENT 1, com toe
  zero.
- Objetivo: CD compensado no CL de projeto, M = 0,85, `aft.avl`.
- Restrição de estol: no método da seção crítica (M = 0,2), a primeira faixa
  a estolar fica para dentro de η = 0,56 (raiz do aileron), com desempate de
  0,2° medido direto no AVL.
- Otimizador: SLSQP (scipy) com diferenças finitas. Cada avaliação roda o AVL.
- Conferência de independência: a otimização roda também em 2 valores de
  `Lc_h` vizinhos. Se a torção ótima mudar pouco, a sequência varredura →
  torção fica justificada. Se mudar muito, a varredura é refeita com a torção
  ótima.

## Etapa 3 — roteiro do Lab 04

Itens 1–8 e as derivadas de estabilidade, com a geometria final. Reaproveita os
módulos das etapas 1 e 2. Tem plano próprio.

## Estrutura (pasta `avl/`)

| arquivo | responsabilidade |
|---|---|
| `aeronave.py` | designTool → dicionário com geometria, pesos, CGs, ponto de projeto e dε/dα para um `Lc_h` |
| `gera_avl.py` | escreve `fwd.avl`/`aft.avl` a partir desse dicionário e de um vetor de torção, com o deslocamento em Z |
| `avl_run.py` | roda o `avl.exe` por stdin e lê `ft`, `fs` e `st` |
| `sombra.py` | critério de sombreamento (geometria pura, testável sem AVL) |
| `convergencia_malha.py` | etapa 0a |
| `verifica_eliptica.py` | etapa 0b |
| `varredura_eh.py` | etapa 1: tabela, gráfico e `Lc_h` escolhido |
| `otimizacao_torcao.py` | etapa 2 |

Testes: `sombra.py` com casos geométricos construídos à mão. `gera_avl.py`
regenera o `fwd.avl` atual (com `Lc_h` = 4,6 e sem torção) e o resultado tem de
ser igual ao do commit `11c6802`. `avl_run.py` reproduz um `ft` conhecido.

## Riscos

- dε/dα pelo Roskam é uma estimativa. O critério usa uma esteira sem
  espessura própria, então os 2° de margem cobrem esse erro. Isso deve ser
  declarado no relatório.
- O `standard_airplane.py` muda com o `Lc_h` escolhido. É uma mudança de
  projeto e deve ser registrada no commit.
