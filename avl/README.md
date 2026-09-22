# Lab 04 -- análise aerodinâmica no AVL

## Como rodar

Tudo roda de dentro desta pasta, nesta ordem:

```
python dados_designtool.py      # ponte com o designTool + verificação da geometria
julia  convergencia_malha.jl    # define e justifica a malha
julia  lab04.jl                 # análises do roteiro
julia  refina_asa.py            # refina a asa para 11 seções (só se mudar a geometria)
julia  otimizacao_torcao.jl     # otimização de torção (estudo próprio)
julia  verifica_torcao.jl       # verificação independente do resultado acima
julia  carga_completa.jl       # carga por superfície e de onde vem o Oswald
```

O `torcao_monotonica.jl` guarda o modelo em `resultados/modelo_torcao.json` e
só o reconstrói se o `aft.avl` for mais novo que o cache. Construir o modelo
custa cerca de 80 rodadas do AVL, rodar os casos em cima dele custa segundos.

O `dados_designtool.py` precisa ser rodado primeiro: ele gera o
`dados_designtool.json`, de onde os scripts Julia leem o ponto de projeto e
os dados de referência. É o único lugar que chama o designTool.

## Modelo

| arquivo | o que é |
|---|---|
| `aft.avl` / `fwd.avl` | aeronave com CG traseiro e dianteiro. Diferem só na linha `Xref` |
| `airfoils/` | perfis otimizados no Lab 03 (família roteiro): raiz, centro, ponta |
| `fuselage_nondim.dat` | contorno da fuselagem (`BODY` do AVL) |
| `737.avl` | modelo de referência da disciplina, para consulta |
| `avl.exe` | executável da disciplina (versão 3.37) |

**Malha adotada:** asa 8 × 20, empenagem horizontal 16 × 10, empenagem
vertical 8 × 6, corpo com 20 painéis. Total de 832 vórtices, cerca de
0,8 s por rodada de cinco pontos. A envergadura foi definida pelo arrasto
e a corda pelos momentos e pela posição efetiva da charneira do profundor
(ver `convergencia_malha.jl`).

**Alturas em Z:** as superfícies foram deslocadas para não cruzar a
fuselagem, que ocupa de z = −2,348 a +3,732. Asa −1,20 m, empenagem
horizontal +1,85 m, empenagem vertical +0,85 m, naceles acompanhando a asa.
Os deslocamentos são verificados pelo `dados_designtool.py`.

## Numeração dos controles

A numeração do AVL segue a ordem de aparição no arquivo e **difere** da do
roteiro do professor (que usa o `b737mod.avl`):

| roteiro | aqui | o que é |
|---|---|---|
| `d4 pm 0.0` | `d2 pm 0.0` | trimagem pelo profundor |
| menu `de`: `2 <valor>` | menu `de`: `1 <valor>` | incidência da empenagem |

Controles: `d1` aileron, `d2` profundor, `d3` leme. Variável de projeto:
`1` = incidência da empenagem (`it`).

## Ponto de projeto

O mesmo do Lab 03, o peso médio de cruzeiro, para que o CL da aeronave
coincida com o CL para o qual os perfis foram otimizados:
M = 0,85, h = 10.668 m, W = 229.669 kgf, **CL = 0,5053**.

Nos comandos do AVL: `m` → `mn 0.85` (ou `mn 0.2` para baixa velocidade) e
`a c 0.5053`.

## Resultados

Tudo que os scripts geram vai para `resultados/`:

| arquivo | de onde vem |
|---|---|
| `conv_corda.png`, `conv_envergadura.png`, `malha_adotada.txt` | convergência de malha |
| `ponto_projeto_carga.png`, `ponto_de_projeto.txt` | ponto de projeto |
| `otimizacao_torcao.png`, `torcao_otimizada.json` | otimização de torção |
| `evolucao_1_asa_isolada.gif`, `evolucao_3_com_estol.gif` | caminho da otimização |
| `modelo_arrasto_cache.json` | cache dos modelos de arrasto, refeito se o `aft.avl` mudar |

## Torção da asa

A asa sem torção estola em η = 0,842, dentro do aileron (η de 0,56 a 0,90),
logo ela não atende a FAR 25.203 e alguma torção é obrigatória. A pergunta de
projeto não é quanto a torção ganha em arrasto, é qual o menor preço da
conformidade.

A CL fixo a distribuição de sustentação é afim nas torções de estação e o
arrasto induzido é exatamente quadrático nelas, de modo que o modelo usado na
otimização é exato e construído com rodadas do AVL. Todo ótimo é conferido
fora do modelo.

### Por que a torção não é livre por estação

Com uma torção livre por estação o ótimo do modelo é serrilhado e salta 9,5
graus entre estações vizinhas. Isso não é ruído do otimizador: resolvendo o
mesmo problema como QP convexo por pontos interiores sai exatamente a mesma
solução, e a Hessiana é definida positiva com razão de condição 59. É o ótimo
verdadeiro do modelo, que usa um pico local de incidência como tira de estol.
O VLM aceita isso; a asa real não.

A torção é portanto uma **spline PCHIP monotônica**: Hermite cúbica por partes
com as inclinações de Fritsch-Carlson, que preserva monotonicidade e não
ultrapassa os valores de controle. Uma cúbica natural pelos mesmos pontos
chega a subir, ou seja faria a barriga que se quer evitar.

A monotonicidade não é imposta por restrição e sim pela variável de projeto:
otimizam-se os decrementos entre nós, com a raiz em zero por gauge,

    t_1 = 0,   t_k = -(s_1 + ... + s_{k-1}),   s_k >= 0

que é não crescente para qualquer s no octante positivo. A raiz fica sendo a
estação de maior incidência, que é o washout clássico. Para evitar a fronteira
da caixa, onde o Fminbox diverge, usa-se s = DEC_MAX/(1+exp(-u)) e o problema
vira irrestrito em u.

### Onde pôr os nós

Este é o parâmetro que mais importa, e só aparece com a restrição ativa. Para
arrasto puro as nove configurações testadas ficam dentro de 1 count umas das
outras; com a exigência de estol o espalhamento é de 79 counts.

| nós | sem restrição | com estol 0,5° |
|---|---|---|
| 0,00 0,25 0,50 0,75 1,00 | -15,63 | +81,53 |
| 0,00 0,22 0,40 0,56 0,75 1,00 | -15,70 | +23,65 |
| 0,00 0,48 0,56 0,75 1,00 | -15,36 | +2,60 |
| **0,00 0,30 0,48 0,56 0,70 0,85 1,00** | **-15,50** | **+2,59** |

O que decide é ter um nó em 0,48, logo antes da raiz do aileron em 0,56. Com
ele a curva fica plana até 0,48 e vira depressa depois, sem vazar washout para
a região que precisa estolar primeiro. Sem ele a mesma exigência custa vinte
counts a mais. A lição vale para além deste caso: estudo de discretização
precisa ser feito com as restrições ativas.

### Resultado

| etapa | CDff | counts | Oswald | estol η |
|---|---|---|---|---|
| sem torção | 0,010982 | - | 0,7377 | 0,842 |
| 1: asa isolada, irrestrito | 0,008251 | -12,00 | 1,0035 | - |
| 2: completa trimada, irrestrito | 0,009411 | -15,71 | 0,8614 | 0,842 |
| 3: projeto, margem 0,5° | 0,011160 | +1,78 | 0,7279 | 0,520 |

A etapa 1 valida o método: a spline monotônica chega à carga elíptica, com
Oswald 1,0035 e arrasto 0,33 count abaixo do piso plano CL²/(π AR), o que é
esperado porque a asa tem 6 graus de diedro e o mínimo de Munk de asa não
plana fica abaixo do elíptico plano. Contra a torção livre por estação a
spline custa apenas 0,97 count, de modo que a liberdade serrilhada valia menos
de um count.

Preço da conformidade, com a parametrização adotada:

| margem | counts | estol η | torção total | maior salto | taxa |
|---|---|---|---|---|---|
| 0,0° | -9,35 | 0,520 | 6,5° | 2,4° | 1,01 °/m |
| 0,5° | +1,78 | 0,520 | 6,6° | 5,1° | 2,13 °/m |
| 1,0° | +31,89 | 0,439 | 8,7° | 6,0° | 2,49 °/m |
| 1,5° | +61,44 | 0,243 | 10,1° | 5,5° | 2,30 °/m |
| 2,0° | +80,05 | 0,243 | 10,8° | 4,8° | 2,01 °/m |
| 3,0° | inviável | | | | |

A margem de estol é a variável cara, a cerca de 45 counts por grau. Isso vem
da distribuição de clmax ser quase uniforme ao longo da envergadura (1,774 na
raiz, 1,7985 no meio, 1,7338 na ponta), o que deixa a asa sem preferência
natural por onde estolar. Washout linear e torção quadrática são inviáveis em
qualquer margem.

A margem nula é melhor nos dois critérios ao mesmo tempo, arrasto e taxa de
torção, de modo que a margem de 0,5 grau compra apenas robustez do critério de
estol, que é linearizado e apoiado num AVL invíscido comparado a clmax de
XFoil. Torção adotada na etapa 3:

| η | 0,000 | 0,101 | 0,220 | 0,320 | 0,398 | 0,480 | 0,560 | 0,700 | 0,820 | 0,900 | 1,000 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| torção [°] | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | -5,13 | -6,55 | -6,55 | -6,55 | -6,55 |

Ressalva: a mudança de 5,13 graus entre as estações 0,48 e 0,56 dá 2,13 °/m,
bem acima dos cerca de 0,3 °/m de um transporte. A solução é monotônica e
lisa, mas concentra a variação de incidência num trecho curto junto à raiz do
aileron, que é onde a restrição age.

Esta solução coincide com a obtida antes por um caminho independente, um QP
convexo sobre as torções de estação com restrição de monotonicidade, o que dá
confiança de que é o ótimo global e não um mínimo local.

### Verificação

O `verifica_torcao.jl` confere o resultado por fora dos modelos que o
produziram. A parte de física checa monotonicidade e tamanho da torção, o erro
do modelo de arrasto em pontos que não o construíram, a coerência entre CDff,
CL e Oswald, e refaz a posição do estol por varredura direta de ângulo de
ataque no AVL em vez de extrapolar dos dois ângulos do modelo afim. A parte de
otimização testa as condições de KKT no espaço dos decrementos.

O teste de KKT é necessário por um motivo específico. A otimização roda numa
reparametrização logística, e ali ds/du tende a zero junto ao batente, de modo
que a estacionariedade em u não distingue um mínimo legítimo de uma parada
prematura: só o sinal de df/ds no espaço original separa os dois casos. Nas
três etapas os multiplicadores dos batentes ativos saem positivos, o que
confirma que apertar mais aqueles decrementos pioraria o arrasto.

Dois pontos que o teste revelou e que valem registro:

A restrição de estol escrita como margem do aileron é um mínimo de mínimos,
portanto não diferenciável. No ótimo da etapa 3 há duas faixas do aileron
empatadas, em η = 0,596 e η = 0,842, e por isso o gradiente por diferenças
finitas da margem é arbitrário. Com a restrição desagregada por faixa, que é a
forma correta, o resíduo de estacionariedade cai de 1,67 para zero e os dois
multiplicadores saem positivos, somando 36,5 counts por grau. Esse valor fica
entre as secantes da tabela de margens, 22,3 counts por grau entre 0 e 0,5 e
60,2 entre 0,5 e 1,0, como se espera da derivada de uma função convexa
crescente.

Nenhuma etapa encostou no batente superior de 6 graus por intervalo, e a
torção total ficou entre 6,5 e 8,8 graus contra um limite de 12. Ou seja,
nenhum resultado está sendo determinado por um limite arbitrário: os únicos
batentes ativos são os de monotonicidade.

### De onde vem o Oswald abaixo de 1 (e um alerta sobre as naceles)

A asa isolada da etapa 1 chega a e = 1,0035, mas a aeronave completa da etapa
2 fica em 0,8614 mesmo no ótimo. O `carga_completa.jl` separa a carga por
superfície e monta a aeronave peça por peça para achar a causa.

A explicação de manual seria arrasto de trimagem, com a empenagem carregando
para baixo e a asa tendo de carregar mais que o peso. **Não é o caso aqui.** No
ótimo da etapa 2 a empenagem carrega CL de -0,0020, praticamente nada, e a asa
carrega 0,4792 contra um CLff de 0,4998, ou seja MENOS que o total. Pelo fator
de trimagem o e subiria acima de 1.

A decomposição, com a torção congelada na da etapa 2 e o CL fixo:

| configuração | CDff | e | Δ counts |
|---|---|---|---|
| só a asa | 0,008379 | 0,9877 | - |
| asa + fuselagem | 0,008330 | 0,9713 | -0,48 |
| asa + fuselagem + naceles | 0,009293 | 0,8725 | **+9,63** |
| tudo, profundor em zero | 0,009299 | 0,8708 | +0,06 |
| tudo, trimado | 0,009411 | 0,8614 | +1,11 |

As naceles respondem por 9,63 dos 11,6 counts. A empenagem custa 1,17 no total
e a fuselagem chega a melhorar um pouco.

Isso pede atenção. A nacele está modelada como SUPERFÍCIE sustentadora em
forma de anel, com cerca de 91 m² de malha cada, quase metade da área da asa,
e carrega 4,8% da sustentação da aeronave. No método de malha de vórtices um
anel em ângulo de ataque gera circulação, esteira e portanto arrasto induzido;
uma nacele real tem escoamento passante e não se comporta assim. A
representação usual no AVL para nacele é BODY, que desloca o escoamento sem
sustentar, como já é feito com a fuselagem.

O efeito não fica só no valor absoluto: com a mesma torção e o mesmo CL, tirar
as naceles muda a carga local da asa em até 11% na estação onde elas ficam,
η = 0,36. Ou seja, elas influenciam a torção ótima e não apenas o arrasto
total. Conclusões comparativas entre etapas continuam válidas, porque as
naceles são as mesmas em todos os casos, mas o valor absoluto do arrasto
induzido e a forma fina da torção dependem dessa escolha de modelagem.
