# Lab 04 -- análise aerodinâmica no AVL

## Como rodar

Tudo roda de dentro desta pasta, nesta ordem:

```
python dados_designtool.py      # ponte com o designTool + verificação da geometria
julia  convergencia_malha.jl    # define e justifica a malha
julia  lab04.jl                 # análises do roteiro
julia  otimizacao_torcao.jl     # otimização de torção (estudo próprio)
julia  verifica_torcao.jl       # verificação independente do resultado acima
julia  carga_completa.jl        # carga por superfície e de onde vem o Oswald
julia  contexto_resultado.jl    # resultado contra o designTool e transportes reais
python refina_asa.py            # refina a asa para 11 seções (só se mudar a geometria)
```

O `dados_designtool.py` precisa ser rodado primeiro: ele gera o
`dados_designtool.json`, de onde os scripts Julia leem o ponto de projeto e os
dados de referência. É o único lugar que chama o designTool.

O `otimizacao_torcao.jl` guarda o modelo de arrasto em
`resultados/modelo_arrasto_cache.json` e só o reconstrói se o `aft.avl` for mais
novo que o cache. Construí-lo custa 66 rodadas do AVL; rodar os casos em cima
dele custa segundos.

## Modelo

| arquivo | o que é |
|---|---|
| `aft.avl` / `fwd.avl` | aeronave com CG traseiro e dianteiro. Diferem só na linha `Xref` |
| `airfoils/` | perfis otimizados no Lab 03 (família roteiro): raiz, centro, ponta |
| `fuselage_nondim.dat` | contorno da fuselagem (`BODY` do AVL) |
| `737.avl` | modelo de referência da disciplina, para consulta |
| `avl.exe` | executável da disciplina (versão 3.37) |

**Malha adotada:** asa 8 × 20, winglet 8 × 8, empenagem horizontal 16 × 10,
empenagem vertical 8 × 6, corpo com 20 painéis. A envergadura da asa foi
definida pelo arrasto e a corda pelos momentos e pela posição efetiva da
charneira do profundor (ver `convergencia_malha.jl`). O winglet tem malha
própria, escolhida no platô: a partir de 6 faixas a variação do CDff fica
abaixo de 0,03 count.

## Geometria v3

A v3 do designTool (`my_airplane`) traz três mudanças, e o `.avl` descreve
exatamente a mesma aeronave: o `dados_designtool.py` confere cada cota, e o
deslocamento entre os dois modelos agora é nulo.

**Alturas.** Fuselagem no mesmo lugar, asa descendo 1,20 m (`zr_w` = −2,5, raiz
abaixo da base da fuselagem, que ocupa z de −2,348 a +3,732), EH e EV subindo
para `zr_h` = `zr_v` = 3,85 (raiz acima do topo). As naceles acompanham a asa
(`z_n` = −4,2). Até a v2 essas alturas eram um deslocamento só do AVL, para
tirar as superfícies de dentro da fuselagem; agora são da aeronave.

**EH recuada**, de `Lc_h` = 4,0 para 4,6. Com o volume de cauda `Cht` fixo, braço
maior dá EH menor. O ganho é monótono, sem ótimo intermediário, de modo que
quem decide é a restrição geométrica:

| Lc_h | S_h [m²] | folga até o fim da fuselagem | W0 [kgf] | CDff compensado | profundor |
|---|---|---|---|---|---|
| 4,0 | 64,5 | 6,10 m | 291 292 | referência | −6,23° |
| 4,4 | 58,7 | 3,38 m | 289 948 | −0,47 count | −6,49° |
| **4,6** | **56,1** | **2,02 m** | **289 362** | **−0,63 count** | **−6,79°** |
| 4,8 | 53,8 | 0,66 m | 288 826 | −1,15 count | −7,78° |

Em 4,6 a raiz da EH termina em x = 59,4 m, onde o cone de cauda ainda tem 2,5 m
de altura. Ir a 4,8 renderia mais 540 kgf e meio count, mas encosta no fim do
cone e leva o profundor de compensação a −7,8°, consumindo autoridade de
controle.

**Winglet**, na geometria que o designTool assume (Raymer, Fig. 7.34): vertical,
altura e corda de raiz iguais à corda da ponta (2,3736 m), afilamento 0,21 e
bordo de fuga reto. No AVL ele é uma superfície própria, e não uma seção a mais
da asa, por dois motivos medidos:

- numa superfície contínua a seção da ponta é compartilhada, e o `Ainc` dela,
  que é a torção da asa, vira ângulo de convergência (toe) do winglet. A
  otimização ganhava assim uma variável escondida: cada grau de washout na
  ponta girava o winglet junto;
- com `COMPONENT 1` o winglet fica no mesmo componente da asa e a junção é
  tratada como folha de vórtices contínua. Em componentes separados o AVL usa
  núcleo finito entre os dois e a quina perde metade do efeito.

| modelagem do winglet (asa sem torção) | efeito no CDff |
|---|---|
| seção a mais da asa (média entre malhas de 20 a 48 faixas) | −9,5 ± 1,3 counts |
| superfície própria, componente próprio | −5,5 counts |
| **superfície própria, COMPONENT 1** | **−8,9 counts** |

O valor adotado cai dentro da banda da versão contínua, e a versão em
componente separado é a discrepante.

**Pendência de trem de pouso, a decidir pelo grupo.** Com a asa 1,20 m mais
baixa e a fuselagem no lugar, as naceles chegam 1,20 m mais perto do chão: a
folga ao solo do designTool passa de +1,02 m para **−0,18 m**. Voltar a nacele
para cima não resolve, porque ela atravessaria a asa. Alongar o trem 1,20 m
recupera a folga e melhora o ângulo de tailstrike (10,1° para 14,0°), mas o
tipback cai de 15,4° para 12,9° e fica abaixo do tailstrike, contra a regra
clássica. Nada disso afeta a aerodinâmica deste lab: peso e CG praticamente
não dependem do comprimento do trem no designTool.

## Numeração dos controles

A numeração do AVL segue a ordem de aparição no arquivo e **difere** da do
roteiro do professor (que usa o `b737mod.avl`):

| roteiro | aqui | o que é |
|---|---|---|
| `d4 pm 0.0` | `d2 pm 0.0` | compensação pelo profundor |
| menu `de`: `2 <valor>` | menu `de`: `1 <valor>` | incidência da empenagem |

Controles: `d1` aileron, `d2` profundor, `d3` leme. Variável de projeto:
`1` = incidência da empenagem (`it`).

## Ponto de projeto

O mesmo do Lab 03, o peso médio de cruzeiro, para que o CL da aeronave coincida
com o CL para o qual os perfis foram otimizados: M = 0,85, h = 10.668 m,
W = 229.669 kgf, **CL = 0,5053**. O peso fica fixo de propósito: a EH mais leve
reduz o W0 em 0,7%, e manter o peso deixa o CL de projeto colado no das seções.

Nos comandos do AVL: `m` → `mn 0.85` (ou `mn 0.2` para baixa velocidade) e
`a c 0.5053`.

## Resultados

Tudo que os scripts geram vai para `resultados/`:

| arquivo | de onde vem |
|---|---|
| `conv_corda.png`, `conv_envergadura.png`, `malha_adotada.txt` | convergência de malha |
| `ponto_projeto_carga.png`, `ponto_de_projeto.txt` | ponto de projeto |
| `otimizacao_torcao.png`, `torcao_otimizada.json` | otimização de torção |
| `evolucao_3_com_estol.gif` | caminho da otimização, a partir da asa sem torção |
| `mapa_restricao_estol.png` | arrasto e restrição de estol no plano das variáveis livres |
| `carga_completa.png` | carga de asa e EH no mesmo eixo de envergadura |
| `modelo_arrasto_cache.json` | cache do modelo de arrasto, refeito se o `aft.avl` mudar |

## Torção da asa

A asa sem torção estola em η = 0,923, dentro do aileron (η de 0,56 a 0,90) e
junto à ponta, que o winglet carrega mais. Ela não atende a FAR 25.203, que exige
comando de rolamento eficaz até e durante o estol, de modo que alguma torção é
obrigatória. A pergunta de projeto não é quanto a torção ganha em arrasto, é
qual o menor preço da conformidade.

A CL fixo a distribuição de sustentação é afim nas torções de estação e o
arrasto induzido é quadrático nelas. O modelo usado na otimização é montado com
rodadas do AVL nesse espaço, e todo ótimo é conferido fora dele.

### Parametrização: spline monotônica

Com uma torção livre por estação o ótimo do modelo é serrilhado, com saltos de
9,5 graus entre estações vizinhas. Não é ruído: resolvido como QP convexo por
pontos interiores sai a mesma solução. O VLM aceita um pico local de incidência
funcionando como tira de estol; a asa real não.

A torção é por isso uma **spline PCHIP monotônica** (Hermite cúbica com as
inclinações de Fritsch-Carlson), que preserva monotonicidade e não ultrapassa os
valores de controle. Uma cúbica natural pelos mesmos pontos chega a subir, isto
é, faria a barriga que se quer evitar. A monotonicidade vem da própria variável
de projeto, os decrementos entre nós com a raiz em zero por gauge,

    t_1 = 0,   t_k = -(s_1 + ... + s_{k-1}),   s_k >= 0

de modo que a raiz é a estação de maior incidência, o washout clássico. Para
evitar a fronteira da caixa, onde o Fminbox diverge, usa-se
s = DEC_MAX/(1+exp(-u)) e o problema vira irrestrito em u.

A validação com a asa isolada foi feita na v2 e não é repetida: a spline chegou
à carga elíptica com Oswald 1,0035 e arrasto 0,33 count abaixo do piso plano,
coerente com os 6 graus de diedro, e custou só 0,97 count contra a torção livre.

### Onde pôr os nós

O parâmetro que mais importa, e ele só aparece com a restrição ativa. Para
arrasto puro as nove configurações ficam dentro de 1,2 count umas das outras;
com a exigência de estol o espalhamento passa de 40 counts.

| nós | sem restrição | com estol 0,5° |
|---|---|---|
| 0,00 0,25 0,50 0,75 1,00 | −12,37 | +53,45 |
| 0,00 0,22 0,40 0,56 0,75 1,00 | −12,52 | +32,52 |
| 0,00 0,48 0,56 0,75 1,00 | −12,11 | +11,03 |
| **0,00 0,30 0,48 0,56 0,70 0,85 1,00** | **−12,26** | **+10,97** |

O que decide é ter um nó em 0,48, logo antes da raiz do aileron em 0,56: a curva
fica plana até ali e vira depressa depois, sem vazar washout para a região que
precisa estolar primeiro. Estudo de discretização precisa ser feito com as
restrições ativas.

### Resultado

| caso | CDff | counts | Oswald | estol η |
|---|---|---|---|---|
| sem torção | 0,010030 | - | 0,7983 | 0,923 |
| 2: completa compensada, sem restrição | 0,008784 | −12,46 | 0,9180 | 0,842 |
| **3: projeto, margem 0,5°** | **0,011007** | **+9,77** | **0,7379** | **0,520** |

Preço da conformidade:

| margem | counts | Oswald | estol η | torção total | maior salto | taxa |
|---|---|---|---|---|---|---|
| 0,0° | −1,66 | 0,8215 | 0,520 | 8,5° | 2,4° | 0,92 °/m |
| 0,5° | +9,77 | 0,7379 | 0,520 | 6,8° | 4,8° | 2,01 °/m |
| 1,0° | +45,21 | 0,5623 | 0,439 | 9,4° | 6,0° | 2,49 °/m |
| 1,5° | +72,46 | 0,4756 | 0,243 | 10,5° | 6,0° | 2,49 °/m |
| 2,0° | +90,62 | 0,4316 | 0,243 | 11,1° | 6,0° | 2,49 °/m |
| 3,0° | inviável | | | | | |

Torção adotada:

| η | 0,000 | 0,101 | 0,220 | 0,320 | 0,398 | 0,480 | 0,560 | 0,700 | 0,820 | 0,900 | 1,000 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| torção [°] | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | −4,84 | −6,81 | −6,81 | −6,81 | −6,81 |

Contra a v2 a aeronave de projeto fica 1,5 count melhor em CDff (0,011160 para
0,011007), mas o winglet encarece a restrição: sem torção a v3 já é 9,5 counts
melhor, e conformidade custa +9,8 counts contra +1,8 na v2. É o compromisso
clássico entre winglet e estol de ponta: o winglet carrega a ponta, o estol
sem torção vai para η = 0,923, e a torção precisa de mais trabalho para trazê-lo
para dentro do aileron. A margem de estol é a variável cara porque a
distribuição de clmax é quase uniforme (1,774 na raiz, 1,7985 no meio, 1,7338 na
ponta): a asa não tem preferência natural por onde estolar.

**O Oswald de 0,74 no projeto é baixo, e é o esperado.** Ele fica abaixo até do
da asa sem torção porque toda a região de arrasto baixo é proibida pela
restrição: nela o estol começaria dentro do aileron. O `mapa_restricao_estol.png`
mostra isso no plano das duas variáveis livres (s3 e s4, com as outras no
batente zero). As curvas de nível do arrasto descem para a região inviável, e o
ótimo é o bico da fronteira, onde as duas faixas críticas do aileron, em
η = 0,596 e 0,842, se tornam ativas juntas.

Ressalvas. A mudança de 4,84 graus entre as estações 0,48 e 0,56 dá 2,01 °/m,
bem acima dos cerca de 0,3 °/m de um transporte. E a margem medida direto no
AVL, por varredura de ângulo de ataque, é 0,39°, e não os 0,50° que o modelo
afim de estol prevê: ele extrapola de dois ângulos e é ligeiramente otimista
perto do estol. Para 0,5° medidos, basta pedir 0,6° ao otimizador.

### Verificação

O `verifica_torcao.jl` confere o resultado por fora dos modelos que o
produziram, e passa em todas as 19 checagens.

- **Física**: monotonicidade e tamanho da torção, coerência entre CDff, CL e
  Oswald, e a posição do estol refeita por varredura direta de ângulo de ataque
  no AVL, sem extrapolar do modelo afim. O estol direto bate com o previsto nas
  duas etapas, em η = 0,842 e 0,520.
- **KKT no espaço dos decrementos**, necessário porque na reparametrização
  logística ds/du tende a zero junto ao batente, e a estacionariedade em u não
  distingue um mínimo legítimo de uma parada prematura. Na etapa 3 o ótimo é
  um vértice: seis variáveis, quatro batentes de monotonicidade e duas faixas
  de estol ativas, todos os multiplicadores positivos (+32,7 e +8,0 counts por
  grau nas faixas de estol) e resíduo zero. A restrição de estol escrita como
  margem do aileron é um mínimo de mínimos e não é diferenciável no empate
  entre as duas faixas, por isso o KKT é testado com ela desagregada por faixa.
- **O erro do modelo de arrasto não mexe no ótimo.** O modelo é montado em
  torno da asa sem torção e erra 1,2 count no ótimo, que tem quase 7 graus de
  washout. Remontando o modelo centrado no próprio ótimo, onde ele passa a ser
  exato, e resolvendo de novo, os decrementos saem idênticos até a terceira
  casa. Faz sentido: o vértice é fixado pelas restrições de estol, e o modelo de
  arrasto só escolhe qual vértice.
- **Caminho independente**: refeita a partir da asa sem torção, sem nenhuma
  semente da busca global, a otimização chega ao mesmo ponto (+0,00 count)
  depois de cerca de 2000 iterações. É esse caminho que o
  `evolucao_3_com_estol.gif` mostra. Num teste avulso, fora deste
  repositório, um CMA-ES livre de derivadas também não achou nada melhor em
  cinco sementes: quatro pararam entre 0,3 e 2,0 counts acima, e uma num
  mínimo local 35 counts acima.

Nenhuma etapa encostou no batente superior de 6 graus por intervalo, e a
torção total ficou abaixo do limite de 12 graus: os únicos batentes ativos são
os de monotonicidade.

### De onde vem o Oswald abaixo de 1 (e um alerta sobre as naceles)

O `carga_completa.jl` separa a carga por superfície e monta a aeronave peça por
peça. A explicação de manual para o e abaixo de 1 seria arrasto de compensação,
com a empenagem carregando para baixo e a asa tendo de carregar mais que o
peso. **Não é o caso aqui**: no ótimo da etapa 2 a empenagem está quase
descarregada.

A decomposição, com a torção congelada na da etapa 2 e o CL fixo:

| configuração | CDff | e | Δ counts |
|---|---|---|---|
| só a asa, sem winglet | 0,008427 | 0,9822 | - |
| asa com winglet | 0,007646 | 1,0775 | −7,81 |
| asa + fuselagem | 0,007608 | 1,0584 | −0,38 |
| asa + fuselagem + naceles | 0,008559 | 0,9429 | **+9,51** |
| tudo, profundor em zero | 0,008671 | 0,9316 | +1,13 |
| tudo, compensado | 0,008784 | 0,9180 | +1,13 |

As naceles comem mais do que o winglet ganha. Elas estão modeladas como
SUPERFÍCIE sustentadora em forma de anel, com cerca de 91 m² de malha cada,
quase metade da área da asa. No VLM um anel em ângulo de ataque gera
circulação, esteira e arrasto induzido; uma nacele real tem escoamento passante
e não se comporta assim. A representação usual no AVL para nacele é BODY, que
desloca o escoamento sem sustentar. O efeito não fica só no valor absoluto: com
a mesma torção e o mesmo CL, tirar as naceles muda a carga local da asa em até
11% em η = 0,36. Comparações entre casos continuam válidas, porque as naceles
são as mesmas em todos, mas o valor absoluto do induzido e a forma fina da
torção dependem dessa escolha de modelagem.
