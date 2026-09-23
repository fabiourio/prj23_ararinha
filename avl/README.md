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

**Alturas em Z, só no AVL.** Asa −1,20 m, EH +1,85 m, EV +0,85 m e naceles
acompanhando a asa, para tirar as superfícies de dentro da fuselagem (que
ocupa z de −2,348 a +3,732) e o VLM não gerar painéis cruzando o corpo. Não é
alteração de projeto: a aeronave do designTool mantém as alturas originais, e
o `dados_designtool.py` confere que o AVL difere dela exatamente por esses
valores.

## Geometria da v3

**EH recuada**, de `Lc_h` = 4,0 para 4,6. Com o volume de cauda `Cht` fixo,
braço maior dá EH menor. O ganho é monótono, sem ótimo intermediário, de modo
que quem decide é a restrição geométrica:

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

A asa sem torção estola primeiro em η = 0,923, dentro do aileron (η de 0,56 a
0,90) e junto à ponta, que o winglet carrega mais. Ela não atende a FAR 25.203,
que exige comando de rolamento eficaz até e durante o estol, de modo que alguma
torção é obrigatória. A pergunta de projeto não é quanto a torção ganha em
arrasto, é qual o menor preço da conformidade.

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

### Restrição de estol: posição, com desempate

A exigência é de posição: **a primeira faixa a estolar tem de estar para dentro
de η_lim = 0,56**, a raiz do aileron, que é o que a FAR pede. Não há margem
angular de projeto.

Há, porém, um desempate de 0,2 grau, e ele não é arbitrário. Sem nenhuma folga
o ótimo sempre cai num empate: a faixa interna e uma faixa do aileron estolam
no mesmo ângulo, porque é ali que o arrasto é menor, e a exigência vale só no
papel. Medido direto no AVL, um projeto assim estolou primeiro em η = 0,842,
dentro do aileron. A causa é o modelo afim de estol, que extrapola de dois
ângulos e é ligeiramente otimista perto do estol: nos casos medidos, ele
superestimou a folga entre 0,11 e 0,19 grau. O desempate é esse erro, para que a
posição do estol valha também na medição direta, que a verificação confere. No
projeto final a folga pelo modelo é 0,200° e a medida direto é +0,094°.

### Onde pôr os nós

O parâmetro que mais importa, e ele só aparece com a restrição ativa. Para
arrasto puro as nove configurações ficam dentro de 1,2 count umas das outras;
com a exigência de estol o espalhamento passa de 70 counts (valores do modelo):

| nós | sem restrição | com estol antes de 0,56 |
|---|---|---|
| 0,00 0,25 0,50 0,75 1,00 | −12,37 | +16,63 |
| 0,00 0,22 0,40 0,56 0,75 1,00 | −12,52 | +8,28 |
| 0,00 0,48 0,56 0,75 1,00 | −12,11 | +0,91 |
| **0,00 0,30 0,48 0,56 0,70 0,85 1,00** | **−12,26** | **+0,91** |

O que decide é ter um nó em 0,48, logo antes da raiz do aileron em 0,56: a curva
fica plana até ali e vira depressa depois, sem vazar washout para a região que
precisa estolar primeiro.

### Resultado

| caso | CDff | counts | Oswald | estol começa em η | CL no estol |
|---|---|---|---|---|---|
| sem torção | 0,010030 | - | 0,7983 | 0,923 | 0,954 |
| 2: completa compensada, sem restrição | 0,008784 | −12,46 | 0,9180 | 0,842 | 1,073 |
| **3: projeto, estol antes de 0,56** | **0,010045** | **+0,15** | **0,8064** | **0,520** | **1,217** |

Torção adotada:

| η | 0,000 | 0,101 | 0,220 | 0,320 | 0,398 | 0,480 | 0,560 | 0,700 | 0,820 | 0,900 | 1,000 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| torção [°] | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | −2,95 | −5,81 | −5,81 | −5,81 | −5,81 |

A aeronave fica conforme por praticamente o mesmo arrasto da asa sem torção
(+0,15 count), e o CL em que o estol começa sobe de 0,954 para 1,217, 28% a
mais, logo abaixo do CLmax de 1,333 que o designTool estima para a aeronave
limpa. O maior salto entre estações vizinhas é 2,95 graus (1,22 °/m), acima dos
cerca de 0,3 °/m de um transporte, porque a variação de incidência se concentra
junto à raiz do aileron, onde a restrição age.

Custo de exigir o início do estol mais para dentro:

| estol antes de | CDff | counts | Oswald |
|---|---|---|---|
| **0,56** | **0,010045** | **+0,15** | **0,806** |
| 0,50 e 0,45 | 0,013411 | +33,80 | 0,609 |
| 0,40 e 0,35 | 0,013564 | +35,34 | 0,603 |

O salto entre 0,56 e 0,50 é grande porque a malha só enxerga os centros das
faixas (0,439, 0,520, 0,596...): exigir 0,50 tira da região interna a faixa em
0,520, que é justamente a que estola primeiro no projeto. Por isso 0,45 e 0,50
dão o mesmo resultado.

### Por que o Oswald fica em 0,81

A cascata abaixo, feita pelo `carga_completa.jl` com a torção da etapa 2 e o CL
fixo, atribui cada queda a uma peça da aeronave:

| configuração | e | Δe |
|---|---|---|
| só a asa, sem winglet | 0,982 | - |
| + winglet | 1,078 | +0,096 |
| + fuselagem | 1,058 | −0,020 |
| **+ naceles** | **0,943** | **−0,115** |
| + empenagens, profundor em zero | 0,932 | −0,011 |
| compensado, ótimo da etapa 2 | 0,918 | −0,014 |
| **+ restrição de estol (projeto)** | **0,806** | **−0,112** |

A otimização sem restrição entrega 0,918, e a asa isolada com winglet chega a
1,08. As duas maiores quedas são as naceles e a restrição de estol, cada uma de
cerca de 0,11.

**Alerta sobre as naceles.** Elas estão modeladas como SUPERFÍCIE sustentadora
em forma de anel, com cerca de 91 m² de malha cada, quase metade da área da
asa. No VLM um anel em ângulo de ataque gera circulação, esteira e arrasto
induzido; uma nacele real tem escoamento passante e não se comporta assim. A
representação usual no AVL para nacele é BODY, que desloca o escoamento sem
sustentar. O efeito não fica só no valor absoluto: com a mesma torção e o mesmo
CL, tirar as naceles muda a carga local da asa em até 11% em η = 0,36.
Comparações entre casos continuam válidas, porque as naceles são as mesmas em
todos, mas o valor absoluto do Oswald depende dessa escolha de modelagem.

### Verificação

O `verifica_torcao.jl` confere o resultado por fora dos modelos que o
produziram, e passa em todas as 20 checagens.

- **Física**: monotonicidade e tamanho da torção, coerência entre CDff, CL e
  Oswald, e a posição do estol refeita por varredura direta de ângulo de ataque
  no AVL, sem extrapolar do modelo afim. Na medição direta o estol do projeto
  começa em η = 0,520, antes da raiz do aileron, com folga de +0,094°.
- **KKT no espaço dos decrementos**, necessário porque na reparametrização
  logística ds/du tende a zero junto ao batente, e a estacionariedade em u não
  distingue um mínimo legítimo de uma parada prematura. O ótimo é um vértice:
  seis variáveis, quatro batentes de monotonicidade e duas faixas de estol
  ativas (η = 0,596 e 0,842), com multiplicadores de +19,6 e +6,7 counts por
  grau nas faixas de estol e resíduo zero. Um dos batentes (s5) tem
  multiplicador de −0,016 count por grau, zero para fins práticos: é o
  vestígio da parada perto do batente que a reparametrização causa, e o ganho
  possível ao sair dele fica abaixo de décimos de count.
- **O erro do modelo de arrasto não mexe no ótimo.** Remontando o modelo
  centrado no próprio ótimo, onde ele passa a ser exato, e resolvendo de novo,
  os decrementos saem idênticos até a terceira casa. O vértice é fixado pelas
  restrições de estol, e o modelo de arrasto só escolhe qual vértice.
- **Caminho independente**: refeita a partir da asa sem torção, sem nenhuma
  semente da busca global, a otimização chega ao mesmo ponto (+0,00 count)
  depois de cerca de 2000 iterações. É esse caminho que o
  `evolucao_3_com_estol.gif` mostra.
- **Otimizar sem restrição e corrigir depois não ganha.** A alternativa óbvia
  seria otimizar só o arrasto e dar washout na ponta até o estol sair do
  aileron. Todo projeto assim é um ponto viável, e o ótimo restrito não pode
  perder para nenhum deles. Nenhum perde:

| etapa 2 + washout manual | washout necessário | contra o otimizador |
|---|---|---|
| rampa linear desde η = 0,48 (melhor rampa) | 11,6° | +48,5 counts |
| degrau até η = 0,56 | 5,2° | +25,8 counts |
| **degrau até η = 0,60 (melhor manual)** | **4,7°** | **+14,1 counts** |
| degrau até η = 0,70 | 5,4° | +19,3 counts |

  A rampa começando antes do aileron descarrega também as faixas que precisam
  estolar primeiro e pede 11° a 18° de washout. O otimizador deixa a asa plana
  até 0,48 e concentra o washout na raiz do aileron, e mesmo o degrau manual com
  essa forma perde 14 counts porque parte da torção da etapa 2 em vez de
  escolher a torção inteira sabendo da restrição.

Nenhuma etapa encostou no batente superior de 6 graus por intervalo, e a
torção total ficou abaixo do limite de 12 graus: os únicos batentes ativos são
os de monotonicidade.
