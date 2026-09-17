# Lab 04 -- análise aerodinâmica no AVL

## Como rodar

Tudo roda de dentro desta pasta, nesta ordem:

```
python dados_designtool.py      # ponte com o designTool + verificação da geometria
julia  convergencia_malha.jl    # define e justifica a malha
julia  lab04.jl                 # análises do roteiro
julia  refina_asa.py            # refina a asa para 11 seções (só se mudar a geometria)
julia  otimizacao_torcao.jl     # otimização de torção (estudo próprio)
julia  suficiencia_nos.jl       # quantos nós de torção são necessários
julia  torcao_monotonica.jl     # que vínculo a torção precisa obedecer
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
| `torcao_decisao.png`, `torcao_recomendada.txt` | escolha do vínculo da torção |
| `custo_conformidade.png`, `custo_conformidade.txt` | preço da margem de estol |
| `torcao_monotonica.png`, `candidatos_torcao.json` | comparação das parametrizações |

## Torção da asa

A asa sem torção estola em η = 0,842, dentro do aileron (η de 0,56 a 0,90),
logo ela não atende a FAR 25.203 e alguma torção é obrigatória. A pergunta de
projeto não é quanto a torção ganha em arrasto, é qual o menor preço da
conformidade.

A CL fixo a distribuição de sustentação é afim na torção, de modo que o
arrasto induzido é exatamente quadrático e o cl de cada faixa exatamente afim.
Cada caso é portanto um QP de restrições lineares, resolvido por pontos
interiores. A Hessiana é definida positiva (autovalores de 1,12e-4 a 1,89e-6,
razão 59), então o problema é convexo e o ótimo é único: não há ambiguidade de
mínimo local.

Com as torções livres o ótimo é serrilhado e salta 9,5 graus entre estações
vizinhas. Isso não é ruído do otimizador, é o ótimo verdadeiro do modelo: ele
usa um pico local de incidência como tira de estol, o que o VLM aceita e a asa
real não. O vínculo que falta não é monotonicidade e sim limite de salto.

Preço da conformidade, em counts de CDff sobre a asa sem torção:

| margem | livre | suave 3°/est | monotônica | cúbica monótona |
|---|---|---|---|---|
| 0,0° | -18,8 | -16,2 | -10,7 | +3,5 |
| 0,5° | -16,0 | -1,4 | +1,8 | inviável |
| 1,0° | -8,4 | +12,0 | +17,7 | inviável |
| 1,5° | -3,3 | +22,3 | +34,4 | inviável |
| 2,0° | +3,9 | +38,3 | +56,6 | inviável |

A margem de estol é a variável cara, não a forma da torção. Isso vem da nossa
distribuição de clmax ser quase uniforme ao longo da envergadura (1,774 na
raiz, 1,7985 no meio, 1,7338 na ponta), o que deixa a asa sem preferência
natural por onde estolar. Washout linear e torção quadrática são inviáveis em
qualquer margem.

Configuração recomendada, monotônica com salto máximo de 2 graus por estação e
margem nula, em `resultados/torcao_recomendada.txt`:

| η | 0,000 | 0,101 | 0,220 | 0,320 | 0,398 | 0,480 | 0,560 | 0,700 | 0,820 | 0,900 | 1,000 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| torção [°] | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | -0,98 | -2,98 | -4,98 | -5,53 | -5,82 | -5,82 |

CDff 0,010126 (-8,56 counts), Oswald 0,8014, estol em η = 0,520, faixa total de
torção 5,8 graus. A alternativa com meio grau de margem de estol custa +1,78
counts, ou seja arrasto neutro.
