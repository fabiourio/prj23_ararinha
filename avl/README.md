# Lab 04 -- análise aerodinâmica no AVL

## Como rodar

Tudo roda de dentro desta pasta, nesta ordem:

```
python dados_designtool.py      # ponte com o designTool + verificação da geometria
julia  convergencia_malha.jl    # define e justifica a malha
julia  lab04.jl                 # análises do roteiro
julia  otimizacao_torcao.jl     # otimização de torção (estudo próprio)
julia  suficiencia_nos.jl       # quantos nós de torção são necessários
```

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
