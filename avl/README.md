# Lab 04 — análise aerodinâmica no AVL (equipe Ararinha)

Todo o Lab 04 é gerado por scripts Python a partir do `designTool`. Nenhum
arquivo do AVL é editado à mão: `fwd.avl` e `aft.avl` são saídas do
`gera_avl.py`.

## Requisitos

- Windows com o `avl.exe` 3.37 da disciplina, que já está nesta pasta.
- Python 3.13 com `numpy`, `scipy`, `matplotlib` e `pytest`
  (`requirements.txt` na raiz).
- O `designTool/` da disciplina na pasta-pai. Os scripts o importam de lá.

O AVL roda sempre com o diretório corrente nesta pasta e com caminhos
relativos, porque o executável não lida com o "º" do caminho absoluto. Os
scripts já fazem isso.

## Como rodar

De dentro de `avl/`, na ordem abaixo:

| # | Script | O que faz | Tempo aprox. |
|---|---|---|---|
| 1 | `python convergencia_malha.py` | convergência de malha e malha adotada (`resultados/malha_adotada.json`) | 10 min |
| 2 | `python verifica_eliptica.py` | verificação do otimizador: asa limpa → carga elíptica | 2 min |
| 3 | `python varredura_eh.py` | varredura da posição da EH e escolha do `Lc_h` | 7 min |
| 4 | `python otimizacao_torcao.py` | torção ótima e geração dos `fwd.avl`/`aft.avl` finais | 30 min |
| 5 | `python ponto_neutro.py` | ponto neutro montado componente a componente, designTool × AVL | 1 min |
| 6 | `python lab04_roteiro.py` | itens 1 a 8 do roteiro (`resultados/roteiro/`) | 15 min |
| 7 | `python derivadas_estabilidade.py` | seção 3: Tabelas 4, 6, 7 e 9 (`resultados/derivadas/`) | 1 min |
| 8 | `python estabilidade_direcional.py` | decomposição de Cnβ por componente | 1 min |
| 9 | `python sensibilidade_estabilidade.py` | correção da margem estática por Cht ou por recuo da asa | 3 min |

Testes: `python -m pytest tests -q` (41 testes).

Material complementar:

- `investigacao_torcao/`: verificação da otimização de torção, com o
  histórico de convergência, múltiplas partidas, margens de estol maiores,
  parametrizações suaves e o caso sem EH. `python figuras_investigacao.py`
  refaz as figuras e o `resumo.csv` a partir de `runs/`.
- `planilha_pn/`: planilha de ponto neutro da disciplina (Torenbeek)
  preenchida com a geometria final, na versão original recebida e na
  corrigida. Veja o README da pasta.

O `python.exe` da Microsoft Store é só um lançador e termina logo; o
processo real é o `python3.13.exe`. Para acompanhar execuções longas, rode
o interpretador real.

## Módulos

| arquivo | o que faz |
|---|---|
| `aeronave.py` | roda o designTool para um `Lc_h` e devolve a geometria, o ponto de projeto e os CGs, com as alturas reais |
| `gera_avl.py` | escreve o `.avl`: deslocamento em Z, torção, malha, winglet, nacele e fuselagem reamostrada |
| `avl_run.py`, `avl_saida.py` | rodam o AVL por stdin e leem `ft`, `fs`, `st` e `sb` |
| `analises.py` | `it` que zera o profundor no cruzeiro |
| `estol.py` | método da seção crítica, com o clmax dos perfis do Lab 03 |
| `sombra.py` | sombreamento da EH pela asa e leme encoberto, pela geometria real |

## Decisões de modelagem

- **Deslocamento em Z só no AVL**, por exigência do professor, para afastar
  os painéis: asa −1,20 m, EH e EV +1,85 m, nacele junto com a asa.
  Nenhuma superfície cruza a fuselagem (folga mínima de 0,15 m, na raiz da
  asa), e EH e EV mantêm entre as raízes a separação real de 1 m do
  designTool. Não é alteração de projeto.
- **Winglet** em superfície própria no COMPONENT 1, com toe zero.
- **Nacele** como anel sustentador em componente próprio. No COMPONENT 1 o
  CDff oscilava ±0,2 count com a malha da asa.
- **Fuselagem** com o contorno reamostrado (`fuselage_reamostrada.dat`). O
  arquivo original fazia a spline do AVL formar um laço na cauda.
- **Estol:** as faixas da raiz e da junção com o winglet ficam fora da
  seção crítica, porque são singularidades de quina do VLM.
- **Arrasto:** sempre `CD = CDp + CDff` (Trefftz).
- **Controles:** `d1` aileron, `d2` profundor, `d3` leme. A variável de
  projeto `1` é o `it`. A numeração difere do `b737mod.avl` do roteiro,
  que usa `d4` e a variável `2`.

## Resultados principais

| | |
|---|---|
| Malha adotada | asa 8×60, winglet 4, EH 8×20, EV 8×6, corpo 80 |
| Verificação do otimizador | e 0,8715 → 0,9994 (asa limpa) |
| `Lc_h` escolhido | 4,8237 (limite da fuselagem); W0 288.764 kgf; S_h 53,52 m² |
| Torção ótima | CDff 104,0 → 94,0 count; estol em η 0,493, α 13,2°; folga de sombreamento 3,9° |
| `it` que zera δe | −3,242° (CG dianteiro); −0,964° (CG traseiro) |
| CLmax (M 0,2) | 1,14 a 1,19 |
| MS com CG traseiro | AVL −2,9% (M 0,85), +4,1% (M 0,2); designTool +7,3% |
| PN com e sem naceles | AVL 36,3% / 42,6%; planilha (Torenbeek) 38,7% / 44,7%; designTool (sem naceles) 46,5% (M 0,85) |
| CG | não alterado; pelo AVL, o CG traseiro teria de avançar 0,2 m (MS = 0) a 0,5 m (MS = 5%) |
