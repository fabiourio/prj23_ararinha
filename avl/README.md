# Lab 04 — AVL (equipe Ararinha)

Os arquivos `fwd.avl` e `aft.avl` são gerados por `gera_avl.py` a partir do
`designTool`.

## Requisitos

- Windows, com o `avl.exe` 3.37 da disciplina (já nesta pasta).
- Python 3.13 com `numpy`, `scipy`, `matplotlib` e `pytest`
  (`requirements.txt` na raiz).
- `designTool/` na pasta-pai.

Rode os scripts de dentro de `avl/`.

## Como rodar

| # | Script | Resultado |
|---|---|---|
| 1 | `python convergencia_malha.py` | malha adotada |
| 2 | `python verifica_eliptica.py` | verificação do otimizador (asa isolada) |
| 3 | `python varredura_eh.py` | varredura da posição da EH |
| 4 | `python otimizacao_torcao.py` | torção ótima e `fwd.avl`/`aft.avl` finais |
| 5 | `python ponto_neutro.py` | ponto neutro por componente |
| 6 | `python lab04_roteiro.py` | itens 1 a 8 do roteiro |
| 7 | `python derivadas_estabilidade.py` | Tabelas 4, 6, 7 e 9 |
| 8 | `python estabilidade_direcional.py` | decomposição de Cnβ |

Testes: `python -m pytest tests -q`.

`investigacao_torcao/`: otimizações de verificação da torção (outras
partidas e restrições). `python figuras_investigacao.py` refaz as figuras.

`planilha_pn/`: planilha de ponto neutro da disciplina com a geometria final.
Na versão corrigida, `xr_h` passou a 56,185 m, os CGs a 15,45% e 39,21% da
CMA, e as células da nacele (C75 a C80) foram preenchidas pela Eq. E-41 de
Torenbeek.

## Convenções do modelo

- Só no AVL, a asa e as naceles foram baixadas 1,20 m e as empenagens subidas
  1,85 m, para que nenhuma superfície cruze a fuselagem.
- Controles: `d1` aileron, `d2` profundor, `d3` leme; variável de projeto `1`
  é o `it`.
- Arrasto: `CD = CDp + CDff` (Trefftz).
