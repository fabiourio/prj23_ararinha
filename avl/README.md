# Lab 04 -- analises no AVL

## Estrutura

| Item | O que e |
|---|---|
| `avl.exe` | Executavel do AVL (3.37), fornecido na disciplina |
| `fwd.avl` / `aft.avl` | Entradas do AVL, CG dianteiro (26,0772 m) e traseiro (27,7190 m). Gerados por `lab04_gera_avl.py` (raiz do repo) -- nao editar na mao |
| `fwd.mass` / `aft.mass` | Massa, CG e inercias do ponto de projeto (so para o comando `mode`) |
| `aerofolios/` | Perfis otimizados do Lab 03 (familia roteiro), reamostrados para 197 pontos (limite IBX deste build) |
| `dados_lab04.md` | Ponto de projeto, CGs, polar e inercias extraidos do designTool (`lab04_dados.py`) |
| `saidas/` | Despejos das rodadas: `fs.txt`, prints do plano de Trefftz, etc. |

## Como rodar

Sempre de dentro desta pasta (os `AFILE` sao relativos a ela):

```
cd avl
avl.exe
load fwd.avl
```

## Diferencas em relacao aos comandos do roteiro (b737mod.avl)

A numeracao dos nossos arquivos e outra:

| Roteiro (b737) | Aqui | O que e |
|---|---|---|
| `d4 pm 0.0` | `d2 pm 0.0` | trimagem pelo profundor |
| menu `de`: `2 <valor>` | menu `de`: `1 <valor>` | incidencia da empenagem (it) |

Controles: `d1` aileron, `d2` elevator, `d3` rudder.

## Metas do ponto de projeto (Lab 03, peso medio de cruzeiro)

- `mn 0.85` (cruzeiro) ou `mn 0.2` (baixa velocidade, metodo da secao critica)
- `a c 0.5053` (CL de projeto)
- CDp ja embutido no cabecalho: 0,01489 (CD0 do designTool)
