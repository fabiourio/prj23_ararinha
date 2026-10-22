# Atividade 06 — Cargas simétricas na asa (equipe Ararinha)

Diagramas de esforços da asa em manobra simétrica, nas condições de início
e fim do cruzeiro. Usa o `designTool/` (pesos, geometria) e o modelo do AVL
do Lab 04 (`avl/fwd.avl`, `avl/avl_run.py`).

## Requisitos

Os mesmos do Lab 04 (`avl/README.md`): Windows, `avl.exe` 3.37 em `avl/`,
Python com `numpy`, `scipy`, `matplotlib` e `pytest`.

Rode os scripts de dentro de `cargas/`.

## Como rodar

| # | Script | Resultado |
|---|---|---|
| 0 | `python clmin_perfis.py` | cl mínimo dos perfis do Lab 03 no XFoil (alguns minutos; o CSV já está em `resultados/01_envelope/`) |
| 1 | `python envelope.py` | CLmin, CLα e diagramas V-n (`resultados/01_envelope/`) |
| 2 | `python simulacao_va.py` | AVL em V_A e n limite, trimado (`resultados/02_avl/`) |
| 3 | `python esforcos.py` | cargas distribuídas, cortante, fletor e torçor (`resultados/03_esforcos/`) |

Testes: `python -m pytest tests -q`.

## Hipóteses

- Caso crítico: CG dianteiro (`fwd.avl`) e V_A ao nível do mar. Os V-n ficam
  na altitude de cruzeiro.
- V_C pelo Mach de cruzeiro (0,85) e V_D por M_D = M_C + 0,07 (25.335(b)(2)).
  Rajada com U_ref de 25.341(a)(5) e fator de alívio de Pratt.
- CLmax (1,139) e CLmin pelo método da seção crítica em M 0,2, compensados.
- Estrutura da asa e combustível proporcionais à corda; combustível só no
  tanque (0 a 90% da semienvergadura). Motor instalado (fator 1,3 do
  designTool) e trem principal concentrados.
- Eixo elástico a 32,5% da corda (centro da caixa entre 12% e 53%).
- A sustentação da nacele (anel do AVL) entra como carga concentrada no
  pilone.
