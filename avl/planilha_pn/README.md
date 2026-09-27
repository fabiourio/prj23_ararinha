# Planilha de ponto neutro da disciplina (Torenbeek)

Planilha do professor (`estudo_cg_pn_deflexoes_aula`), preenchida pelo Rafael
com os dados da configuração final (`Lc_h` = 4,8237).

| arquivo | conteúdo |
|---|---|
| `estudo_cg_pn_original_rafael.xls` | como o Rafael preencheu |
| `estudo_cg_pn_corrigida.xls` | a mesma planilha com as correções abaixo |

## Correções na aba "Estimativa de PN"

1. **`xr_h` (C44): de 30,37 m para 56,185 m.** É o bordo de ataque da raiz da EH
   no designTool. Com 30,37 m, a alavanca da EH (`L_h`) saía 7,5 m em vez de
   33,3 m, e o dε/dα subia para 0,60. Isso derrubava o ponto neutro para
   19% da CMA.
2. **CG dianteiro e traseiro (C91, C92): 15,45% e 39,21% da CMA,** os limites
   do designTool.

## Pendente

As células da nacele (C75 a C79: `ln_front`, `ln_back`, `ln`, `k_n`,
`dxac_n`) estão com 0 fixo, sem fórmula, então a correção de Torenbeek
(Eq. E-41) para as naceles não entra no cálculo. As fórmulas precisam ser
restauradas a partir da planilha original do professor.

## Resultado (planilha corrigida, M 0,85)

| | sem a nacele (como está) | com a nacele (estimativa) |
|---|---|---|
| x_np | 28,10 m = 44,7% da CMA | ≈ 27,69 m ≈ 38,7% da CMA |
| MS com CG traseiro | +5,5% | ≈ −0,5% |

A estimativa com a nacele usa k_n = −4 (nacele à frente da asa), b_n = 4 m,
l_n = x_m,w − x_n = 7,0 m e 2 naceles, o que dá Δx_ac = −0,44 m. Com a
nacele, o método semi-empírico fica perto do AVL (36,3% da CMA no
cruzeiro), e o designTool (46,5%), que ignora as naceles, fica como o valor
otimista.
