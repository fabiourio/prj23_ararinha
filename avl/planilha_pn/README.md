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

3. **Células da nacele reconstruídas (C75 a C80).** Na cópia recebida elas
   estavam com 0 fixo, sem fórmula. Foram reconstruídas pela Eq. E-41 de
   Torenbeek, seguindo os rótulos da planilha:
   - `ln_front` = `xm_w − x_n`: bordo de ataque da CMA menos a face frontal
     da nacele;
   - `ln_back` = `x_n + L_n − xm_w`;
   - `ln` = `ln_front` para motor na asa e `ln_back` para motor na cauda;
   - `k_n` = −4 para nacele à frente da asa e −2,5 para nacele na cauda;
   - `dxac_n` = `num_eng · k_n · D_n² · ln / (S_w · CLα_wf)`;
   - `xac_wfn` = `xac_wf + dxac_n`.
4. **`x_np` (C95) passa a usar `xac_wfn` (C80)** em vez de `xac_wf` (C74).
   Antes, a nacele ficava de fora mesmo com as células preenchidas.

A planilha foi salva com Mach 0,85.

## Resultado (planilha corrigida)

| | sem a nacele | com a nacele |
|---|---|---|
| x_np, M 0,85 | 28,10 m = 44,7% da CMA | 27,69 m = 38,7% da CMA |
| MS com CG traseiro, M 0,85 | +5,5% | −0,5% |
| x_np, M 0,2 | 28,31 m = 47,8% da CMA | 27,77 m = 40,0% da CMA |
| MS com CG traseiro, M 0,2 | +8,6% | +0,8% |

A nacele desloca o centro aerodinâmico do conjunto asa-fuselagem em
−0,44 m. Com ela, o método semi-empírico fica perto do AVL (36,3% da CMA no
cruzeiro e 43,3% em M 0,2). O designTool (46,5%), que ignora as naceles,
fica como o valor otimista.
