# Dados do Lab 04 -- gerados por lab04_dados.py

Ponto de projeto: o mesmo do Lab 03 (peso medio de cruzeiro).

## Tabela 1 -- ponto de projeto (item 2)

| Parametro | Valor |
|---|---|
| W0 [N] | 2857571.4 (291291.7 kgf) |
| W [N] | 2253055.8 (229669.3 kgf) |
| h [m] | 10668.0 |
| rho_inf [kg/m3] | 0.38046 |
| a_inf [m/s] | 296.587 |
| M | 0.85 |
| V [m/s] | 252.1 |
| CL | 0.5053 |
| Sref [m2] | 368.833 |

## Cabecalho dos arquivos .avl (item 1)

| Parametro | Valor | Uso no .avl |
|---|---|---|
| Sref [m2] | 368.833 | linha Sref Cref Bref |
| Cref = MAC [m] | 6.8995 | linha Sref Cref Bref |
| Bref [m] | 60.1518 | linha Sref Cref Bref |
| Xref (CG dianteiro) [m] | 26.0772 | fwd.avl (15.4 %MAC) |
| Xref (CG traseiro) [m] | 27.7190 | aft.avl (39.2 %MAC) |
| CDp = CD0 designTool | 0.01489 | linha CDp |

## CG, ponto neutro e margem estatica (designTool, item 8)

| Parametro | [m] | [%MAC] |
|---|---|---|
| xcg_fwd | 26.0772 | 15.4 |
| xcg_aft | 27.7190 | 39.2 |
| xcg no ponto de projeto | 27.0124 | 29.0 |
| xnp | 28.1467 | 45.4 |

SM_fwd = 30.0 %MAC, SM_aft = 6.2 %MAC
(valores do designTool -- o item 8 pede os equivalentes via AVL)

## Polar do designTool no ponto de projeto (itens 1 e 5)

| Parametro | Valor |
|---|---|
| CD0 (parasita, config. limpa) | 0.01489 |
| CDind | 0.00867 |
| CDwave | 0.00051 |
| CD total | 0.02407 |
| K (fator de arrasto induzido) | 0.03398 |
| e (Oswald) | 0.7958 |
| CLmax limpa (estimativa designTool) | 1.333 |

## Momentos de inercia -- Tabelas 4 e 9 do MVO

Carregamento: fuel_frac = 0.4392 (reproduz W = 229669.3 kgf), payload_frac = 1.0. Referencia: CG deste carregamento.

| Parametro | Valor |
|---|---|
| m [kg] | 229669.3 |
| Ixx [kg m2] | 1.5083e+07 |
| Iyy [kg m2] | 4.2843e+07 |
| Izz [kg m2] | 5.5780e+07 |
| Ixz [kg m2] | 1.2091e+06 |
| Ixy [kg m2] | -1.1642e-10 |
| Iyz [kg m2] | 0.0000e+00 |

## Demais itens da Tabela 4 (motor)

| Parametro | Valor | Observacao |
|---|---|---|
| ip [deg] | 0.0 | designTool nao define incidencia de motor |
| xp [m] | 18.000 | face frontal da nacele (x_n) |
| zp [m] | -3.000 | inverter o sinal para o MVO |
| Tmax por motor [N] | 513000 | n_engines = 2 |
| T0 (thrust matching) [N] | 1026000 | tracao total de decolagem |

