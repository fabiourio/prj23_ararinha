# Handoff — Lab 04 (PRJ-23, equipe Ararinha): fechar a entrega e montar o relatório no Clause

Prazo: **hoje, 27/09/2026** (entrega no Google Classroom, zip `Ararinha_PRJ23_Lab04.zip`).
Repositório: `C:\Users\gilva\Desktop\ITA\4º PROF\2º SEMESTRE\PRJ-23\prj23_ararinha`, branch `main`.

## O que fazer nesta sessão, em ordem

1. **Conferir o estado.** Rode `git log --oneline -25` e `git status --short`.
   - Procure o commit *"Lab 04: CDq pela mesma definicao de arrasto (Trefftz) e notas do roteiro revistas"*.
   - Se ele não existir, aplique a correção (ver "Pendências" abaixo).
   - As exclusões não staged em `avl/` (arquivos antigos: `.jl`, `README.md` velho, `resultados/*` antigos) são **intencionais**: o Lab 04 foi refeito do zero. Nunca use `git add -A`, `git add .` ou `git commit -a`; adicione arquivos pelo nome. As exclusões podem entrar num commit próprio ("remove versao antiga do Lab 04") se o usuário concordar.
2. **Montar o relatório no Clause** (conector MCP `clauaw`, https://mcp.clause.studio/mcp).
   - Estrutura, conteúdo e números estão abaixo.
   - As figuras e tabelas estão em `avl/resultados/**` e em `relatorio_lab04/images|tables` (cópias).
   - Há um rascunho LaTeX parcial em `relatorio_lab04/` (`main.tex` e seções 01–03). Use-o como texto-base se servir.
3. **Escrever `avl/README.md`**: como rodar tudo, na ordem abaixo.
4. **Montar o zip** `Ararinha_PRJ23_Lab04.zip`, fora do git (`*.zip` está no `.gitignore`).
   - **Entra:** `avl/` (código, `fwd.avl`, `aft.avl`, `airfoils/`, `fuselage_*.dat`, `avl.exe`, `resultados/` sem `_tmp/` e sem `varredura/`, `tests/`, `README.md`), `designTool/` (sem `__pycache__`), `requirements.txt` e o PDF do relatório exportado do Clause.
   - **Não entra:** `737.avl`, `__pycache__`, `_tmp`.
5. **Push na `main`**, sem PR (fluxo da equipe). A conta ativa do `gh` costuma ser `devRWGilvan`, que **não** tem acesso ao repositório. Receita:
   ```
   gh auth switch --hostname github.com --user JoseJunior4144
   git -c "credential.helper=!gh auth git-credential" push origin main
   gh auth switch --hostname github.com --user devRWGilvan
   ```
   Os commits terminam com `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Pendências possivelmente abertas

**CDq na Tabela 9.** Deve usar a diferença finita de **Trefftz** (CDp + CDff): **≈ 0,1786 /rad**. É a mesma definição de CD0, CDα, CDα², CDit e CDδe.
- O valor de campo próximo, 0,503, fica só como sensibilidade. A razão campo próximo/Trefftz é ≈ 2,7× também para CDit e CDδe.
- Lugar da correção: `avl/derivadas_estabilidade.py`, perto da linha 205 (`fd['CDq_trefftz']`).
- Arredondar a Tabela 9 para 3–4 algarismos significativos.

**Notas do roteiro** (`avl/lab04_roteiro.py` → `resultados/roteiro/notas_roteiro.txt`):
- **Item 4:** o estol começa em η 0,493, dentro do aileron como exigido, mas só ~2 m para dentro, e toda a asa externa fica a 0,02–0,05 do clmax. Escrever "estol quase simultâneo na asa externa; critério atendido com margem pequena".
- **Item 7:** com o CG traseiro em M 0,85, dδe/dCL = +1,07° (gradiente invertido). A aeronave é estaticamente instável no cruzeiro pelo AVL.
- **Polares:** o VLM é linear, sem buffet nem compressibilidade em CL alto (CL até ~1,19 em M 0,85).

## Contexto do trabalho (o que o professor pediu e o que fizemos)

**Roteiro:** `Aulas/Documentos/lab04_prj23_2026.pdf`, com os itens 1–8 e a seção 3 (Tabelas 1, 4, 6, 7, 9).

**Pedido extra do professor**, depois de ver o AVL da equipe: a asa sombreava a EH em alto α. Ele pediu para recuar a EH até o limite da fuselagem, aproveitando o braço para reduzir o MTOW, e para otimizar a torção.

**Pipeline** (Python, em `avl/`; rodar de dentro de `avl/`; testes com `python -m pytest avl/tests -q`, cerca de 40 passando):

| Etapa | Script | O que faz |
|---|---|---|
| 0a | `convergencia_malha.py` | Convergência de malha |
| 0b | `verifica_eliptica.py` | Verificação do otimizador: asa limpa → carga elíptica |
| 1 | `varredura_eh.py` | Varredura da posição da EH (`Lc_h`) |
| 2 | `otimizacao_torcao.py` | Otimização de torção (~30 min; 3 casos em paralelo) |
| — | `ponto_neutro.py` | PN designTool × AVL, decomposto por componente |
| — | `lab04_roteiro.py` | Itens 1–8 |
| — | `derivadas_estabilidade.py` | Seção 3 |
| — | `estabilidade_direcional.py` | Decomposição de Cnβ |
| — | `sensibilidade_estabilidade.py` | Correção da margem: Cht ou posição da asa |

Módulos de apoio: `aeronave.py` (designTool → dicionário), `gera_avl.py`, `avl_run.py`, `avl_saida.py`, `analises.py`, `estol.py`, `sombra.py`, `estilo.py`.

**Cuidado no Windows:** o `python.exe` da Microsoft Store é só um launcher e termina logo; o processo real é o `python3.13.exe`. Para acompanhar processos longos, use o interpretador real: `/c/Users/gilva/AppData/Local/Microsoft/WindowsApps/PythonSoftwareFoundation.Python.3.13_qbz5n2kfra8p0/python.exe`.

### Decisões de modelagem (justificar no relatório)

- **Deslocamento em Z só no AVL** (asa −1,20 m, EH +1,85 m, EV +0,85 m, nacele junto com a asa). É exigência do professor, para manter os painéis afastados. O sombreamento é calculado com as alturas **reais** do designTool.
- **Winglet** como superfície própria no COMPONENT 1, com toe fixo em 0.
- **Nacele** como anel sustentador em **componente próprio**. No COMPONENT 1, o CDff oscilava ±0,2 count com a malha da asa. Foi decisão do usuário.
- **Fuselagem** com contorno **reamostrado**. O arquivo original, de 41 pontos com cauda rombuda, fazia a spline do AVL formar um laço, e o Nbody nunca convergia.
- **Estol pela seção crítica**, com o clmax dos perfis do Lab 03. As faixas da junção com o winglet e da raiz ficam excluídas: são singularidades de vórtice de quina, com área → 0 e cl ilimitado, e um avião real tem carenagem nesses pontos.
- **Malha adotada:** asa 8×60, winglet 4, EH 8×20, corpo 80. Tolerâncias: 0,1 count no CDff, 0,02° no `it`, 0,15° no α de estol. Incerteza residual do CDff ao dobrar a malha: ~±0,2 count.
- **Numeração dos controles**, diferente da do b737mod do professor: d1 aileron, d2 profundor, d3 leme; a variável de projeto 1 é o `it`.
- **Arrasto sempre CD = CDp + CDff (Trefftz).** O CDtot de campo próximo dá induzido irrealisticamente baixo com as naceles em anel.
- **Ponto de projeto:** M 0,85, h 10.668 m, peso de projeto com a fração de combustível do Lab 03 (0,451) e 100% de carga paga.

### Resultados principais

**Verificação do otimizador** (asa limpa plana):

| Partida | Oswald | Desvio RMS da elipse |
|---|---|---|
| Sem torção | 0,8715 → 0,9994 | 0,7% |
| Washout de −4° | 0,8715 → 0,9982 | 1,3% |

As duas chegam ao mesmo CDff, com diferença de 0,04 count.

**Posição da EH:**
- Varredura de `Lc_h` 4,0 → 4,8237 (fim da fuselagem − 0,5 m), com Cht fixo. W0 e CD caem de forma monótona.
- **`Lc_h` escolhido = 4,8237.**
  - W0: 291.292 (`Lc_h` 4,0) → 289.362 (4,6) → **288.764 kgf**.
  - S_h: 56,13 → **53,52 m²**.
  - O designTool já foi atualizado.
- **α_saída da sombra** cai com o recuo (12,0° → 9,3°).
- O **leme encoberto** pelas linhas de 60°/30° é informativo.
- A **margem estática não entra na escolha**: com Cht fixo, o PN praticamente não se move.

**Torção** (10 variáveis, SLSQP; objetivo CDff compensado no CL de projeto; restrição: estol começando para dentro de η 0,56 com desempate de 0,2°):

| η | 0 | 0,10 | 0,22 | 0,32 | 0,40 | 0,48 | 0,56 | 0,70 | 0,82 | 0,90 | 1,00 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| torção [°] | 0 | +3,00 | +1,13 | −0,26 | +0,94 | +3,00 | −3,32 | −3,28 | −5,02 | −4,65 | −3,75 |

- **CDff:** 104,0 → **94,0 counts** (CG traseiro). **Oswald:** 0,754 → 0,839.
- **Estol:** começa em η 0,493, com α 13,2° e CLmax 1,19 (CG traseiro) / 1,14 (dianteiro).
- **Sombreamento:** folga de **3,9°**, acima do mínimo de 2°.
- **Independência:** a torção muda no máximo 0,42° entre `Lc_h` 4,82, 4,62 e 4,42.
- Duas variáveis ficam no limite de +3°. O degrau entre η 0,48 e 0,56 é forçado pela restrição de estol.

**Item 2 (Tabela 1):**

| W0 | W | h | ρ | a | M | V | CL | Sref |
|---|---|---|---|---|---|---|---|---|
| 2.832.779 N | 2.249.127 N | 10.668 m | 0,380455 | 296,59 m/s | 0,85 | 252,10 m/s | 0,50439 | 368,833 m² |

**Item 3 (`it` que zera o δe no cruzeiro):**

| CG | it | α | CD | e |
|---|---|---|---|---|
| Dianteiro | −3,242° | 4,15° | 0,024707 | 0,787 |
| Traseiro | −0,964° | 3,94° | 0,024084 | 0,839 |

**Item 4 (estol, M 0,2):**

| Caso | α_max | CLmax | δe |
|---|---|---|---|
| Dianteiro, sem compensação | 13,18° | 1,173 | 0 |
| Dianteiro, compensado | 13,20° | 1,139 | −7,55° |
| Traseiro, sem compensação | 13,18° | 1,193 | 0 |
| Traseiro, compensado | 13,18° | 1,193 | +0,17° |

O início do estol fica em η 0,493 nos quatro casos.

**Itens 5–7:**
- CD no ponto de projeto, AVL contra designTool (0,023830): **+8,8 counts** (CG dianteiro) e **+2,5 counts** (traseiro).
- Profundor com o CG dianteiro: δe de −4,7° a +6,9° ao longo da faixa.
- Com o CG traseiro, o gradiente de profundor é invertido (instabilidade no cruzeiro).

**Item 8: ACHADO PRINCIPAL. Os pontos neutros divergem bastante.**

| Fonte | PN [% CMA] | MS com o CG traseiro (39,2% CMA) |
|---|---|---|
| AVL, M 0,85 | 36,5 | **−2,9%** |
| AVL, M 0,2 | 43,4 | +4,1% |
| designTool | 46,3 | +7,3% |

Decomposição por componente (variação do PN, em % CMA, AVL / designTool):

| Componente | AVL | designTool |
|---|---|---|
| Empenagens | +34/+38 | +31/+37 (concordam) |
| Asa isolada (c.a. em relação a c/4) | +9/+14 | 0 (fixo em c/4) |
| Fuselagem | −24/−34 | −13/−10 |
| Naceles | −7/−6 | não modela |

A fuselagem do AVL é teoria de corpo esbelto, com momento de Munk sem correção viscosa e amplificado pela correção de Prandtl-Glauert. Ela provavelmente superestima a instabilidade; a verdade deve estar entre as duas ferramentas.

**Fluxo do CG no designTool:**
- Cinco cenários de carregamento.
- O CG traseiro vem do cenário vazio + tripulação + carga paga, somado a 2% da CMA.
- O dianteiro vem do cenário com combustível sem carga paga, menos 2% da CMA.

**Sensibilidade**, a correção para a próxima iteração (`avl/resultados/sensibilidade/`):

| Estratégia | MS ≥ 5% (AVL, M 0,2) | MS ≥ 0 (AVL, M 0,85) |
|---|---|---|
| Aumentar o Cht | Cht 0,72, +339 kgf | Cht 0,77, +1.293 kgf |
| **Recuar asa + motores + trem principal** (empenagens fixas) | Δ 0,11 m, +82 kgf | **Δ 0,34 m, +267 kgf** |

- Com o recuo de 0,34 m, a margem pelo designTool sobe para +10%. O trem de pouso continua coerente (tipback acima de tailstrike) e a carga no trem de nariz melhora (6,3% → 7,3%).
- **Decisão:** o projeto entregue atende o critério da disciplina (designTool +7,3%). A divergência com o AVL está identificada, explicada e dimensionada, e a correção fica recomendada para a próxima iteração.

**Seção 3 (CG traseiro, ponto de projeto, sem compensação, `it` = −0,964°):**

| Tabela | Valores |
|---|---|
| 4 | Sref 368,833; cref 6,8995; bref 60,152; m 229.269 kg; Ixx 1,512e7; Iyy 4,297e7; Izz 5,595e7; Ixz 1,207e6 kg·m²; ip 0; xp 21,65 m (a partir do nariz; xp − xcg_aft = −6,07 m); zp +3,0 m; Tmax 195.552 N (2 motores na altitude); V 252,10 m/s; h 10.668 m |
| 6 | CL0 0,0367; CM0 −0,0660 |
| 7 (ajuste de polar 5b, α em rad) | CD0 0,016734; CDα −0,003425; CDα² 1,6007 |

Tabela 9:
- **Sustentação:** CLα 6,920; CLq 6,866; CLit 0,6430; CLδe 0,3495.
- **Arrasto:** CDq **0,179 (Trefftz)**; CDit 0,01112; CDδe 0,00670.
- **Arfagem:** CMα **+0,203**; CMq −36,42; CMit −2,984; CMδe −1,675.
- **Força lateral:** CYβ 0,760; CYp 0,0698; CYr 0,0228; CYδr −0,2325.
- **Rolamento:** Clβ −0,2158; Clp −0,6029; Clr 0,1326; Clδa 0,2409; Clδr 0,0305.
- **Guinada:** Cnβ **−0,0811**; Cnp −0,0511; Cnr −0,2073; Cnδa 0,0160; Cnδr −0,1259.

As conversões seguem o roteiro: controles ×180/π; inverter o sinal de CYβ, CYp, CYr, Clδa, Clδr, Cnδa e Cnδr; derivadas látero-direcionais de momento pelo `sb`.

Observações para o texto da seção 3:
- As inércias são em torno do CG do ponto de projeto (x = 27,005 m), não do CG traseiro. Isso dá +0,27% em Iyy/Izz. Ixx/Iyy são tomados em torno de z = 0.
- A origem de coordenadas usada pelo professor de MVO **ainda precisa ser confirmada**.
- O CD0 da Tabela 7 já inclui ~CDit·it ≈ −1,9 counts, porque o ajuste é feito no `it` do item 3.
- **CMα > 0 e Cnβ < 0:** as derivadas descrevem uma condição instável pelo AVL.
- **Decomposição de Cnβ** (`avl/resultados/derivadas/decomposicao_cnb.csv`), M 0,2 / 0,85:

  | Modelo | Cnβ |
  |---|---|
  | Completo | +0,016 / −0,081 |
  | Sem fuselagem | +0,133 / +0,152 |
  | Só a fuselagem | −0,130 / −0,243 |

  A instabilidade vem do momento de Munk da fuselagem e dos anéis das naceles; só com a EV a aeronave seria estável. Recomenda-se conferir com métodos empíricos (DATCOM/Roskam) e rever o Cvt.

### Dados pedidos por um membro da equipe (para MVO)

| Grandeza | Valor |
|---|---|
| Posição da espessura máxima dos perfis da asa | raiz 40,5%, centro 35,8%, ponta 32,8% da corda (t/c 21,7%, 17,6%, 10,8%) |
| Incidência da asa | 0° na raiz (define o eixo do corpo) |
| α de sustentação nula da asa (AVL, asa isolada com torção) | −0,62° em M 0,85; −0,56° em M 0,2 |
| Cm dos perfis da asa (c/4, cl de projeto, M 0,85, Euler Lab 03) | raiz −0,057; centro −0,097; ponta −0,109 (subsônico: −0,050) |
| Torção da asa | tabela acima; ponta −3,75° em relação à raiz |
| Área da EH | 53,52 m² |
| Espessura máxima do perfil da EH | 30% da corda (NACA 0010) |
| Fração da corda ocupada pelo profundor | 30% (charneira em 70%) |

## Estrutura sugerida do relatório (português, para o professor Ney Sêcco)

1. Introdução: objetivos, pedido extra do professor, conteúdo do zip.
2. Modelo no AVL e ponto de projeto (itens 1–2).
3. Verificações numéricas: malha e otimizador.
4. Posição da EH: varredura e sombreamento.
5. Otimização de torção.
6. Incidência da EH e plano de Trefftz (item 3).
7. Estol pelo método da seção crítica (item 4).
8. Polares, CL×α e CL×δe (itens 5–7).
9. Ponto neutro e margem estática (item 8) + comparação designTool × AVL + sensibilidade e recomendação.
10. Derivadas de estabilidade (Tabelas 4, 6, 7, 9) e decomposição de Cnβ.
11. Conclusões e recomendações.

Toda figura deve ser referenciada e discutida, e toda escolha justificada em 1–2 frases.

### Figuras

| Pasta | Arquivos |
|---|---|
| `avl/resultados/` | `convergencia_malha.png`, `verifica_eliptica.png`, `varredura_eh.png`, `ponto_neutro.png`, `otimizacao_torcao.png`, `carga_torcao.png` |
| `avl/resultados/roteiro/` | `trefftz_fwd.png`, `trefftz_aft.png`, `trefftz_*_avl.png` (hardcopy do próprio AVL), `estol_cl_y.png`, `polar_cd_cl.png`, `cl_alpha.png`, `cl_delta_e.png` |
| `avl/resultados/derivadas/` | `ajuste_cd_alpha.png`, `decomposicao_cnb.png` |
| `avl/resultados/sensibilidade/` | `sensibilidade_ms.png`, `sensibilidade_w0.png` |

Duas legendas precisam de observação:
- em `cl_delta_e.png`, os quatro pontos de projeto coincidem em (0°; 0,504);
- em `estol_cl_y.png`, as distribuições com e sem compensação coincidem.
