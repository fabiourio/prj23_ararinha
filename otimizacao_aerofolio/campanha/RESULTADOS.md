# Lab 03 — Resultados consolidados

Equipe Ararinha | Projeto de perfil transônico para a aeronave B

Este documento reúne **só os resultados**. O raciocínio por trás de cada
escolha está em `DEFINICAO_DO_PROBLEMA.md`; como rodar está no `README.md`.

---

## 1. Ponto de projeto

Aeronave **B** (joelho da frente de Pareto do Lab 02), no cruzeiro, com peso
médio de missão.

| | |
|---|---|
| altitude | 35.000 ft (10.668 m) |
| M∞ | 0,850 |
| peso (médio de cruzeiro) | 229.669 kgf |
| CL da aeronave | 0,5053 |
| enflechamento a c/2 | 32,158° |
| **M_n = M·cos Λ** | **0,7196** |

O perfil **não** vê M = 0,85. A transformação de enflechamento é a mesma que
a equação de Korn do `designTool` já usa internamente, então adotá-la mantém
o Lab 03 coerente com o modelo que produziu a aeronave.

Distribuição de sustentação adotada: **elíptica** — é o alvo da futura
otimização de torção, então projetar para ela é projetar para a asa que
pretendemos ter.

### As três estações

| estação | η | corda [m] | (t/c)ₙ | c_ℓ,ₙ | Re,ₙ |
|---|---|---|---|---|---|
| raiz (junção asa-fuselagem) | 0,101 | 9,130 | 0,2179 | 0,600 | 4,34×10⁷ |
| meio (MAC — a do roteiro) | 0,398 | 6,899 | 0,1772 | 0,732 | 3,28×10⁷ |
| ponta (fim do tanque) | 0,900 | 3,125 | 0,1082 | 0,768 | 1,49×10⁷ |

---

## 2. Os perfis otimizados — resultado principal

Formulação: `min c_d` sujeito a `c_ℓ = c_ℓ,ref`, `(t/c) ≥ (t/c)ref`, bordo de
fuga fino e a restrição de sustentação máxima. SLSQP com gradiente adjunto.

| | raiz | meio | ponta |
|---|---|---|---|
| c_d de partida | 0,071561 | 0,055979 | 0,024822 |
| **c_d otimizado** | **0,019857** | **0,011644** | **0,008872** |
| ganho | −72,3 % | −79,2 % | −64,3 % |
| ganho reavaliado em malha fina | −71,0 % | −85,9 % | −81,7 % |
| c_ℓ atingido | 0,5994 | 0,7319 | 0,7677 |
| α do ótimo | 3,457° | 2,389° | 1,573° |
| iterações / avaliações | 11 / 14 | 17 / 31 | 21 / 120 |

### Geometria dos ótimos

| | x da espessura máxima | arqueamento máx. | x do arq. máx. |
|---|---|---|---|
| raiz | 0,410 | −0,0300 | 0,221 |
| meio | 0,362 | −0,0152 | 0,249 |
| ponta | 0,334 | +0,0240 | 0,561 |
| *(partida: NACA 4 dígitos)* | *0,300* | *+0,0101* | *0,418* |

O otimizador **redescobriu o perfil supercrítico** partindo de um NACA de 4
dígitos: espessura máxima recuada, arqueamento dianteiro negativo e carga
jogada para trás. Na distribuição de Mach, trocou um choque normal forte
(M = 1,51 despencando para 0,63) por um **platô supersônico** em M ≈ 1,15-1,20
com recompressão gradual. Quanto mais grossa a seção, mais atrás ele precisou
jogar a espessura — a assinatura fica mais pronunciada onde o problema é pior.

---

## 3. A restrição de sustentação máxima

O código de Euler é não-viscoso e não enxerga cl_max. Sem amarra, nada impede
o otimizador de afiar o bordo de ataque para reduzir arrasto de onda e
destruir a sustentação máxima.

**Ela é ativa nas três estações.** Comparando as duas rodadas na MAC:

| | c_d | cl_max no XFoil |
|---|---|---|
| com a restrição | 0,011644 | **2,036** |
| sem a restrição | 0,010664 | **1,042** |

**Custo: +9,2 % de arrasto. Compra: 0,99 de cl_max.**

Sem ela, o `CLmaxTO` da aeronave desabaria e os motores ficariam no limite
para decolar nos 2.900 m de pista — o orçamento de sustentação máxima põe o
esgotamento de empuxo em `clmax_w = 1,20`.

### Como ela entrou na formulação

O cl_max do XFoil **não serve** como restrição de gradiente: medindo a
derivada por diferença central em cinco passos, ela espalha **1.617 %** e
troca de sinal (o adjunto do Euler, no mesmo teste, varia 0,04 %). O ruído é
intrínseco, da quantização do `max` sobre a grade de α.

A saída foi uma **substituta geométrica**: `t₀₁/√(t/c) ≥ limiar`, onde `t₀₁` é
a espessura a 1 % da corda. Foi o descritor de maior correlação com cl_max num
DOE de 277 perfis — melhor que o raio de bordo de ataque e que o Δy de Abbott.
É linear nos coeficientes CST, então a derivada é exata e de graça. Limiar
calibrado por estação: 0,1220 na raiz, 0,0938 na MAC, 0,0978 na ponta.

---

## 4. Verificação

Nada foi aceito sem conferência independente.

| o quê | resultado |
|---|---|
| **Convergência de malha** | o nível padrão superestima c_d em **70 %** (Richardson, ordem observada 2,43, c_d extrapolado 0,007308) |
| **O ganho sobrevive ao refinamento?** | sobrevive e **cresce** (−79 % → −86 % na MAC) |
| **Gradientes adjuntos** | erro de 0,00 a 1,8 % contra diferenças finitas; pior componente `dc_ℓ/dA_u0` com 4,3 % |
| **Viés do XFoil em cl_max** | **+0,24** acima do experimental (5 perfis NACA, Abbott & von Doenhoff); α de estol 19-20,5° medidos contra 14-17° reais |

Aplicando o viés, os cl_max **reais** ficam em 1,72 na raiz, **1,80 na MAC**
(no alvo) e 1,90 na ponta.

Por que o ganho cresce com o refinamento: a dissipação numérica adiciona um
piso de arrasto quase **absoluto**. No perfil de partida, dominado por choque
forte, ele pesa pouco em proporção; no otimizado, de arrasto baixo, ele
domina. Otimizar na malha padrão portanto **subestima** o benefício.

---

## 5. O que isso diz sobre a aeronave

### 5.1 O `k_korn = 0,95` está correto

| estação | Korn previa | Euler mediu | k_korn implícito |
|---|---|---|---|
| raiz | 0,011601 | 0,015218 | 0,939 |
| meio | 0,005321 | 0,007145 | 0,940 |
| ponta | 0,000301 | 0,003858 | 0,894 |

A física de seção do modelo está certa. Ele apostava que a equipe entregaria
bons perfis, e o Lab 03 mostra que é alcançável.

### 5.2 Mas a agregação subestima o arrasto de onda em ~8×

| | CDwave |
|---|---|
| integração das seções ao longo da envergadura | **0,00403** |
| fórmula do `designTool` | 0,00051 |

A causa é aritmética: `(t/c)ₘ = 0,25·0,196 + 0,75·0,080 = 0,109`, que na
variação linear de espessura é a espessura da estação **η = 0,750** — onde a
corda já caiu para 43 % da raiz. A fórmula representa a asa inteira por uma
estação a três quartos da semi-envergadura, pulando a região interna, que é
onde a asa é grossa **e** tem corda longa.

**Ressalva:** nossa integração usa seções 2D na raiz, e a raiz real é aliviada
pela fuselagem e pela carenagem. O valor verdadeiro está entre 0,00051 e
0,00403, provavelmente mais perto do meio. A direção é sólida; a magnitude
tem incerteza.

### 5.3 A asa é grossa demais para M = 0,85

Margem até a divergência de arrasto (M_dd − M_n), variando só a espessura de
raiz:

| t/c raiz | raiz | MAC | ponta | CDwave |
|---|---|---|---|---|
| **0,196 (a nossa)** | **−0,048** | **−0,020** | +0,045 | 0,00244 |
| 0,160 | −0,009 | +0,006 | +0,050 | 0,00093 |
| 0,150 | +0,002 | +0,013 | +0,051 | 0,00068 |
| 0,140 | +0,012 | +0,020 | +0,052 | 0,00049 |

Raiz e MAC operam **além da divergência de arrasto**. Com margem negativa,
nenhum perfil elimina arrasto de onda — ele é imposto pela espessura e pelo
Mach, e a otimização só escolhe *como* gastá-lo.

E a distribuição de espessura é atípica: **19,6 % na raiz** contra 14-16 % da
indústria, e **8,0 % na ponta** contra 9-11 %. Esse padrão invertido é o que a
ponderação `0,25·tcr + 0,75·tct` recompensa — engordar a raiz custa um quarto
do peso na conta e afinar a ponta rende três quartos do alívio.

### 5.4 Recomendação

A equipe **decidiu não reotimizar o Lab 02** — a cascata (MTOW, frente de
Pareto, escolha da aeronave) não cabe no escopo, e a espessura de raiz serve a
um propósito estrutural real. A recomendação é **documentar e corrigir o
modelo numa iteração futura**: trocar a média ponderada por uma integração ao
longo da envergadura é mais barato e mais correto do que afinar `tcr_w`, e não
mexe na aeronave.

---

## 6. Item 9 — alterar a definição do problema

### 6.1 Os ótimos estavam no batente, não no interior

Com os batentes sugeridos pelo roteiro havia **cinco batentes ativos**:

| estação | batentes ativos |
|---|---|
| raiz | `Al4 = −0,05` e `Au2 = +0,05` |
| meio | `Al4 = −0,05` |
| ponta | `Al3 = −0,05` e `Al4 = −0,05` |

O `Al4` é o coeficiente que constrói o **cusp côncavo do bordo de fuga**, a
essência do carregamento traseiro supercrítico. O RAE2822 tem esse coeficiente
em **+0,052**; o batente do roteiro (`Al ≤ −0,05`) proíbe essa forma.

Soltando apenas ele (`Al4 ≤ +0,30`):

| estação | c_d com batente | c_d com cusp | ganho |
|---|---|---|---|
| raiz | 0,019857 | 0,009197 | **−53,7 %** |
| meio | 0,011644 | 0,008215 | **−29,4 %** |
| ponta | 0,008872 | 0,008429 | −5,0 % |

Os ótimos passam a ser **interiores**. E o `Al4` da ponta convergiu para
**+0,050**, praticamente o valor do RAE2822 — o otimizador chegou sozinho ao
mesmo número de um supercrítico real.

### 6.2 Mas o ganho não sobrevive à aeronave

O momento de arfagem quase dobra, e a nossa função objetivo é 2D: ela não
enxerga a empenagem.

| | arrasto de onda da asa | c_m médio | CD de trimagem |
|---|---|---|---|
| batentes do roteiro | 0,004013 | −0,086 | 0,000965 |
| cusp liberado | 0,002598 | −0,184 | **0,002644** |

```
ganho de arrasto de onda    +0,001415
custo de compensação        −0,001680
─────────────────────────────────────
saldo líquido               −0,000264     ← PIOR
```

**O perfil isolado melhora 29 % e a aeronave piora.** Momento maior exige mais
download da empenagem, que a asa compensa com mais sustentação, que volta como
arrasto induzido — na asa e na própria empenagem.

É o caso de livro de **otimização de subsistema sem acoplamento**: o
otimizador gastou livremente uma moeda que não paga. Para o cusp valer, a
formulação precisa de restrição de c_m.

*(Cálculo de ordem de grandeza: equilíbrio rígido, sem contribuição da
fuselagem, arrasto induzido parabólico.)*

### 6.3 A polar tem poço de arrasto estreito

No ponto de projeto o c_d é 0,0118; em c_ℓ = 0,54 ele sobe para 0,0185 —
**57 % maior abaixo do ponto de projeto**. É a assinatura da otimização
mono-ponto, e importa porque o CL de cruzeiro da aeronave varia de **0,612 a
0,417** ao longo da missão. Justifica uma rodada multiponto.

---

## 7. Lições de método

Cada uma custou tempo e vale para o próximo laboratório.

**Verifique se o ótimo é interior antes de discutir o custo das restrições.**
Gastamos tempo comparando formulações da restrição de cl_max (que vale 9 %)
enquanto a limitação dominante eram os batentes (28 a 54 %), invisível até
alguém verificar.

**`ftol` precisa ficar acima do piso de ruído.** Com `ftol = 1e-6`, a raiz
convergiu na avaliação 16 e gastou mais 123 sem parar — o c_d oscila 1,5×10⁻⁵
depois de convergido, quinze vezes o `ftol`.

**Restrição redundante não é inofensiva.** Somar `mint ≥ 0` ao `mint ≤ 0,01`
põe duas linhas antiparalelas na jacobiana e a deixa deficiente em posto *por
construção* — o SLSQP morreu com "Singular matrix E in LSQ subproblem", mesmo
com as restrições inativas.

**Ótimos planos existem.** A ponta convergiu na avaliação 21 e gastou mais 99
passeando: o t/c varia de 0,1094 a 0,1111 com o c_d constante na sétima casa.

**Não anote número à mão.** O primeiro balanço de trimagem usou `c_m` copiados
de rodadas parciais e concluía que o cusp valia a pena (+13 % de saldo). Com
os valores lidos do Euler nos ótimos convergidos, o saldo é **−19 %**.

---

## 8. Defeitos encontrados no material fornecido

1. **`export_airfoil` mata o XFoil em silêncio.** Com o padrão
   `close_te=True` repete o ponto do bordo de fuga que o `cstfoil` já fechou;
   o terceiro ponto coincidente vira painel de comprimento zero e o XFoil
   encerra com código de retorno 0 e sem mensagem. Afeta o `single_run.py`.

2. **`cstfoil` devolve posição errada da espessura máxima** — 0,117 em vez de
   0,30 para o NACA 1411. A linha 216 filtra o vetor de espessuras e a 219
   indexa o de abscissas completo. **Afeta a coluna `x_t/c,max` da Tab. 2.**

3. **`cstfoil(Au, Al, …)` e `run_cst(Al, Au, …)` têm ordem invertida.** Trocar
   gera perfil de dentro para fora; o erro que chega ao Python é um
   `FileNotFoundError` que não diz nada sobre a causa.

---

## 9. Figuras

Por estação (`raiz`, `meio`, `ponta`):

| arquivo | item do roteiro |
|---|---|
| `otim_<est>_convergencia.png` | 4 |
| `otim_<est>_geometria.png` | 5 |
| `otim_<est>_cp_mach.png` | 5 e 6 |
| `polar_transonica_<est>.png` | 7 |
| `subsonico_<est>.png` | 8 |
| `item9_<est>.png` | 9 |

Do DOE de sustentação máxima: `doe_clmax_corte.png`,
`doe_clmax_dispersao.png`, `doe_clmax_limiar.png`.
