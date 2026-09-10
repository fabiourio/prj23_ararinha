# Projeto de perfil para a ararinha_23 — definição do problema

Lab 03 do PRJ-23 (Maj. Ney Rafael Sêcco), entrega em 13/09/2026.
Equipe Ararinha. Documento escrito em 10/09/2026.

O roteiro pede a otimização de **uma** seção. Este documento define uma
campanha maior, que usa a seção do roteiro como primeiro caso e estende para
três estações da asa, com o objetivo declarado de **realimentar o Lab 02**:
decidir a distribuição de espessura da asa com base em aerodinâmica de maior
fidelidade, em vez das duas constantes que a equipe assumiu.

---

## 1. Por que este laboratório existe, do ponto de vista do projeto

O `designTool` usa o perfil aerodinâmico em exatamente dois lugares, e em
ambos por meio de uma **constante assumida**:

| Constante | Valor no Lab 02 | Onde entra | Consequência |
|---|---|---|---|
| `k_korn` | 0,95 | `aerodynamics.py:269`, equação de Korn | arrasto de onda de cruzeiro |
| `clmax_w` | 1,8 | `aerodynamics.py:282` | `CLmax_clean`, daí `CLmaxTO`, daí distâncias de decolagem e pouso, daí `S_w` |

A cadeia da segunda constante é direta e quantificada:

```
CLmax_clean = 0,9 · clmax_w · cos Λ = 0,9 · 1,8 · 0,82263 = 1,3327
CLmaxTO     = CLmax_clean + ΔCL_flape + ΔCL_slat = 2,2196   (designTool)
```

Na configuração de decolagem, verificado: `1,3327 + 0,7174 + 0,1696 = 2,2196`.
(Cuidado ao ler esses números: o dicionário que fica em
`airplane['aerodynamics']` depois do `analyze` é o da **última** chamada de
`aerodynamics`, que é outra configuração — é preciso chamar explicitamente
com `highlift_config='takeoff'`.)

O que importa é a **sensibilidade**: as contribuições de flape e slat não
dependem de `clmax_w`, então cada unidade de `cl_max` do perfil vale
`0,9 · cos Λ = 0,7404` de `CLmaxTO`, um-para-um — conferido numericamente
varrendo `clmax_w` de 1,4 a 2,0.

### 1.1 O orçamento de cl_max

Vale a pena saber **quanto** de sustentação máxima o perfil pode perder antes
de a aeronave precisar mudar, porque isso transforma a restrição de um
"1,8 ou falhou" numa escala com significado físico. Varrendo `clmax_w` na
aeronave B:

| `clmax_w` | `CLmaxTO` | `T0req/T0` | `deltaS_wlan` |
|---|---|---|---|
| 1,80 (assumido) | 2,2196 | 0,800 | 106,6 |
| 1,50 | 1,9975 | 0,889 | 81,2 |
| 1,35 | 1,8865 | 0,942 | 66,6 |
| **1,20** | **1,7754** | **1,000** | **50,4** |
| 1,00 | 1,6273 | 1,091 | 25,9 |

Dois resultados:

- **A restrição ativa é a decolagem, por empuxo — não o pouso.** A folga de
  área de asa para pouso (`deltaS_wlan`) permanece positiva mesmo em
  `clmax_w = 1,0`. Quem limita é o empuxo requerido para os 2.900 m de pista.
- Com a geometria da aeronave B **congelada**, o empuxo esgota em
  `clmax_w ≈ 1,20`.

### 1.2 Por que o alvo é 1,80 e não 1,20

O 1,20 da tabela acima é real, mas **não serve como alvo de projeto**, por
três razões que se acumulam:

1. **O ótimo do Lab 02 foi obtido com `clmax_w = 1,8`.** A varredura acima
   congela a geometria da aeronave B e mexe numa constante — responde "esta
   aeronave aguenta?", não "que aeronave o problema produz?". Com
   `clmax_w = 1,20` desde o início, o NSGA-II teria convergido para outra
   aeronave, provavelmente de asa maior. Só a segunda pergunta importa para
   projeto.
2. **`T0` é resultado do casamento de empuxo, não dado.** Quando `T0req`
   excede `Tmax`, a resposta de projeto não é "inviável", é motor maior ou
   asa maior. O 1,20 marca o fim *desta escolha de motor*, e o
   `Tmax = 513 kN` é herança do PRJ-22.
3. **1,80 é o número com o qual a aeronave foi otimizada.** Segurar o perfil
   em 1,80 mantém o Lab 02 auto-consistente. Descer abaixo disso gasta,
   calado, margem que pertence a decisões anteriores.

**Decisão: o alvo é `cl_max = 1,80`.** O 1,20 entra no relatório apenas como
contexto — mostra que não estamos à beira de um precipício, e quantifica o
que se perde caso o alvo se revele caro demais em arrasto de onda.

**O Lab 03 mede essas duas constantes.** É esse o produto.

### A hipótese que motiva a campanha

Na equação de Korn a espessura entra como média ponderada
`(t/c)_m = 0,25·tcr + 0,75·tct`, com `tct` fixo em 0,08 (não era variável de
projeto no Lab 02). O otimizador do Lab 02 engordou a raiz de 0,18 para
0,196 porque isso alivia a estrutura, e o modelo cobrou pouco por isso:
no ponto de cruzeiro o arrasto de onda previsto é de apenas **2,1 % do CD
total** (`CDwave` = 0,000508 contra `CD` = 0,024069), com
`Mach_dd = 0,887 > M = 0,85`.

**Hipótese a testar:** a seção da raiz, com espessura relativa de 21,8 % no
plano normal ao enflechamento e `M_n = 0,72`, apresenta choque bem mais
forte do que o modelo de Korn sugere, e o ótimo do Lab 02 se apoiou num
ponto cego do modelo. O código de Euler decide.

---

## 2. Aeronave de referência

Aeronave **B** — o joelho da frente de Pareto, escolhida pela equipe no
Lab 02 pelo critério de taxa marginal de troca.

| | |
|---|---|
| W0 | 291.291,1 kgf |
| Wf | 109.891,3 kgf |
| S_w | 368,833 m² |
| AR_w geométrico / efetivo | 9,810 / **11,772** (winglet = `True`) |
| afilamento | 0,240 |
| enflechamento c/4 | 34,651° |
| enflechamento c/2 | 32,158° |
| envergadura | 60,152 m |
| corda de raiz / ponta / MAC | 9,890 / 2,374 / 6,899 m |
| MAC em | η = 0,398 |
| tcr / tct | 0,19623 / 0,08000 |

Observação registrada porque contraria a suposição inicial da equipe: a
aeronave **tem winglet**, e o `designTool` traduz isso em `AR_eff = 1,2·AR_w`.
Isso já estava embutido no arrasto induzido que produziu a aeronave B.

---

## 3. Ponto de projeto

### 3.1 Condição da aeronave

| Parâmetro | Valor |
|---|---|
| h | 10.668 m (35.000 ft) |
| M∞ | 0,850 |
| a∞ | 296,587 m/s → V = 252,10 m/s |
| ρ∞ | 0,38046 kg/m³ |
| q∞ | 12.089,69 Pa |

**Peso escolhido: peso médio de cruzeiro.** A fração de combustível de
cruzeiro é 0,681 (alcance de 14.816 km), então o CL varia muito ao longo do
voo:

| Instante | Peso [kgf] | CL |
|---|---|---|
| início do cruzeiro | 278.385,6 | 0,612 |
| **médio: √(W_ini · W_fim)** | **229.669,3** | **0,505** |
| fim do cruzeiro: W_ini · Mf_cru | 189.478,1 | 0,417 |

O `designTool` usa o peso de **início** de cruzeiro para calcular L/D
(`weight.py:378`), mas o perfil passa a maior parte do voo perto do peso
médio. Adotamos o médio para o caso mono-ponto e cobrimos a faixa inteira na
rodada multiponto (§7), o que transforma essa escolha de aposta em resultado.

### 3.2 Transformação para a condição 2D

O perfil não vê M = 0,85. Pela teoria de asa enflechada, a seção normal vê

```
M_n     = M · cos Λ
cl_n    = cl / cos²Λ
(t/c)_n = (t/c) / cos Λ
c_n     = c · cos Λ
```

**Usamos Λ a 50 % da corda (32,158°), não a 1/4.** A justificativa é de
consistência interna: a própria equação de Korn do `designTool` usa
`sweep_50` e embute exatamente essa transformação —

```
Mach_dd = k_korn/cos Λ₅₀ − (t/c)_m/cos²Λ₅₀ − CL/(10 cos³Λ₅₀)
```

Adotar outro enflechamento aqui produziria um perfil incoerente com o modelo
que gerou a aeronave. Resulta **M_n = 0,7196** para as três estações.

### 3.3 Distribuição de sustentação

Adotamos distribuição **elíptica**, e não Schrenk nem a distribuição real da
asa sem torção. O motivo é de projeto, não de conveniência: a otimização de
torção é o passo seguinte do projeto da asa e tem a elíptica como alvo, de
modo que projetar os perfis para a elíptica é projetá-los para a asa que
pretendemos ter. Com `c·cl = K√(1−η²)` e `K = 4·CL·S/(π·b)`.

Consequência não óbvia, e importante: com afilamento 0,24 o cl local
**máximo fica em η = 0,76**, não na raiz — a corda encurta mais rápido do que
a carga. As estações externas são as mais carregadas.

---

## 4. As três estações

| Estação | η | c [m] | (t/c) | (t/c)_n | cl | cl_n | Re_n |
|---|---|---|---|---|---|---|---|
| raiz (junção asa-fuselagem) | 0,101 | 9,130 | 0,1845 | 0,2179 | 0,430 | 0,600 | 4,34×10⁷ |
| meio (MAC — a seção do roteiro) | 0,398 | 6,899 | 0,1500 | 0,1772 | 0,524 | 0,732 | 3,28×10⁷ |
| ponta | 0,900 | 3,125 | 0,0916 | 0,1082 | 0,550 | 0,768 | 1,49×10⁷ |

Justificativa das estações escolhidas:

- **Raiz em η = 0,101, não em η = 0.** A raiz geométrica está dentro da
  fuselagem e nunca vê escoamento livre. η = 0,101 é a lateral da fuselagem
  (`D_f/2 = 3,04 m`).
- **Meio na MAC.** É a seção que o roteiro pede, com `c_ref = c_MAC`.
- **Ponta em η = 0,90, não em η = 1.** Em carregamento elíptico o cl na ponta
  geométrica tende a zero (cl = 0,074 em η = 0,999) — é uma estação sem
  informação. Além disso η = 0,90 é o fim do tanque de combustível
  (`b_tank_b_w_end = 0,9`), o que dá sentido físico à restrição de espessura.

---

## 5. Formulação da otimização

Para cada estação *i*:

```
min   c_d(A_l, A_u, α)                    Euler, gradiente adjunto
s.a.  c_l(A_l, A_u, α) = cl_n,i           condição de projeto da seção
      (t/c)_max        ≥ (t/c)_n,i        função KS, já suave no airfoil_mod
      (t/c)_min        ≤ 0,01             bordo de fuga fino
      (t/c)_min        ≥ 0                superfícies não se cruzam
      t_01/√(t/c)      ≥ b_i              substituta de cl_max — §6
      batentes em A_l, A_u, α
```

Variáveis: 4 coeficientes CST do intradorso, 4 do extradorso, e α — nove no
total. Otimizador SLSQP, como no roteiro. `M_n = 0,7196` fixo.

Nota sobre a restrição de bordo de fuga: `mint ≤ 0,01` é a formulação do
professor e mantém o bordo fino; ela **não** impede o cruzamento das
superfícies (`mint < 0` a satisfaz), por isso acrescentamos `mint ≥ 0`.

### 5.1 Por que restrição rígida, e não penalidade

Consideramos três formas de tratar a sustentação máxima:

| | como | XFoil no laço | complexidade |
|---|---|---|---|
| **A** | batente na substituta geométrica | não | quase nula |
| B | substituta recalibrada por XFoil a cada iteração | 1 chamada/iteração | média |
| C | cl_max do XFoil direto como `g(x)` | sim | inviável (§6) |

**Adotamos A.** O argumento é de economia de meios: se no ótimo a restrição
sair **inativa**, ela não custou nada e não havia máquina a construir; se
sair ativa e o arrasto parecer penalizado, aí sim vale montar B. Não se
constrói máquina antes de saber se ela é necessária.

Para medir isso, cada estação é otimizada **duas vezes** — com e sem a
restrição de cl_max (`--sem-bluntez`). Se o ótimo irrestrito já atender o
limiar, a restrição era inativa. Se não atender, a diferença de `c_d` entre
as duas rodadas é o **preço exato** de manter `clmax_w = 1,80`, no mesmo
formato em que o Lab 02 reportou que "as restrições de realismo custam
+1,9 t".

Penalidade quadrática (conteúdo da Aula 03) foi considerada e descartada:
ela exigiria escolher um peso ρ, e não há taxa de câmbio derivável entre
`c_d` e `cl_max` pela aeronave — no projeto atual `dW0/dclmax_w = 0`, porque
a restrição de empuxo não está ativa (§1.1). ρ seria uma preferência
arbitrária disfarçada de número.

---

## 6. Como o cl_max entra sem quebrar o método de gradiente

Nada na formulação impede o otimizador de afiar o bordo de ataque para
reduzir arrasto de onda e destruir a sustentação máxima no caminho — o corte
controlado mostrou variação de 0,888 a 2,132 mexendo só no nariz. Então a
sustentação máxima precisa entrar de alguma forma.

A pergunta é se ela pode entrar **direto**, com o XFoil dentro do laço. O
custo não é o obstáculo: uma chamada de XFoil custa ~7 s contra ~98 s de uma
avaliação do Euler com adjunto, e nove chamadas por gradiente apenas dobrariam
o tempo por iteração. O obstáculo é a qualidade do gradiente.

**Medimos, em vez de supor** (`testa_gradiente_clmax.py`). Derivada de cl_max
por diferença central, em cinco passos:

| passo | d(cl_max)/dA_u0 | d(cl_max)/dA_u1 |
|---|---|---|
| 3×10⁻² | 0,918 | −1,237 |
| 1×10⁻² | 1,050 | −1,190 |
| 3×10⁻³ | 2,717 | 0,083 |
| 1×10⁻³ | 7,750 | 4,400 |
| 3×10⁻⁴ | **−7,667** | −0,167 |

Espalhamento de **1.617 %** e **1.491 %**, com troca de sinal. Para comparar,
o adjunto do eulerblock variou **0,04 %** entre passos.

Dois detalhes importam:

- **Não são falhas de convergência.** As 21 avaliações reportaram estol
  capturado (`queda`); a lógica de repetição de `xfoil_runner.py` já eliminou
  as falhas. O ruído é **intrínseco**, e vem da quantização do `max` sobre a
  grade discreta de α. Piso medido: ~0,005 em cl_max, contra um sinal de
  0,0006 no passo de 3×10⁻⁴ — ruído oito vezes maior que o sinal.
- **Passo grande não salva.** Os dois maiores passos concordam dentro de 13 %,
  mas 3×10⁻² num coeficiente de 0,19 é uma perturbação de 16 %: isso é uma
  secante sobre uma corda enorme, não uma derivada, e mentiria exatamente
  perto das fronteiras de restrição, que é onde o otimizador opera.

**Conclusão medida: o cl_max do XFoil não serve como `g(x)` para o SLSQP.**
É preciso um intermediário suave.

Solução em três camadas:

1. **Substituta suave** `ĉl_max(A_u, A_l)`, analítica nos coeficientes CST —
   diferenciável e de custo desprezível. Calibrada por um DOE de XFoil, e
   usada dentro da penalidade quadrática de §5.
2. **Verificação a posteriori** com o XFoil de verdade no perfil ótimo.
3. **Laço externo de correção**, se a verificação discordar da substituta:
   reajusta a substituta com o novo ponto e re-otimiza com partida quente
   (converge em poucas iterações).

### 6.1 O DOE que escolhe a substituta

Não pré-comprometemos o preditor. O DOE varre as 8 variáveis CST e registra
**doze descritores geométricos candidatos** (`descritores.py`), deixando a
regressão decidir:

| Categoria | Descritores |
|---|---|
| bordo de ataque | `r_LE` extradorso e intradorso, Δy de Abbott, espessura em x/c = 1 % e 5 % |
| espessura | t/c máximo e posição, t/c mínimo |
| arqueamento | arqueamento máximo e posição |
| carregamento traseiro | arqueamento em x/c = 0,85, inclinação da linha média no BF |

Dois deles merecem nota:

- **`r_LE/c = A[0]²/2`** — na função de classe do CST (N1 = 0,5) a superfície
  tende a `y → A[0]·√x` perto do bordo, e a parábola `y² = A²x` tem raio
  `A²/2` no vértice. Confere com o valor teórico do NACA de 4 dígitos.
- **Δy de Abbott & von Doenhoff** — diferença de ordenada do extradorso entre
  x/c = 6 % e 0,15 %, o correlator de handbook (DATCOM, Roskam) para cl_max e
  tipo de estol. Vem da mesma literatura que o `designTool` usa para flape,
  slat e fator de forma.

**O que o DOE precisa entregar.** Como a formulação usa penalidade (§5.1), e
não um batente, o produto do DOE é um **preditor do valor** de `cl_max`, não
um classificador de limiar:

```
ĉl_max(A_u, A_l) = superfície de resposta sobre k descritores
```

**Critério de aceitação**, em ordem:

1. **R² de validação cruzada ≥ 0,85** com o menor k que já sature. O R² de
   treino sempre sobe com mais termos; só o de validação diz se generaliza.
2. **Resíduo máximo < 0,15 em cl_max** na faixa que interessa (perto de
   1,80). Um preditor com viés grande ali empurraria o otimizador para o
   lugar errado sem que a penalidade percebesse.
3. **Monotonicidade no sentido físico** para o descritor dominante — um
   preditor que recompense afiar o bordo de ataque estaria invertido, e a
   verificação em XFoil pegaria isso tarde demais.

Se nenhuma superfície atingir (1), a penalidade passa a usar o descritor
dominante isolado, normalizado — pior em precisão, mas ainda no sentido
certo, e a verificação em XFoil no fim continua sendo o juiz.

A correlação de posto e a análise de limiar continuam sendo calculadas, mas
como **diagnóstico** — para saber qual descritor domina e se a separação é
limpa — e não mais como a restrição em si.

Condição do DOE: Re = 4,2×10⁷ e M = 0,266 — decolagem na MAC
(`V₂ = 1,2·V_stall = 90,59 m/s`, com `CLmaxTO = 2,2196`).

### 6.2 Resultado — corte controlado

Variando **só o nariz** a partir do NACA 1411, `cl_max` vai de **0,888 a
2,132** — fator de 2,4 — e cruza o alvo de 1,8 em **Δy = 2,77 % da corda**.
O batente do roteiro (`Au_lower = 0,05`) permite `r_LE/c = 0,125 %`, um bordo
praticamente afiado: não protege nada. Figura: `doe_clmax_corte.png`.

### 6.3 Resultado — LHS, e a substituta escolhida

256 perfis amostrados, 233 avaliados (o resto descartado por geometria
degenerada antes de gastar XFoil), **199 com estol capturado de forma
confiável**. `cl_max` de 0,685 a 2,257, mediana 1,832; 56,3 % atingem o alvo.

Correlação de posto com `cl_max`:

| descritor | ρ |
|---|---|
| `t_01` (espessura em x/c = 1 %) | **+0,756** |
| `t_05` (espessura em x/c = 5 %) | +0,736 |
| `r_LE_sup` | +0,695 |
| `delta_y` (Abbott) | +0,667 |
| `x_tmax` | −0,619 |
| demais (arqueamento, carregamento traseiro, espessura máxima) | ≤ 0,36 |

Dois resultados que contrariam a expectativa e merecem registro:

- **O vencedor é `t_01`, não o raio de bordo de ataque nem o Δy de
  handbook.** Faz sentido a posteriori: `t_01` mede a espessura efetiva do
  nariz somando as duas superfícies, enquanto `r_LE` enxerga apenas o
  primeiro coeficiente CST do extradorso e o Δy apenas o extradorso. Foi bom
  não ter pré-comprometido o preditor.
- **O Δy não admite limiar utilizável** para o alvo de 1,8 com 95 % de
  precisão, apesar de ser o correlator clássico da literatura. Ele continua
  bom como *descritor* (ρ = 0,67), mas não como *critério*.

Superfície de resposta quadrática, R² de validação cruzada por número de
descritores:

| k | 1 | 2 | 3 | 4 | 5 | **6** | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| R²_cv | 0,528 | 0,542 | 0,738 | 0,800 | 0,834 | **0,871** | 0,923 | 0,964 |
| resíduo perto de 1,8 | 0,150 | 0,146 | 0,151 | 0,130 | 0,119 | **0,099** | 0,070 | 0,041 |

**Escolha: k = 6** — `t_01`, `t_05`, `r_LE_sup`, `delta_y`, `x_tmax`,
`r_LE_inf`. É o menor k que satisfaz os dois critérios de §6.1
(R²_cv ≥ 0,85 e resíduo < 0,15 perto do alvo), com 28 termos para 199
pontos — folgado o bastante para não superajustar, e o R² de validação
cruzada confirma que generaliza. Não usamos k = 7 ou 8 porque, com apenas 8
variáveis CST no total, tantos descritores passam a reconstruir as próprias
variáveis de projeto, e o modelo perde a robustez de extrapolação fora da
caixa amostrada.

Como diagnóstico, o limiar `t_01 ≥ 0,0368` tem **95,2 % de precisão e 52,7 %
de cobertura** para o alvo de 1,8 (figura `doe_clmax_limiar.png`).

### 6.4 A substituta que ficou: `t_01/√(t/c)`, por estação

O `t_01` puro exige da ponta fina a mesma bluntez absoluta que da raiz
grossa, e a cobertura despenca para 25 % nos perfis finos. Normalizar pela
espessura inteira (`t_01/(t/c)`) é pior ainda — ρ cai de 0,756 para 0,619 —
porque o que governa o pico de sucção é a bluntez **absoluta**, não a
relativa. A raiz quadrada é o meio-termo que os dados escolheram:
ρ = 0,742, precisão 95,8 %, cobertura 60,7 %.

Continua diferenciável: `t_01` é linear nos coeficientes CST e `t/c` vem da
função KS do `airfoil_mod`, suave por construção, com gradiente já disponível
em `grads['maxt']`. O gradiente da razão é a regra do quociente, conferido
contra diferenças finitas em `test_otimiza_secao.py`.

**Limiar não-viesado por estação**, resolvendo `E[cl_max] = 1,80` na regressão
`cl_max ~ f(bluntez, t/c)` ajustada aos 277 perfis dos dois DOEs
(t/c de 0,076 a 0,270):

| estação | (t/c)_n | limiar | partida |
|---|---|---|---|
| raiz | 0,2179 | **0,1220** | 0,1301 (atende) |
| meio | 0,1772 | **0,0938** | 0,1173 (atende) |
| ponta | 0,1082 | **0,0978** | 0,0917 (abaixo) |

O limiar **não é monótono na espessura**: a raiz grossa exige nariz mais
gordo que a MAC. É a mesma física que derrubou a mediana de cl_max no DOE de
perfis grossos (1,676 contra 1,832 do LHS geral) — perfil muito espesso perde
sustentação máxima por separação de bordo de fuga, e precisa compensar no
nariz. Um limiar único, como propusemos de início, estaria errado nas duas
pontas.

**Escolha do critério: não-viesado, e não conservador.** O limiar de 95 % de
precisão rejeita ~40 % dos perfis bons — inclusive o próprio NACA 1411, cujo
cl_max medido é 1,98 — e cada perfil bom rejeitado é arrasto pago à toa. Como
o cl_max do ótimo é verificado no XFoil ao final, errar centrado no alvo é
melhor que errar sempre para o lado caro. (Sobre os dados combinados o
critério conservador sequer existe: nenhum limiar atinge 95 % de precisão,
porque a correlação cai de 0,742 para 0,518 ao incluir perfis grossos.)

**Ressalva honesta:** o R² da superfície `cl_max ~ f(bluntez, t/c)` é 0,51 —
dispersão grande. O limiar acerta a **média**, não cada perfil. Por isso a
verificação em XFoil no fim é essencial, não opcional, e o laço de correção
(apertar o limiar e re-otimizar com partida quente) faz parte do método.

---

## 7. Divisão de trabalho entre as ferramentas

Cada código faz só o que sabe fazer:

| | eulerblock (Euler 2D, adjunto) | XFoil (painéis + viscoso) |
|---|---|---|
| regime | transônico, M_n = 0,72 | baixa velocidade, M = 0,266 |
| entrega | c_d, c_l, c_m, choque, gradientes | c_l_max, polar viscosa, transição |
| papel | dentro do laço de otimização | calibração e verificação, fora do laço |
| custo | 68 s sem adjunto, 98 s com | 4–6 s por perfil |

---

## 8. Verificação da ferramenta antes de confiar nela

Antes de qualquer otimização longa:

### 8.1 Convergência de malha — feito, e o resultado é grave

NACA 1411 na condição de projeto (M_n = 0,7196, α = 2°):

| nível | células | CL | CD | erro de CD vs extrapolado |
|---|---|---|---|---|
| 0,50 | 360 | 0,5243 | 0,033519 | +359 % |
| 0,75 | 792 | 0,5500 | 0,018428 | +152 % |
| **1,00** (padrão do professor) | 1440 | 0,5626 | **0,012426** | **+70,0 %** |
| 1,25 | 2220 | 0,5688 | 0,010334 | +41,4 % |
| 1,50 | 3240 | 0,5723 | 0,009170 | +25,5 % |
| 2,00 | 5760 | 0,5777 | 0,008279 | +13,3 % |

Ajustando `CD = CD_∞ + C·h^p` aos níveis ≥ 1,0: **ordem observada p = 2,43**
(compatível com o esquema de 2ª ordem) e **CD_∞ = 0,007308**.

**O CD não está convergido nem no nível 2,0.** A malha padrão do professor
superestima o arrasto em **70 %**. Fisicamente é o esperado: numa malha
grosseira a dissipação numérica age como viscosidade artificial e gera
entropia — logo arrasto — espúria. O CL sofre muito menos (nível 1,0 erra
−5,5 %, `CL_∞ = 0,595`), porque é uma integral de pressão dominada pela
circulação global.

Três consequências para a campanha:

1. **O CD absoluto do nível 1,0 não serve para calibrar `k_korn`.** Ou se
   extrapola por Richardson, ou se roda o ponto final numa malha fina. Sem
   isso, entregaríamos ao Lab 02 um arrasto de onda 70 % maior que o real.
2. **Para a otimização, o que importa é a diferença entre projetos na mesma
   malha.** Isso costuma sobreviver ao erro de malha, mas não é garantido
   aqui: a dissipação numérica depende da intensidade do choque, que é
   exatamente o que estamos otimizando.
3. **Mitigação obrigatória:** otimizar no nível 1,0 (é o que cabe no tempo) e,
   ao final, reavaliar a linha de base **e** o ótimo nos níveis 1,5 e 2,0,
   confirmando que a melhoria Δc_d sobrevive ao refinamento. Se não
   sobreviver, o resultado da otimização é artefato de malha e precisa ser
   refeito mais fino.

### 8.2 Ainda pendente

1. **Ruído numérico.** Verificar que `c_d(α)` é suave; o gradiente exige isso.
2. **Verificação do adjunto** contra diferenças finitas em `dc_d/dα` e
   `dc_d/dA_u`. É o conteúdo da Aula 03 (FD, CS, AD, AM) aplicado, e é o tipo
   de evidência que sustenta o resto do relatório.

---

## 9. O produto sobre t/c

Uma otimização por estação responde "qual o melhor perfil **dado** este t/c".
O Lab 02 precisa de outra coisa: a curva **c_d,mín(t/c)** por estação — o
preço aerodinâmico da espessura, para pesar contra o ganho estrutural.

A restrição de espessura fica **ativa** no ótimo (mais fino é sempre melhor
para onda), então seu multiplicador de Lagrange é `d c_d,mín / d(t/c)`. Uma
otimização por estação dá o ponto e a inclinação; duas ou três
re-otimizações com partida quente em valores vizinhos de t/c fecham a curva a
custo baixo, porque cada uma parte do ótimo anterior.

---

## 10. Ordem de execução

1. Verificação do solver (§8) — ~35 min de máquina.
2. DOE de cl_max no XFoil (§6) → fixa o limiar da substituta.
3. **MAC, mono-ponto** — entregável do roteiro e diagnóstico do que o
   otimizador faz.
4. **Raiz e ponta, mono-ponto**, em paralelo.
5. **Varredura em t/c** nas três estações, com partida quente (§9).
6. **Multiponto** na faixa de CL do cruzeiro (0,417 a 0,612), que responde os
   itens 7 e 9 do roteiro.
7. **Realimentação**: `k_korn` e `clmax_w` medidos, de volta ao `designTool`.

Cada otimização roda em **pasta isolada**. Isso não é organização, é
necessidade: o `eulerblock.exe` escreve `grid.xyz`, `settings.txt`,
`wall.dat`, `solution.vtk`, `derivatives.dat` e `opt_results.pickle` com
nomes fixos no diretório corrente (`euler_mod.py:47`). Duas otimizações na
mesma pasta se corrompem mutuamente. A máquina tem 10 núcleos físicos;
usamos 6 processos simultâneos.

---

## 11. Defeitos encontrados no material fornecido

Registrados aqui porque afetam resultados do roteiro e valem menção ao
professor.

1. **`export_airfoil` quebra o XFoil em modo batch.** Com o padrão
   `close_te=True`, a função repete o ponto do bordo de fuga. O `cstfoil` já
   devolve o contorno fechado (primeiro ponto = último, que é a convenção do
   XFoil para BF afiado), então a cópia extra vira um **terceiro** ponto
   coincidente — um painel de comprimento zero. Diante disso o XFoil encerra
   **em silêncio**, com código de retorno 0 e sem mensagem de erro, logo após
   imprimir "Clockwise ordering". Afeta o `single_run.py` do próprio
   professor. Contornado em `xfoil_runner.py`, com teste de regressão.

2. **A ordem dos argumentos é invertida entre duas funções da mesma API.**
   `cstfoil(Au, Al, x)` recebe o **extradorso** primeiro;
   `run_cst(Al, Au, ...)` recebe o **intradorso** primeiro. Passar na ordem
   errada gera um perfil de dentro para fora: o gerador de malha avisa
   "negative-area cells", o solver encerra sem escrever `wall.dat` e o erro
   que chega ao Python é um `FileNotFoundError` — que não diz nada sobre a
   causa. Caímos nisso. Mitigado com uma checagem de espessura máxima logo
   após cada chamada, em `verificacao_malha.py`.

3. **`cstfoil` devolve posição errada da espessura máxima.** A linha 216
   filtra o vetor de espessuras para `x > 0,05`, mas a linha 219 indexa o
   vetor de abscissas **completo** com o índice do vetor filtrado. Para o
   NACA 1411 devolve `x/c = 0,117` quando o valor correto é 0,30. O valor de
   `max_thickness` está certo; só a posição está errada. **Afeta diretamente
   a coluna `x_t/c,max` da Tab. 2 do roteiro.** O arqueamento não sofre do
   problema, porque o vetor `cc` não é filtrado. Contornado em
   `descritores.py`, com teste que afere contra a geometria conhecida do
   NACA 1411.

---

## 12. Estrutura dos arquivos

```
otimizacao_aerofolio/
  eulerblock/ airfoil_mod/ xfoil/     material do professor, intacto
  Lab03_PRJ23_2026.pdf
  campanha/                           trabalho da equipe
    xfoil_runner.py                   interface batch com o XFoil
    descritores.py                    descritores geométricos do perfil
    doe_clmax_xfoil.py                DOE de cl_max
    analise_doe_clmax.py              escolhe a restrição substituta
    figuras_doe_clmax.py              figuras do DOE
    estilo.py                         paleta e estilo, herdados do Lab 02
    test_xfoil_runner.py              testes de regressão
    resultados/                       CSVs e figuras
```

A paleta é a mesma do Lab 02 (`#2a78d6`, `#eb6834`, `#1baf7a`), verificada
para daltonismo: pior par adjacente com ΔE 9,2 (deutan) e 27,6 (visão
normal).

---

## 13. Fora do escopo

- Otimização de torção da asa (é o passo seguinte; a hipótese de distribuição
  elíptica em §3.3 é justamente a preparação para ele).
- Otimização 3D da asa.
- Reotimização completa da aeronave no `designTool` com as constantes
  medidas — o Lab 03 entrega os números; a reotimização é decisão de projeto
  posterior.
