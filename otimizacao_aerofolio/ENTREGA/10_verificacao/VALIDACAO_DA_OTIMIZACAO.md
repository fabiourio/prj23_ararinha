# Validacao da otimizacao da MAC -- Lab 03

Duas verificacoes que o Lab 02 fez para a aeronave e que faltavam aqui.

## 1. O otimo e global? (multistart)

Seis pontos de partida geometricamente distintos, todos viaveis: espessura
maxima de x/c = 0,285 a 0,376, arqueamento de +0,008 a +0,038, c_d inicial
variando por um fator 2,5.

| partida | c_d inicial | c_d final | c_l final | avaliacoes |
|---|---|---|---|---|
| ms0 (NACA 1411, o do roteiro) | 0,055979 | 0,011644 | 0,73189 | 31 |
| ms1 | 0,059049 | 0,011643 | 0,73188 | 136 |
| ms2 | 0,039666 | 0,011645 | 0,73190 | 26 |
| ms3 | 0,097431 | 0,011644 | 0,73189 | 17 |
| ms4 | 0,091078 | 0,011642 | 0,73190 | 19 |
| ms5 | 0,087456 | 0,011643 | 0,73190 | 20 |

Espalhamento entre os seis: **2,6x10^-6**, seis vezes MENOR que o piso de
ruido do solver (1,5x10^-5). Nao e "parecido" -- e o mesmo ponto, dentro da
precisao que o Euler entrega. Evidencia de bacia unica dominante.

**Uma ressalva de leitura.** O c_d so pode ser comparado entre pontos que
estejam na MESMA sustentacao. Tomando como otimo de cada rodada o melhor
ponto com |c_l - 0,7319| < 5x10^-3 -- a tolerancia larga --, o espalhamento
sobe para 2,3x10^-5, e ms1 e ms2 parecem 1,5x10^-5 melhores que os demais.
Isso nao e otimo diferente: e o c_d medido num c_l ate 0,4% mais baixo.
Apertando a tolerancia para 1x10^-3 os seis colapsam nos 2,6x10^-6 da tabela.
A coluna ``cd_final_tol_5e-3`` do CSV guarda as duas leituras.

## 2. Por que o otimo para ali? (cortes 1-a-1)

Cada um dos oito coeficientes CST deslocado de +-0,06 em torno do otimo, com
o alpha re-trimado para segurar o c_l no alvo. 48 direcoes, 44 avaliadas.

**Nenhuma direcao melhora o arrasto sem violar restricao.**

| o que barra a direcao | direcoes |
|---|---|
| espessura minima | 24 |
| bluntez (substituta de cl_max) | 13 |
| nada -- simplesmente piora o arrasto | 7 |

A restricao de ESPESSURA e quem manda na forma: barra mais da metade das
direcoes. A de sustentacao maxima e a segunda. O perfil obtido e o minimo de
arrasto possivel dada a espessura estrutural que a asa exige.

Dados: ``cortes_meio.csv``, ``multistart_meio.csv``

## 3. O poco de arrasto e estreito de verdade? (varredura fina + malha)

A polar do item 7 mostra c_d caindo 40% num unico ponto -- o de projeto. Tres
verificacoes para separar fisica de artefato:

**Varredura com passo de 0,1 grau** (5x mais fino): o poco existe e os pontos
que o cercam convergiram com residuo ~1e-9. A +0,2 grau o c_d ja sobe 26%; a
-0,2 grau, 50%.

**Residuo dos pontos suspeitos**: 6 de 21 pontos da varredura estouraram as
20000 iteracoes. Reavaliados com CFL menor ate convergir, o c_d muda no
maximo 0,04% -- o residuo estagnou em ciclo-limite em torno do choque, mas
sobre o valor certo. A polar da entrega esta correta.

**Refinamento de malha** (nivel 1,5): o poco sobrevive e APROFUNDA (razao
c_d_max/c_d_min na faixa de +-0,4 grau: 1,60x -> 1,89x), e o minimo continua
no ponto de projeto. O poco e fisica de otimizacao mono-ponto, nao artefato.

Leitura de projeto: o perfil e excelente NO ponto de projeto e degrada rapido
fora dele -- assinatura classica de otimizacao mono-ponto. A proxima iteracao
natural e a otimizacao multiponto. E o c_l de projeto depende da hipotese de
distribuicao eliptica, que o Lab 04 (torcao) precisa confirmar.

Dados: ``poco_de_arrasto_meio.csv``, ``poco_malha_meio.csv``

## 4. Alpha de trimagem da partida (item 3 do roteiro)

NACA 1411 reescalado para t/c = 0,1772, trimado ao c_l de projeto por Newton
com a derivada do adjunto (3 avaliacoes):

    alpha = 3,567 graus   c_l = 0,7319   c_d = 0,06705   c_m = -0,0680

Na MESMA sustentacao, o ganho da otimizacao e 0,06705 -> 0,01164 = -82,6%
(o -79,2% citado antes comparava a partida em alpha = 3 graus, com c_l menor).
