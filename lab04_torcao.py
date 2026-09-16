'''
Lab 04 -- otimizacao da torcao geometrica da asa.

Problema: escolher as torcoes das estacoes eta = 0,398, 0,56, 0,90 e 1,0
(raiz fixa em zero como referencia) para MINIMIZAR o arrasto induzido da
aeronave trimada no ponto de projeto (M = 0,85, CL = 0,5053, arfagem
compensada pelo profundor), com a RESTRICAO de que o estol pelo metodo da
secao critica em M = 0,2 comece em eta <= 0,50, antes do aileron (que
ocupa de 0,56 a 0,90), com folga de 0,5 grau de alpha sobre a regiao
externa.

A estrutura do VLM e explorada para evitar otimizacao cara:
  - a distribuicao de cl_norm e AFIM nas torcoes e no alpha, entao o
    modelo do estol e construido com 2 rodadas de base e 1 por variavel;
  - o arrasto induzido e QUADRATICO nas torcoes, entao o modelo exato da
    funcao objetivo sai de 15 rodadas (base, +-1 grau por variavel e os
    pares) no ponto de projeto.
O problema resultante e resolvido por SLSQP com multipartida.

Gera avl/saidas/torcao.json, que o lab04_gera_avl.py le automaticamente
ao regenerar os arquivos oficiais.

Rodar da raiz do repo:  python lab04_torcao.py
(depois: lab04_gera_avl.py e o restante do pipeline, nesta ordem)
'''

# IMPORTS
import itertools
import json
import os
import re

import numpy as np
from scipy.optimize import minimize

from avl_batch import roda_avl, pega, pega_todos
import lab04_gera_avl as gera

#=========================================

ETAS_VAR = [0.398, 0.56, 0.90, 1.0]
ETA_SEGURA = 0.50                        # estol deve comecar antes daqui
FOLGA_ALPHA = 0.5                        # [graus] de folga da regiao externa
LIMITES = [(-6.0, 3.0)]*len(ETAS_VAR)

MACH_CRU = 0.85
CL_PROJ = 0.5053
MACH_BAIXO = 0.2
ALFAS_BASE = (8.0, 14.0)
SEMI_ENV = 30.0759

ETA_LIM = np.array([0.1011, 0.398, 0.90])
CLMAX_LIM = np.array([1.774, 1.7985, 1.7338])

RE_LINHA = re.compile(r'^\s*\d+\s+(' + r'[-\dEe.+]+\s+'*11 +
                      r'[-\dEe.+]+)\s*$', re.M)

VARIANTE = 'avl/saidas/tw.avl'

with open('avl/saidas/trim.json', encoding='ascii') as f:
    TRIM = json.load(f)


def torcoes_de(theta):
    d = {0.0: 0.0, 0.1011: 0.0}
    d.update({eta: float(t) for eta, t in zip(ETAS_VAR, theta)})
    return d


def cdff_cruzeiro(theta):
    '''CDff (Trefftz) da aeronave trimada no ponto de projeto.'''
    gera.gera_variante(VARIANTE, torcoes_de(theta), cg='aft')
    saida = roda_avl('load saidas/tw.avl\noper\n'
                     f'm\nmn {MACH_CRU}\n\n'
                     f'de\n1 {TRIM["aft"]["it_deg"]}\n\n'
                     f'd2 pm 0\na c {CL_PROJ}\nx\n\nquit\n')
    return pega(saida, r'CDff\s*=\s*([-\d.Ee+]+)')


def estol_amostras(theta, alfas):
    '''cl_norm por faixa da asa direita e CLtot em cada alpha, M = 0,2.'''
    gera.gera_variante(VARIANTE, torcoes_de(theta), cg='fwd')
    cmds = ('load saidas/tw.avl\noper\n'
            f'm\nmn {MACH_BAIXO}\n\n'
            f'de\n1 {TRIM["fwd"]["it_deg"]}\n\nd2 d2 0\n')
    for alfa in alfas:
        cmds += f'a a {alfa}\nx\nfs\n\n'
    cmds += '\nquit\n'
    saida = roda_avl(cmds)
    dists, cls = [], []
    for bloco in saida.split('Vortex Lattice Output')[1:]:
        cls.append(pega_todos(bloco, r'CLtot =\s+([-\d.]+)')[0])
        parte = bloco.split('Surface # 1     Wing')[1].split('Surface # 2')[0]
        dados = np.array([[float(v) for v in m.group(1).split()]
                          for m in RE_LINHA.finditer(parte)])
        dists.append((dados[:, 0]/SEMI_ENV, dados[:, 5]))
    return dists, cls


# ------------------------------------------------------------------
# MODELO AFIM DO ESTOL: cl_norm_i = a_i + b_i alpha + sum theta_j c_ji
print('Construindo o modelo afim do estol (M = 0,2)...')
(d8, d14), (cl8, cl14) = estol_amostras(np.zeros(4), ALFAS_BASE)
ETA_STRIP = d8[0]
B_STRIP = (d14[1] - d8[1])/(ALFAS_BASE[1] - ALFAS_BASE[0])
A_STRIP = d8[1] - B_STRIP*ALFAS_BASE[0]
B_CL = (cl14 - cl8)/(ALFAS_BASE[1] - ALFAS_BASE[0])
A_CL = cl8 - B_CL*ALFAS_BASE[0]

C_STRIP = np.zeros((len(ETAS_VAR), len(ETA_STRIP)))
D_CL = np.zeros(len(ETAS_VAR))
for j in range(len(ETAS_VAR)):
    theta = np.zeros(4)
    theta[j] = 1.0
    (dj,), (clj,) = estol_amostras(theta, ALFAS_BASE[:1])
    C_STRIP[j] = dj[1] - d8[1]
    D_CL[j] = clj - cl8
    print(f'  sensibilidade da torcao em eta = {ETAS_VAR[j]}')

LIMITE_STRIP = np.interp(ETA_STRIP, ETA_LIM, CLMAX_LIM)
INTERNA = ETA_STRIP <= ETA_SEGURA
EXTERNA = ~INTERNA


def alfas_de_estol(theta):
    '''alpha em que cada faixa alcanca o limite de clmax.'''
    cln0 = A_STRIP + C_STRIP.T @ np.asarray(theta)
    return (LIMITE_STRIP - cln0)/B_STRIP


def restricao(theta):
    '''>= 0 quando a regiao interna estola antes da externa com folga.'''
    alfas = alfas_de_estol(theta)
    return np.min(alfas[EXTERNA]) - np.min(alfas[INTERNA]) - FOLGA_ALPHA


def clmax_previsto(theta):
    alfas = alfas_de_estol(theta)
    a_estol = np.min(alfas)
    return A_CL + B_CL*a_estol + np.asarray(theta) @ D_CL, a_estol


# ------------------------------------------------------------------
# MODELO QUADRATICO DO CDff NO CRUZEIRO
print('Construindo o modelo quadratico do CDff (M = 0,85)...')
n = len(ETAS_VAR)
f0 = cdff_cruzeiro(np.zeros(n))
print(f'  CDff da asa sem torcao: {f0:.6f}')
f_mais, f_menos = np.zeros(n), np.zeros(n)
for j in range(n):
    e = np.zeros(n)
    e[j] = 1.0
    f_mais[j] = cdff_cruzeiro(e)
    f_menos[j] = cdff_cruzeiro(-e)
    print(f'  sensibilidade da torcao em eta = {ETAS_VAR[j]}')
H = np.zeros((n, n))
G = (f_mais - f_menos)/2
for j in range(n):
    H[j, j] = f_mais[j] + f_menos[j] - 2*f0
for j, k in itertools.combinations(range(n), 2):
    e = np.zeros(n)
    e[j] = e[k] = 1.0
    fjk = cdff_cruzeiro(e)
    H[j, k] = H[k, j] = fjk - f_mais[j] - f_mais[k] + f0
    print(f'  termo cruzado eta = {ETAS_VAR[j]} x {ETAS_VAR[k]}')


def objetivo(theta):
    t = np.asarray(theta)
    return f0 + G @ t + 0.5*t @ H @ t


# ------------------------------------------------------------------
# OTIMIZACAO COM MULTIPARTIDA
print('Otimizando...')
rng = np.random.default_rng(23)
melhor = None
partidas = [np.zeros(n), np.array([0.0, -0.5, -2.0, -3.0])]
partidas += [rng.uniform(-4, 1, n) for _ in range(6)]
for x0 in partidas:
    r = minimize(objetivo, x0, method='SLSQP', bounds=LIMITES,
                 constraints=[{'type': 'ineq', 'fun': restricao}],
                 options={'maxiter': 200, 'ftol': 1e-10})
    if r.success and restricao(r.x) > -1e-6:
        if melhor is None or r.fun < melhor.fun:
            melhor = r
if melhor is None:
    raise RuntimeError('nenhuma partida convergiu com a restricao ativa')
theta_otimo = melhor.x

# Verificacao com o AVL de verdade (fora dos modelos)
cdff_ver = cdff_cruzeiro(theta_otimo)
clmax_mod, a_estol_mod = clmax_previsto(theta_otimo)
os.remove(VARIANTE)

print('\nTorcoes otimas [graus]:')
for eta, t in zip(ETAS_VAR, theta_otimo):
    print(f'  eta = {eta}: {t:+.3f}')
print(f'CDff sem torcao : {f0:.6f}')
print(f'CDff otimizado (modelo): {objetivo(theta_otimo):.6f}')
print(f'CDff otimizado (AVL)   : {cdff_ver:.6f}')
print(f'restricao (folga em alpha) = {restricao(theta_otimo):+.3f} graus')
print(f'alpha de estol previsto = {a_estol_mod:.2f} graus, '
      f'CLmax previsto = {clmax_mod:.3f}')

resultado = {'torcoes': {'0.0': 0.0, '0.1011': 0.0},
             'CDff_sem_torcao': round(float(f0), 6),
             'CDff_otimizado_avl': round(float(cdff_ver), 6),
             'folga_restricao_graus': round(float(restricao(theta_otimo)), 3)}
for eta, t in zip(ETAS_VAR, theta_otimo):
    resultado['torcoes'][str(eta)] = round(float(t), 3)
with open('avl/saidas/torcao.json', 'w', encoding='ascii') as f:
    json.dump(resultado, f, indent=2)
print('\ngravado: avl/saidas/torcao.json')
print('agora rode: lab04_gera_avl.py, lab04_trim.py, lab04_secao_critica.py,')
print('lab04_polares.py, lab04_derivadas.py e as figuras')
