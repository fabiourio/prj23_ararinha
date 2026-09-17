'''
Lab 04 -- otimizacao da torcao geometrica da asa.

A torcao e parametrizada pelos valores em quatro pontos de controle
(eta = 0,3, 0,56, 0,8 e 1,0), interpolados por uma spline cubica natural
que passa pela raiz com torcao nula, avaliada nas treze estacoes da asa
do modelo. Polinomios globais (cubico e quartico) foram testados antes e
descartados: ou nao alcancam a folga de estol exigida ou precisam ondular
tanto que pioram o arrasto. Restricoes de forma completam a
parametrizacao: wash-in limitado a +2,5 graus, washout ate -10 graus e
monotonia da estacao 0,398 para fora. O problema e escolher os
coeficientes para MINIMIZAR o arrasto induzido da aeronave trimada no
ponto de projeto (M = 0,85, CL = 0,5053, arfagem compensada pelo
profundor), com duas familias de RESTRICOES: o estol pelo metodo da secao
critica em M = 0,2 deve comecar em eta <= 0,50, antes do aileron, com
folga de 0,5 grau de alpha sobre a regiao externa; e a torcao em qualquer
estacao fica na faixa de -10 a +10 graus.

A estrutura do VLM e explorada para evitar otimizacao cara:
  - a distribuicao de cl_norm e AFIM nos coeficientes e no alpha, entao o
    modelo do estol e construido com 2 rodadas de base e 1 por coeficiente;
  - o arrasto induzido e QUADRATICO nos coeficientes, entao o modelo exato
    da funcao objetivo sai de 10 rodadas no ponto de projeto.
O problema resultante e resolvido por SLSQP com multipartida e o otimo e
verificado com rodadas do AVL fora dos modelos.

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

N_COEF = 4
ETA_SEGURA = 0.50                        # estol deve comecar antes daqui
FOLGA_ALPHA = 0.5                        # [graus] de folga da regiao externa
TORCAO_MIN = -10.0                       # [graus] washout maximo por estacao
TORCAO_POS = 2.5                         # [graus] wash-in maximo por estacao
ETA_MONO = 0.398                         # daqui para fora, washout monotono
EPS_H = 4.0                              # [graus] passo do modelo quadratico

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
ETAS = np.array(gera.ETAS_ASA)

# As variaveis de projeto sao os VALORES de torcao em quatro pontos de
# controle, em graus. A curva e a spline cubica natural pelos cinco nos
# (raiz fixa em zero mais os quatro controles), avaliada em todas as
# estacoes pela matriz MAPA. A spline e linear nos valores, o que preserva
# os modelos exatos, e e local o bastante para combinar washout suave no
# meio com washout forte na ponta sem as ondulacoes de um polinomio
# global, que foram testadas e pioravam o arrasto.
from scipy.interpolate import CubicSpline  # noqa: E402

ETAS_CTRL = np.array([0.3, 0.56, 0.8, 1.0])
NOS = np.concatenate([[0.0], ETAS_CTRL])
MAPA = np.column_stack([
    CubicSpline(NOS, np.concatenate([[0.0], e]), bc_type='natural')(ETAS)
    for e in np.eye(N_COEF)])

with open('avl/saidas/trim.json', encoding='ascii') as f:
    TRIM = json.load(f)


def torcoes_de(t):
    valores = MAPA @ np.asarray(t)
    return {float(eta): float(v) for eta, v in zip(ETAS, valores)}


def cdff_cruzeiro(c):
    '''CDff (Trefftz) da aeronave trimada no ponto de projeto.'''
    gera.gera_variante(VARIANTE, torcoes_de(c), cg='aft')
    saida = roda_avl('load saidas/tw.avl\noper\n'
                     f'm\nmn {MACH_CRU}\n\n'
                     f'de\n1 {TRIM["aft"]["it_deg"]}\n\n'
                     f'd2 pm 0\na c {CL_PROJ}\nx\n\nquit\n')
    return pega(saida, r'CDff\s*=\s*([-\d.Ee+]+)')


def estol_amostras(c, alfas):
    '''cl_norm por faixa da asa direita e CLtot em cada alpha, M = 0,2.'''
    gera.gera_variante(VARIANTE, torcoes_de(c), cg='fwd')
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


CACHE = 'avl/saidas/torcao_modelo.npz'

# ------------------------------------------------------------------
# MODELO AFIM DO ESTOL: cl_norm_i = a_i + b_i alpha + sum c_j f_ji
if os.path.exists(CACHE):
    print(f'Usando modelos do cache {CACHE}')
    _m = np.load(CACHE)
    ETA_STRIP, A_STRIP, B_STRIP = _m['eta'], _m['a'], _m['b']
    F_STRIP, A_CL, B_CL, D_CL = _m['f'], _m['acl'], _m['bcl'], _m['dcl']
    f0, G, H = float(_m['f0']), _m['g'], _m['h']
else:
    print('Construindo o modelo afim do estol (M = 0,2)...')
    (d8, d14), (cl8, cl14) = estol_amostras(np.zeros(N_COEF), ALFAS_BASE)
    ETA_STRIP = d8[0]
    B_STRIP = (d14[1] - d8[1])/(ALFAS_BASE[1] - ALFAS_BASE[0])
    A_STRIP = d8[1] - B_STRIP*ALFAS_BASE[0]
    B_CL = (cl14 - cl8)/(ALFAS_BASE[1] - ALFAS_BASE[0])
    A_CL = cl8 - B_CL*ALFAS_BASE[0]

    F_STRIP = np.zeros((N_COEF, len(ETA_STRIP)))
    D_CL = np.zeros(N_COEF)
    for j in range(N_COEF):
        c = np.zeros(N_COEF)
        c[j] = 1.0
        (dj,), (clj,) = estol_amostras(c, ALFAS_BASE[:1])
        F_STRIP[j] = dj[1] - d8[1]
        D_CL[j] = clj - cl8
        print(f'  sensibilidade do coeficiente c{j + 1}')

LIMITE_STRIP = np.interp(ETA_STRIP, ETA_LIM, CLMAX_LIM)
INTERNA = ETA_STRIP <= ETA_SEGURA
EXTERNA = ~INTERNA


def alfas_de_estol(c):
    cln0 = A_STRIP + F_STRIP.T @ np.asarray(c)
    return (LIMITE_STRIP - cln0)/B_STRIP


def folga_estol(c):
    '''Quanto a regiao externa estola depois da interna [graus].'''
    alfas = alfas_de_estol(c)
    return np.min(alfas[EXTERNA]) - np.min(alfas[INTERNA])


FOLGA_USADA = FOLGA_ALPHA                # pode ser adaptada na fase de
                                         # factibilidade, com aviso


def restricao_estol(c):
    '''>= 0 quando a regiao interna estola antes da externa com folga.'''
    return folga_estol(c) - FOLGA_USADA


K_MONO = int(np.searchsorted(ETAS, ETA_MONO))


def restricao_faixa(t):
    '''>= 0 quando a torcao de toda estacao esta na faixa admissivel e o
    washout e monotono da estacao ETA_MONO para fora.'''
    valores = MAPA @ np.asarray(t)
    monotonia = valores[K_MONO:-1] - valores[K_MONO + 1:]
    return np.concatenate([TORCAO_POS - valores, valores - TORCAO_MIN,
                           monotonia])


def clmax_previsto(c):
    alfas = alfas_de_estol(c)
    a_estol = np.min(alfas)
    return A_CL + B_CL*a_estol + np.asarray(c) @ D_CL, a_estol


if not os.path.exists(CACHE):
    # --------------------------------------------------------------
    # MODELO QUADRATICO DO CDff NO CRUZEIRO. O passo EPS_H largo dilui o
    # ruido de leitura da saida do AVL sem erro de modelo, ja que o CDff
    # e exatamente quadratico nas torcoes.
    print('Construindo o modelo quadratico do CDff (M = 0,85)...')
    f0 = cdff_cruzeiro(np.zeros(N_COEF))
    print(f'  CDff da asa sem torcao: {f0:.6f}')
    f_mais, f_menos = np.zeros(N_COEF), np.zeros(N_COEF)
    for j in range(N_COEF):
        e = np.zeros(N_COEF)
        e[j] = EPS_H
        f_mais[j] = cdff_cruzeiro(e)
        f_menos[j] = cdff_cruzeiro(-e)
        print(f'  sensibilidade da variavel t{j + 1}')
    H = np.zeros((N_COEF, N_COEF))
    G = (f_mais - f_menos)/(2*EPS_H)
    for j in range(N_COEF):
        H[j, j] = (f_mais[j] + f_menos[j] - 2*f0)/EPS_H**2
    for j, k in itertools.combinations(range(N_COEF), 2):
        e = np.zeros(N_COEF)
        e[j] = e[k] = EPS_H
        fjk = cdff_cruzeiro(e)
        H[j, k] = H[k, j] = (fjk - f_mais[j] - f_mais[k] + f0)/EPS_H**2
        print(f'  termo cruzado t{j + 1} x t{k + 1}')

    np.savez(CACHE, eta=ETA_STRIP, a=A_STRIP, b=B_STRIP, f=F_STRIP,
             acl=A_CL, bcl=B_CL, dcl=D_CL, f0=f0, g=G, h=H)
    print(f'modelos gravados em {CACHE}')


def objetivo(c):
    t = np.asarray(c)
    return f0 + G @ t + 0.5*t @ H @ t


# ------------------------------------------------------------------
# FASE DE FACTIBILIDADE: qual a maior folga que a parametrizacao alcanca?
rng = np.random.default_rng(23)
melhor_folga, x_folga = -np.inf, None
for x0 in [np.zeros(N_COEF)] + [rng.uniform(-9, 5, N_COEF)
                                for _ in range(20)]:
    r = minimize(lambda c: -folga_estol(c), x0, method='SLSQP',
                 bounds=[(TORCAO_MIN, TORCAO_POS)]*N_COEF,
                 constraints=[{'type': 'ineq', 'fun': restricao_faixa}],
                 options={'maxiter': 300, 'ftol': 1e-9})
    if (restricao_faixa(r.x).min() > -1e-6
            and folga_estol(r.x) > melhor_folga):
        melhor_folga, x_folga = folga_estol(r.x), r.x.copy()
print(f'Folga maxima alcancavel pela parametrizacao: '
      f'{melhor_folga:+.3f} graus')
if melhor_folga < FOLGA_ALPHA:
    FOLGA_USADA = max(0.2, melhor_folga - 0.05)
    print(f'AVISO: alvo de folga adaptado de {FOLGA_ALPHA} para '
          f'{FOLGA_USADA:.2f} graus')

# ------------------------------------------------------------------
# OTIMIZACAO COM MULTIPARTIDA (objetivo escalado para drag counts)
print('Otimizando...')
melhor_x, melhor_f = None, np.inf
partidas = [np.zeros(N_COEF), x_folga]
partidas += [rng.uniform(-8, 3, N_COEF) for _ in range(10)]
vinculos = [{'type': 'ineq', 'fun': restricao_estol},
            {'type': 'ineq', 'fun': restricao_faixa}]
for x0 in partidas:
    r = minimize(lambda c: 1e4*objetivo(c), x0, method='SLSQP',
                 bounds=[(TORCAO_MIN, TORCAO_POS)]*N_COEF,
                 constraints=vinculos,
                 options={'maxiter': 500, 'ftol': 1e-8})
    viavel = (restricao_estol(r.x) > -1e-6
              and restricao_faixa(r.x).min() > -1e-6)
    print(f'  partida {np.round(x0, 1)}: f = {objetivo(r.x):.6f}, '
          f'viavel = {viavel}, sucesso = {r.success} ({r.message})')
    if viavel and objetivo(r.x) < melhor_f:
        melhor_x, melhor_f = r.x.copy(), objetivo(r.x)
if melhor_x is None:
    raise RuntimeError('nenhuma partida terminou viavel')
c_otimo = melhor_x

# Verificacao com o AVL de verdade (fora dos modelos)
cdff_ver = cdff_cruzeiro(c_otimo)
clmax_mod, a_estol_mod = clmax_previsto(c_otimo)
os.remove(VARIANTE)

torcoes = torcoes_de(c_otimo)
print('\nTorcoes nos pontos de controle: '
      + ', '.join(f'eta {e} = {v:+.3f}'
                  for e, v in zip(ETAS_CTRL, c_otimo)))
print('Torcoes por estacao [graus]:')
for eta in gera.ETAS_ASA:
    print(f'  eta = {eta}: {torcoes[eta]:+.3f}')
print(f'CDff sem torcao : {f0:.6f}')
print(f'CDff otimizado (modelo): {objetivo(c_otimo):.6f}')
print(f'CDff otimizado (AVL)   : {cdff_ver:.6f}')
print(f'restricao de estol (folga em alpha) = '
      f'{restricao_estol(c_otimo):+.3f} graus')
print(f'alpha de estol previsto = {a_estol_mod:.2f} graus, '
      f'CLmax previsto = {clmax_mod:.3f}')

resultado = {'controle': {str(e): round(float(v), 3)
                          for e, v in zip(ETAS_CTRL, c_otimo)},
             'torcoes': {str(eta): round(torcoes[eta], 3)
                         for eta in gera.ETAS_ASA},
             'CDff_sem_torcao': round(float(f0), 6),
             'CDff_otimizado_avl': round(float(cdff_ver), 6),
             'folga_restricao_graus':
                 round(float(restricao_estol(c_otimo)), 3)}
with open('avl/saidas/torcao.json', 'w', encoding='ascii') as f:
    json.dump(resultado, f, indent=2)
print('\ngravado: avl/saidas/torcao.json')
print('agora rode: lab04_gera_avl.py, lab04_trim.py, lab04_secao_critica.py,')
print('lab04_polares.py, lab04_derivadas.py e as figuras')
