'''
Lab 04 -- otimizacao da torcao geometrica da asa.

A torcao e parametrizada pelos valores em quatro pontos de controle
(eta = 0,3, 0,56, 0,8 e 1,0), interpolados por uma spline cubica natural
que passa pela raiz com torcao nula, avaliada nas treze estacoes da asa
do modelo. Polinomios globais (cubico e quartico) foram testados antes e
descartados: ou nao alcancam a folga de estol exigida ou precisam ondular
tanto que pioram o arrasto. Restricoes de forma completam a parametrizacao:
wash-in limitado a +2,5 graus, washout ate -10 graus e monotonia da estacao
0,398 para fora.

O problema e escolher as torcoes que MINIMIZAM o arrasto induzido da
aeronave trimada no ponto de projeto (M = 0,85, CL = 0,5053, arfagem
compensada pelo profundor), com a RESTRICAO de que o estol pelo metodo da
secao critica em M = 0,2 comece em eta <= 0,50, antes do aileron.

Como a incidencia de empenagem que trima a aeronave (it) depende da torcao,
e o estol e o arrasto dependem de it, o problema e um ponto fixo: para cada
torcao candidata, a aeronave e re-trimada antes de avaliar as restricoes e
o objetivo. O laco externo alterna otimizacao e re-trimagem ate a torcao
convergir, de modo que o otimo entregue e consistente com a trimagem final.

A estrutura do VLM e explorada para evitar otimizacao cara: para uma
trimagem fixa, a distribuicao de cl_norm e afim nas torcoes e o arrasto
induzido e quadratico, entao cada iteracao externa constroi modelos exatos
com poucas rodadas e resolve por SLSQP com multipartida.

Gera avl/saidas/torcao.json, que o lab04_gera_avl.py le automaticamente.

Rodar da raiz do repo:  python lab04_torcao.py
(depois: lab04_gera_avl.py e o restante do pipeline, nesta ordem)
'''

# IMPORTS
import itertools
import json
import os
import re

import numpy as np
from scipy.interpolate import CubicSpline
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
LOAD_REL = 'saidas/tw.avl'
ETAS = np.array(gera.ETAS_ASA)
ETAS_CTRL = np.array([0.3, 0.56, 0.8, 1.0])
NOS = np.concatenate([[0.0], ETAS_CTRL])
MAPA = np.column_stack([
    CubicSpline(NOS, np.concatenate([[0.0], e]), bc_type='natural')(ETAS)
    for e in np.eye(N_COEF)])
K_MONO = int(np.searchsorted(ETAS, ETA_MONO))
LIMITE_STRIP = None                      # definido ao construir o modelo


def torcoes_de(t):
    valores = MAPA @ np.asarray(t)
    return {float(eta): float(v) for eta, v in zip(ETAS, valores)}


#=========================================
# ROTINAS DO AVL (parametrizadas pela incidencia de empenagem)

def trima(cg, ctrl):
    '''Incidencia de empenagem (it) que anula a deflexao de profundor em
    cruzeiro para a torcao dada, pelo mesmo metodo do lab04_trim.'''
    def delta_e(it):
        gera.gera_variante(VARIANTE, torcoes_de(ctrl), cg=cg)
        saida = roda_avl(f'load {LOAD_REL}\noper\n'
                         f'm\nmn {MACH_CRU}\n\n'
                         f'a c {CL_PROJ}\nd2 pm 0\n'
                         f'de\n1 {it}\n\nx\n\nquit\n')
        return pega(saida, r'elevator\s+=\s+([-\d.]+)')
    d0, d1 = delta_e(0.0), delta_e(-3.0)
    return 0.0 - d0*(-3.0)/(d1 - d0)


def cdff_cruzeiro(ctrl, it_aft):
    '''CDff (Trefftz) da aeronave trimada no ponto de projeto.'''
    gera.gera_variante(VARIANTE, torcoes_de(ctrl), cg='aft')
    saida = roda_avl(f'load {LOAD_REL}\noper\n'
                     f'm\nmn {MACH_CRU}\n\n'
                     f'de\n1 {it_aft}\n\n'
                     f'd2 pm 0\na c {CL_PROJ}\nx\n\nquit\n')
    return pega(saida, r'CDff\s*=\s*([-\d.Ee+]+)')


def estol_amostras(ctrl, alfas, it_fwd):
    '''cl_norm por faixa da asa direita e CLtot em cada alpha, M = 0,2.'''
    gera.gera_variante(VARIANTE, torcoes_de(ctrl), cg='fwd')
    cmds = (f'load {LOAD_REL}\noper\n'
            f'm\nmn {MACH_BAIXO}\n\n'
            f'de\n1 {it_fwd}\n\nd2 d2 0\n')
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


#=========================================
# MODELOS EXATOS PARA UMA TRIMAGEM FIXA

def constroi_modelos(it_fwd, it_aft):
    '''Modelo afim do estol e quadratico do arrasto para os it dados.'''
    (d8, d14), _ = estol_amostras(np.zeros(N_COEF), ALFAS_BASE, it_fwd)
    eta = d8[0]
    b_strip = (d14[1] - d8[1])/(ALFAS_BASE[1] - ALFAS_BASE[0])
    a_strip = d8[1] - b_strip*ALFAS_BASE[0]
    f_strip = np.zeros((N_COEF, len(eta)))
    for j in range(N_COEF):
        c = np.zeros(N_COEF)
        c[j] = 1.0
        (dj,), _ = estol_amostras(c, ALFAS_BASE[:1], it_fwd)
        f_strip[j] = dj[1] - d8[1]

    f0 = cdff_cruzeiro(np.zeros(N_COEF), it_aft)
    f_mais, f_menos = np.zeros(N_COEF), np.zeros(N_COEF)
    for j in range(N_COEF):
        e = np.zeros(N_COEF)
        e[j] = EPS_H
        f_mais[j] = cdff_cruzeiro(e, it_aft)
        f_menos[j] = cdff_cruzeiro(-e, it_aft)
    H = np.zeros((N_COEF, N_COEF))
    G = (f_mais - f_menos)/(2*EPS_H)
    for j in range(N_COEF):
        H[j, j] = (f_mais[j] + f_menos[j] - 2*f0)/EPS_H**2
    for j, k in itertools.combinations(range(N_COEF), 2):
        e = np.zeros(N_COEF)
        e[j] = e[k] = EPS_H
        fjk = cdff_cruzeiro(e, it_aft)
        H[j, k] = H[k, j] = (fjk - f_mais[j] - f_mais[k] + f0)/EPS_H**2
    return {'eta': eta, 'a': a_strip, 'b': b_strip, 'f': f_strip,
            'f0': f0, 'G': G, 'H': H}


def otimiza(mod):
    '''Resolve o problema com restricao sobre os modelos dados.'''
    eta, a, b, f = mod['eta'], mod['a'], mod['b'], mod['f']
    f0, G, H = mod['f0'], mod['G'], mod['H']
    limite = np.interp(eta, ETA_LIM, CLMAX_LIM)
    interna, externa = eta <= ETA_SEGURA, eta > ETA_SEGURA

    def folga(t):
        alfas = (limite - a - f.T @ np.asarray(t))/b
        return np.min(alfas[externa]) - np.min(alfas[interna])

    def faixa(t):
        v = MAPA @ np.asarray(t)
        mono = v[K_MONO:-1] - v[K_MONO + 1:]
        return np.concatenate([TORCAO_POS - v, v - TORCAO_MIN, mono])

    def obj(t):
        t = np.asarray(t)
        return 1e4*(f0 + G @ t + 0.5*t @ H @ t)

    # folga alcancavel pela parametrizacao
    rng = np.random.default_rng(23)
    melhor_folga, x_folga = -np.inf, np.zeros(N_COEF)
    for x0 in [np.zeros(N_COEF)] + [rng.uniform(-9, 2.5, N_COEF)
                                    for _ in range(20)]:
        r = minimize(lambda t: -folga(t), x0, method='SLSQP',
                     bounds=[(TORCAO_MIN, TORCAO_POS)]*N_COEF,
                     constraints=[{'type': 'ineq', 'fun': faixa}],
                     options={'maxiter': 300, 'ftol': 1e-9})
        if faixa(r.x).min() > -1e-6 and folga(r.x) > melhor_folga:
            melhor_folga, x_folga = folga(r.x), r.x.copy()
    folga_alvo = FOLGA_ALPHA if melhor_folga >= FOLGA_ALPHA \
        else max(0.2, melhor_folga - 0.05)

    vinc = [{'type': 'ineq', 'fun': lambda t: folga(t) - folga_alvo},
            {'type': 'ineq', 'fun': faixa}]
    melhor_x, melhor_f = None, np.inf
    for x0 in [np.zeros(N_COEF), x_folga] + [rng.uniform(-8, 2.5, N_COEF)
                                             for _ in range(10)]:
        r = minimize(obj, x0, method='SLSQP',
                     bounds=[(TORCAO_MIN, TORCAO_POS)]*N_COEF,
                     constraints=vinc, options={'maxiter': 500, 'ftol': 1e-8})
        if (folga(r.x) - folga_alvo > -1e-6 and faixa(r.x).min() > -1e-6
                and obj(r.x) < melhor_f):
            melhor_x, melhor_f = r.x.copy(), obj(r.x)
    if melhor_x is None:
        raise RuntimeError('nenhuma partida terminou viavel')
    return melhor_x, folga_alvo, melhor_folga


#=========================================
# LACO DE PONTO FIXO: torcao <-> trimagem

with open('avl/saidas/trim.json', encoding='ascii') as f:
    trim = json.load(f)
it_fwd, it_aft = trim['fwd']['it_deg'], trim['aft']['it_deg']

MAX_ITER = 6
TOL = 0.05                                # graus, convergencia da torcao
ctrl = np.zeros(N_COEF)
historico = []
for it_ext in range(MAX_ITER):
    print(f'\n===== iteracao externa {it_ext} '
          f'(it_fwd = {it_fwd:.3f}, it_aft = {it_aft:.3f}) =====')
    mod = constroi_modelos(it_fwd, it_aft)
    ctrl_novo, folga_alvo, folga_max = otimiza(mod)
    tors_novo = MAPA @ ctrl_novo

    it_fwd = trima('fwd', ctrl_novo)
    it_aft = trima('aft', ctrl_novo)

    passo = np.max(np.abs(MAPA @ ctrl_novo - MAPA @ ctrl))
    f0 = mod['f0']
    obj_avl = cdff_cruzeiro(ctrl_novo, it_aft)
    print(f'  torcoes de controle: {np.round(ctrl_novo, 3)}')
    print(f'  folga alvo = {folga_alvo:.2f} (max alcancavel '
          f'{folga_max:.2f}), CDff (AVL) = {obj_avl:.6f}')
    print(f'  maior mudanca de torcao vs iteracao anterior: {passo:.3f} graus')
    historico.append({'it_fwd': it_fwd, 'it_aft': it_aft,
                      'CDff_avl': obj_avl, 'passo': passo})
    ctrl = ctrl_novo
    if passo < TOL:
        print('  convergiu.')
        break

if os.path.exists(VARIANTE):
    os.remove(VARIANTE)

torcoes = torcoes_de(ctrl)
print('\n===== OTIMO CONVERGIDO =====')
print('Torcoes nos pontos de controle: '
      + ', '.join(f'eta {e} = {v:+.3f}' for e, v in zip(ETAS_CTRL, ctrl)))
print('Torcoes por estacao [graus]:')
for eta in gera.ETAS_ASA:
    print(f'  eta = {eta}: {torcoes[eta]:+.3f}')
print(f'CDff sem torcao : {mod["f0"]:.6f}')
print(f'CDff otimizado (AVL) : {historico[-1]["CDff_avl"]:.6f}')
print(f'it de cruzeiro: fwd = {it_fwd:.3f}, aft = {it_aft:.3f}')

resultado = {'controle': {str(e): round(float(v), 3)
                          for e, v in zip(ETAS_CTRL, ctrl)},
             'torcoes': {str(eta): round(torcoes[eta], 3)
                         for eta in gera.ETAS_ASA},
             'CDff_sem_torcao': round(float(mod['f0']), 6),
             'CDff_otimizado_avl': round(float(historico[-1]['CDff_avl']), 6),
             'it_cruzeiro': {'fwd': round(float(it_fwd), 3),
                             'aft': round(float(it_aft), 3)},
             'iteracoes_externas': len(historico)}
with open('avl/saidas/torcao.json', 'w', encoding='ascii') as f:
    json.dump(resultado, f, indent=2)
print('\ngravado: avl/saidas/torcao.json')
print('agora rode: lab04_gera_avl.py, lab04_trim.py, lab04_secao_critica.py,')
print('lab04_polares.py, lab04_derivadas.py e as figuras')
