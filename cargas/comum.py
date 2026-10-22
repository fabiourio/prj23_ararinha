'''
Base comum da Atividade 06 (cargas simetricas na asa).

- caminhos e acesso aos modulos do Lab 04 (avl/) e ao designTool;
- condicoes de peso no inicio e no fim do cruzeiro (fracoes de missao do
  designTool, 100% de carga paga);
- pesos que a asa carrega: estrutura, combustivel, motores e trem principal;
- geometria da asa do modelo do AVL (corda e bordo de ataque lineares);
- parametros do modelo do AVL usados aqui (arquivo, it, CLmax).

Caso critico adotado: CG dianteiro do designTool (fwd.avl, Xref = xcg_fwd)
e V_A ao nivel do mar. O V-n e construido na altitude de cruzeiro
(envelope.py).
'''

import csv
import functools
import json
import os
import re
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, '..'))
AVL_DIR = os.path.join(RAIZ, 'avl')
RESULTADOS = os.path.join(AQUI, 'resultados')
for _p in (RAIZ, AVL_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import estilo  # noqa: E402,F401  (estilo das figuras do Lab 04)
from aeronave import analisa_mod  # noqa: E402
from designTool.auxiliary import atmosphere  # noqa: E402
from designTool.constants import gravity  # noqa: E402

ARQ_AVL = 'fwd.avl'            # CG dianteiro (caso critico)
N_LIM_NEG = -1.0               # 14 CFR 25.337(c)

# Posicoes na corda (fracao da corda local, a partir do bordo de ataque)
X_ESTRUTURA_C = 0.40           # CG da estrutura da asa (designTool: xm_w + 0,4 cm_w)


def _dir(nome):
    d = os.path.join(RESULTADOS, nome)
    os.makedirs(d, exist_ok=True)
    return d


DIR_ENVELOPE = _dir('01_envelope')
DIR_AVL = _dir('02_avl')
DIR_ESFORCOS = _dir('03_esforcos')


def _linha_csv(caminho, **filtro):
    with open(caminho, encoding='utf-8') as f:
        for l in csv.DictReader(f):
            if all(l[k] == v for k, v in filtro.items()):
                return l
    raise KeyError(f'{filtro} nao encontrado em {caminho}')


ROTEIRO = os.path.join(AVL_DIR, 'resultados', 'roteiro')
IT_FWD = float(_linha_csv(os.path.join(ROTEIRO, 'tab_incidencia.csv'), cg='fwd')['it_graus'])
CLMAX = float(_linha_csv(os.path.join(ROTEIRO, 'tab_estol.csv'), cg='fwd', trimagem='com')['CLmax'])

ATM0 = atmosphere(0.0)
RHO0, A0 = ATM0['density'], ATM0['speed_of_sound']


def xref_arquivo(arquivo=ARQ_AVL):
    with open(os.path.join(AVL_DIR, arquivo), encoding='utf-8') as f:
        txt = f.read()
    return float(re.search(r'^\s*([-\d.]+)\s+[-\d.]+\s+[-\d.]+\s+#\s*Xref', txt, re.M).group(1))


@functools.lru_cache(maxsize=None)
def designtool():
    return analisa_mod()


@functools.lru_cache(maxsize=None)
def aeronave():
    '''Dados do designTool que a Atividade 06 usa (pesos em N, cotas em m).'''
    ap = designtool()
    inp, geo, tm = ap['inputs'], ap['geometry'], ap['thrust_matching']
    ew, mf, bal = ap['empty_weight'], ap['fuel_weight']['Mf_hist'], ap['balance']

    W0 = tm['W0']
    W_ini = W0*mf['start']*mf['taxi']*mf['takeoff']*mf['climb']
    W_fim = W_ini*mf['cruise']
    queimado = {'inicio': W0 - W_ini, 'fim': W0 - W_fim}

    h = inp['altitude_cruise']
    atm = atmosphere(h)

    xcg = bal['xcg_fwd']
    if abs(xcg - xref_arquivo()) > 1e-3:
        raise RuntimeError(f'Xref de {ARQ_AVL} difere do xcg_fwd do designTool ({xcg:.4f})')

    n_mlg = 2   # dois trens principais (um em cada semiasa)
    return {
        'W0': W0, 'W_empty': tm['W_empty'], 'W_fuel_missao': tm['W_fuel'],
        'W_payload': inp['W_payload'], 'W_crew': inp['W_crew'],
        'condicoes': {
            'inicio': {'nome': 'Início do cruzeiro', 'W': W_ini,
                       'W_comb': tm['W_fuel'] - queimado['inicio']},
            'fim': {'nome': 'Fim do cruzeiro', 'W': W_fim,
                    'W_comb': tm['W_fuel'] - queimado['fim']},
        },
        'Mf': dict(mf),
        'h_cruzeiro': h, 'M_cruzeiro': inp['Mach_cruise'],
        'rho_cruzeiro': atm['density'], 'a_cruzeiro': atm['speed_of_sound'],
        'S': inp['S_w'], 'c_ref': geo['cm_w'], 'b': geo['b_w'],
        'xcg': xcg,
        # pesos que a asa carrega
        'W_asa': ew['W_w'],
        'W_motor_seco': inp['engine']['weight'],
        'W_motor': ew['W_eng']/inp['n_engines'],          # instalado (fator 1,3 do designTool)
        'W_tdp': ew['W_mlg']/n_mlg,
        'y_motor': inp['y_n'], 'x_motor': ew['xcg_eng'],
        'y_tdp': inp['y_mlg'], 'x_tdp': inp['x_mlg'],
        'tanque': {'eta0': inp['b_tank_b_w_start'], 'eta1': inp['b_tank_b_w_end'],
                   'x0_c': inp['x_tank_c_w'], 'dc': inp['c_tank_c_w']},
        'W_comb_max': bal['W_maxfuel'],
        # geometria da asa (trapezoidal, origem no nariz, x para tras)
        'asa': {'xr': inp['xr_w'], 'cr': geo['cr_w'], 'xt': geo['xt_w'],
                'ct': geo['ct_w'], 'yt': geo['yt_w'],
                'diedro': inp['dihedral_w'], 'R_fus': inp['D_f']/2},
    }


def corda(y, asa):
    return asa['cr'] + (asa['ct'] - asa['cr'])*np.asarray(y)/asa['yt']


def x_ba(y, asa):
    '''x do bordo de ataque.'''
    return asa['xr'] + (asa['xt'] - asa['xr'])*np.asarray(y)/asa['yt']


def x_corda(y, frac, asa):
    '''x do ponto a `frac` da corda local.'''
    return x_ba(y, asa) + frac*corda(y, asa)


def e_eixo_elastico(av):
    '''
    Eixo elastico (fracao da corda): centro da caixa de torcao, entre as
    longarinas a 12% e 53% da corda (limites do tanque do designTool).
    '''
    t = av['tanque']
    return t['x0_c'] + t['dc']/2


def salva_json(obj, caminho):
    def conv(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.generic):
            return o.item()
        raise TypeError(type(o))
    with open(caminho, 'w', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, default=conv)


def le_json(caminho):
    with open(caminho, encoding='utf-8') as f:
        return json.load(f)


def salva_csv(linhas, caminho):
    with open(caminho, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        for l in linhas:
            w.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in l.items()})


def virgula(x, casas):
    '''Numero com virgula decimal, para os textos das figuras.'''
    return f'{x:.{casas}f}'.replace('.', ',')


GRAVIDADE = gravity
