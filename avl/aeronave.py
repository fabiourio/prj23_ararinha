'''
Ponte entre o designTool e o Lab 04.

aeronave(Lc_h) roda o designTool com o braco da EH pedido e devolve um
dicionario com tudo o que os scripts do AVL usam: ponto de projeto,
referencias, CGs, dε/dα e a geometria de asa, EH, EV, nacele e fuselagem.
As alturas (z) sao as REAIS do designTool; o deslocamento em Z exigido
no AVL e aplicado so no gera_avl.py.

O peso no ponto de projeto segue a definicao do Lab 03 (peso medio de
cruzeiro, 229.669,3 kgf com Lc_h = 4,6). Para outros Lc_h mantem-se a mesma
fracao de combustivel, com 100% de carga paga, e o peso e recalculado com o
W_empty e o W_fuel do designTool.

aeronave_mod(overrides) faz o mesmo para inputs arbitrarios do designTool
(estudos de sensibilidade, ex.: sensibilidade_estabilidade.py), sem mudar
o comportamento de aeronave(Lc_h).
'''

import copy
import functools
import os
import sys

RAIZ_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if RAIZ_REPO not in sys.path:
    sys.path.insert(0, RAIZ_REPO)

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.geometry import geometry
from designTool.aerodynamics import aerodynamics
from designTool.auxiliary import atmosphere
from designTool.constants import gravity

LC_H_BASE = 4.6
W_PROJETO_BASE_KGF = 229669.3
FOLGA_FUSELAGEM = 0.5      # [m] bordo de fuga da raiz da EH ate o fim da fuselagem


def _analisa(Lc_h):
    ap = standard_airplane('my_airplane')
    ap['inputs']['Lc_h'] = float(Lc_h)
    analyze(ap, print_log=False, plot=False)
    return ap


@functools.lru_cache(maxsize=None)
def fracao_combustivel():
    '''Fracao de combustivel que reproduz o peso de projeto do Lab 03.'''
    ap = _analisa(LC_H_BASE)
    tm, inp = ap['thrust_matching'], ap['inputs']
    W = W_PROJETO_BASE_KGF*gravity
    return (W - tm['W_empty'] - inp['W_payload'] - inp['W_crew'])/tm['W_fuel']


def _dicionario(ap):
    '''Dicionario da aeronave a partir de um designTool ja analisado.'''
    inp, geo = ap['inputs'], ap['geometry']
    tm, bal = ap['thrust_matching'], ap['balance']

    h, M = inp['altitude_cruise'], inp['Mach_cruise']
    atm = atmosphere(h)
    rho, a_inf = atm['density'], atm['speed_of_sound']
    V = M*a_inf
    q = 0.5*rho*V**2
    S_w = inp['S_w']

    W = (tm['W_empty'] + inp['W_payload'] + inp['W_crew']
         + fracao_combustivel()*tm['W_fuel'])
    CL = W/(q*S_w)
    _, _, drag = aerodynamics(ap, Mach=M, altitude=h, CL=CL)

    return {
        'Lc_h': float(inp['Lc_h']),
        'W0': tm['W0'], 'W': W, 'W_kgf': W/gravity, 'W0_kgf': tm['W0']/gravity,
        'fuel_frac': fracao_combustivel(),
        'M': M, 'h': h, 'rho': rho, 'a': a_inf, 'V': V, 'CL': CL,
        'CD0': drag['CD0'],
        'Sref': S_w, 'Cref': geo['cm_w'], 'Bref': geo['b_w'],
        'xcg_fwd': bal['xcg_fwd'], 'xcg_aft': bal['xcg_aft'],
        'xnp_dt': bal['xnp'], 'deda': bal['deda'],
        'asa': {'xr': inp['xr_w'], 'zr': inp['zr_w'], 'cr': geo['cr_w'],
                'ct': geo['ct_w'], 'xt': geo['xt_w'], 'yt': geo['yt_w'],
                'zt': geo['zt_w']},
        'EH': {'S': geo['S_h'], 'xr': geo['xr_h'], 'zr': inp['zr_h'],
               'cr': geo['cr_h'], 'ct': geo['ct_h'], 'xt': geo['xt_h'],
               'yt': geo['yt_h'], 'zt': geo['zt_h'], 'tc': inp['tcr_h']},
        'EV': {'xr': geo['xr_v'], 'zr': inp['zr_v'], 'cr': geo['cr_v'],
               'ct': geo['ct_v'], 'xt': geo['xt_v'], 'zt': geo['zt_v']},
        'nacele': {'x': inp['x_n'], 'y': inp['y_n'], 'z': inp['z_n'],
                   'L': inp['L_n'], 'D': inp['D_n']},
        'fuselagem': {'L': inp['L_f'], 'D': inp['D_f']},
    }


@functools.lru_cache(maxsize=None)
def _aeronave_cache(Lc_h):
    '''Dicionario da aeronave para um Lc_h (cru, compartilhado pelo cache).'''
    return _dicionario(_analisa(Lc_h))


def aeronave(Lc_h=LC_H_BASE):
    '''Dicionario da aeronave para um Lc_h. Copia independente do cache: pode
    ser modificado livremente sem afetar chamadas futuras.'''
    return copy.deepcopy(_aeronave_cache(float(Lc_h)))


def analisa_mod(overrides=None):
    '''
    Roda o designTool com standard_airplane('my_airplane') e os inputs
    trocados por `overrides` (dict {nome_do_input: valor}); devolve o
    dicionario completo do designTool. Sem overrides e a aeronave entregue.
    '''
    ap = standard_airplane('my_airplane')
    ap['inputs'].update(copy.deepcopy(overrides or {}))
    analyze(ap, print_log=False, plot=False)
    return ap


def aeronave_mod(overrides=None):
    '''
    Como aeronave(), mas para inputs arbitrarios do designTool (estudos de
    sensibilidade: Cht, xr_w, x_n, ...). Mesma definicao de peso de projeto
    (fracao de combustivel do Lab 03, 100% de carga paga). Sem cache: cada
    chamada roda o designTool e devolve um dicionario novo.
    '''
    return _dicionario(analisa_mod(overrides))


def _fim_eh(Lc_h, overrides=None):
    ap = standard_airplane('my_airplane')
    ap['inputs'].update(copy.deepcopy(overrides or {}))
    ap['inputs']['Lc_h'] = float(Lc_h)
    geometry(ap)
    return ap['geometry']['xr_h'] + ap['geometry']['cr_h'], ap['inputs']['L_f']


def lc_h_maximo(folga=FOLGA_FUSELAGEM, lo=4.0, hi=6.0, tol=1e-5, overrides=None):
    '''Maior Lc_h com o bordo de fuga da raiz da EH a `folga` do fim da
    fuselagem. `overrides`: inputs do designTool trocados (ex.: xr_w, Cht).'''
    while hi - lo > tol:
        mid = 0.5*(lo + hi)
        fim, L_f = _fim_eh(mid, overrides)
        if fim <= L_f - folga:
            lo = mid
        else:
            hi = mid
    return lo
