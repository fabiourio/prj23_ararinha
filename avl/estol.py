'''
Metodo da secao critica.

O cl_norm de cada faixa da asa e linear em alfa no VLM (inclusive com o
profundor de compensacao, que tambem varia linearmente). Duas rodadas dao o
alfa em que cada faixa alcanca o clmax do seu perfil; a menor define o
estol. Uma terceira rodada, no alfa encontrado, confere o resultado.

clmax dos perfis do Lab 03 (XFoil corrigido, otimizacao_aerofolio/campanha/
resultados/revisao_plano_clmax.csv, plano "normal"), interpolado entre as
estacoes de referencia e constante fora delas.
'''

import numpy as np

from avl_run import caso

ETA_CLMAX = [0.1011, 0.398, 0.90]
CLMAX = [1.774, 1.7985, 1.7338]
MACH_BAIXO = 0.2
ALFAS_BASE = (8.0, 14.0)


def clmax_local(eta):
    return np.interp(eta, ETA_CLMAX, CLMAX)


def modelo_estol(arquivo, it, trim, semi, empenagem=True):
    '''alfa de estol de cada faixa da asa direita (modelo linear).'''
    a0, a1 = ALFAS_BASE
    r0, r1 = (caso(arquivo, MACH_BAIXO, alfa=a, it=it, trim=trim,
                   empenagem=empenagem, faixas=True) for a in ALFAS_BASE)
    eta = r0['faixas']['y']/semi
    c0, c1 = r0['faixas']['cl_norm'], r1['faixas']['cl_norm']
    b = (c1 - c0)/(a1 - a0)
    a = c0 - b*a0
    return {'eta': eta, 'alfa_i': (clmax_local(eta) - a)/b}


def estol(arquivo, it, trim, semi, empenagem=True):
    m = modelo_estol(arquivo, it, trim, semi, empenagem)
    i = int(np.argmin(m['alfa_i']))
    alfa = float(m['alfa_i'][i])
    r = caso(arquivo, MACH_BAIXO, alfa=alfa, it=it, trim=trim,
             empenagem=empenagem, faixas=True)
    eta = r['faixas']['y']/semi
    lim = clmax_local(eta)
    return {'alfa': alfa, 'CL': r['CL'], 'de': r['de'] if trim else 0.0,
            'eta_crit': float(m['eta'][i]), 'eta': eta,
            'cl_norm': r['faixas']['cl_norm'], 'clmax': lim,
            'excesso_max': float(np.max(r['faixas']['cl_norm'] - lim)),
            'alfa_i': m['alfa_i']}
