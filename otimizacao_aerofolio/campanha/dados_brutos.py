'''
Junta os resultados do Lab 03 num CSV so, para o pessoal do relatorio.

Uma linha por perfil (partida, otimizado, otimizado sem restricao de clmax,
otimizado com cusp liberado), em cada estacao. Todas as colunas em unidades
do SI ou adimensionais, ponto como separador decimal, sem formatacao.

Rodar com:  python dados_brutos.py
Saida:      ENTREGA/dados_lab03.csv
'''

import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import descritores as dsc
import otimiza_secao as ot
from analisa_otimos import carrega

DEST = os.path.normpath(os.path.join(AQUI, '..', 'ENTREGA'))
EST = ('raiz', 'meio', 'ponta')

# corda e Reynolds no plano normal ao enflechamento
GEOM = {'raiz': (0.101, 9.130, 7.729, 4.343e7),
        'meio': (0.398, 6.899, 5.841, 3.282e7),
        'ponta': (0.900, 3.125, 2.646, 1.487e7)}

# cl_max medido no XFoil (Re de decolagem, M = 0,266) e o vies aferido
CLMAX_XFOIL = {('raiz', ''): 1.9591, ('meio', ''): 2.0357,
               ('ponta', ''): 2.1420, ('meio', '_sem_bluntez'): 1.0424,
               ('raiz', '_cusp'): 2.2895, ('meio', '_cusp'): 2.1581,
               ('ponta', '_cusp'): 2.0307}
VIES_XFOIL = 0.24

VARIANTES = [('', 'otimizado'),
             ('_sem_bluntez', 'otimizado_sem_restricao_clmax'),
             ('_cusp', 'otimizado_cusp_liberado')]

COLS = ['estacao', 'variante', 'eta', 'corda_m', 'corda_normal_m',
        'Mach_normal', 'Reynolds_normal',
        'cl_projeto', 'tc_requerido',
        'alpha_deg', 'cl', 'cd', 'cm',
        'tc_max', 'x_tc_max', 'camber_max', 'x_camber_max',
        't_01', 'bluntez', 'bluntez_limiar',
        'clmax_xfoil', 'clmax_corrigido',
        'Al1', 'Al2', 'Al3', 'Al4', 'Au1', 'Au2', 'Au3', 'Au4',
        'n_avaliacoes']


def linha(estacao, variante_suf, variante_nome, xx, cd, cl, cm, n_aval):
    Al, Au, alpha = ot.desmonta(xx)
    g = dsc.descritores(Au, Al)
    eta, c, cn, re = GEOM[estacao]
    est = ot.ESTACOES[estacao]
    t01 = ot.t01_de(Al, Au)
    bl = ot.bluntez(t01, g['t_max'])
    clm = CLMAX_XFOIL.get((estacao, variante_suf), np.nan)

    d = {
        'estacao': estacao, 'variante': variante_nome,
        'eta': eta, 'corda_m': c, 'corda_normal_m': cn,
        'Mach_normal': ot.MACH_N, 'Reynolds_normal': re,
        'cl_projeto': est['cl_ref'], 'tc_requerido': est['tc_ref'],
        'alpha_deg': np.degrees(alpha), 'cl': cl, 'cd': cd, 'cm': cm,
        'tc_max': g['t_max'], 'x_tc_max': g['x_tmax'],
        'camber_max': g['c_max'], 'x_camber_max': g['x_cmax'],
        't_01': t01, 'bluntez': bl,
        'bluntez_limiar': ot.limiar_bluntez(estacao),
        'clmax_xfoil': clm,
        'clmax_corrigido': clm - VIES_XFOIL if np.isfinite(clm) else np.nan,
        'n_avaliacoes': n_aval,
    }
    for i in range(4):
        d[f'Al{i+1}'] = Al[i]
        d[f'Au{i+1}'] = Au[i]
    return d


def main():
    from analisa_otimos import roda_euler

    os.makedirs(DEST, exist_ok=True)
    linhas = []

    for e in EST:
        base = carrega(f'otim_{e}')
        if base is None:
            print(f'{e}: sem resultado')
            continue
        cfl = ot.ESTACOES[e].get('cfl', 0.20)

        # ponto de partida
        Al, Au, a = ot.desmonta(base['xx_ini'])
        r = roda_euler(Al, Au, a, 1.0, cfl=cfl)
        linhas.append(linha(e, '_partida', 'partida_NACA1411_reescalado',
                            base['xx_ini'], r['CD'], r['CL'], r['CM'], 0))
        print(f'  {e}/partida', flush=True)

        for suf, nome in VARIANTES:
            d = carrega(f'otim_{e}{suf}')
            if d is None:
                continue
            i = d['i_otimo']
            Al, Au, a = ot.desmonta(d['xx_otimo'])
            r = roda_euler(Al, Au, a, 1.0, cfl=cfl)
            linhas.append(linha(e, suf, nome, d['xx_otimo'],
                                d['hist']['CD'][i], d['hist']['CL'][i],
                                r['CM'], d['n_aval']))
            print(f'  {e}/{nome}', flush=True)

    caminho = os.path.join(DEST, 'dados_lab03.csv')
    with open(caminho, 'w') as fid:
        fid.write(','.join(COLS) + '\n')
        for l in linhas:
            fid.write(','.join(
                '' if l.get(c) is None
                else (f'{l[c]:.8g}' if isinstance(l[c], (int, float,
                                                         np.floating))
                      else str(l[c]))
                for c in COLS) + '\n')

    print(f'\n{len(linhas)} linhas em {caminho}')
    print('\ncolunas:', ', '.join(COLS))


if __name__ == '__main__':
    main()
