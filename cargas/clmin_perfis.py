'''
cl minimo (estol negativo) dos perfis do Lab 03, para o CLmin do diagrama V-n.

Mesmo procedimento do cl_max usado no Lab 04 (avl/estol.py, que le
otimizacao_aerofolio/campanha/resultados/revisao_plano_clmax.csv, plano
"normal"): perfil otimo de cada estacao na secao normal, XFoil com Ncrit = 9,
Reynolds de decolagem de cada estacao e Mach 0,266, e o vies aferido do XFoil
(0,24, afere_xfoil.py) descontado do modulo. A unica diferenca e a varredura
de alfa, de 0 a -30 graus.

Saida: resultados/01_envelope/clmin_perfis.csv
Rodar com:  python clmin_perfis.py   (alguns minutos)
'''

import csv
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
CAMPANHA = os.path.join(AQUI, '..', 'otimizacao_aerofolio', 'campanha')
sys.path.insert(0, os.path.abspath(CAMPANHA))

import xfoil_runner as xr
import otimiza_secao as ot
from analisa_otimos import carrega
from doe_clmax_xfoil import MACH_DECOLAGEM
from revisao_plano_do_clmax import RE_CORRENTE, VIES_XFOIL

ALPHA_SEQ = (0.0, -30.0, -0.5)
SAIDA = os.path.join(AQUI, 'resultados', '01_envelope', 'clmin_perfis.csv')


def clmin_estacao(estacao):
    d = carrega(f'otim_{estacao}')
    Al, Au, _ = ot.desmonta(d['xx_otimo'])
    r = xr.clmax_cst(Au, Al, Re=RE_CORRENTE[estacao], Mach=MACH_DECOLAGEM,
                     alpha_seq=ALPHA_SEQ, timeout=900)
    p = r['polar']
    i = int(np.argmin(p['CL']))
    clmin = float(p['CL'][i])
    # estol capturado: depois do minimo o cl volta a subir (em modulo, cai 2%)
    capturado = bool(np.any(p['CL'][i + 1:] > 0.98*clmin))
    return {'estacao': estacao, 'eta': ot.ESTACOES[estacao]['eta'],
            'Re': RE_CORRENTE[estacao], 'Mach': MACH_DECOLAGEM,
            'clmin_xfoil': clmin, 'alfa_clmin': float(p['alpha'][i]),
            'clmin_corrigido': clmin + VIES_XFOIL, 'estol_capturado': capturado,
            'n_pontos': len(p['CL'])}


def main():
    linhas = []
    for e in ('raiz', 'meio', 'ponta'):
        l = clmin_estacao(e)
        print(f'{e:6s} eta={l["eta"]:.3f}  clmin XFoil = {l["clmin_xfoil"]:.4f} '
              f'(alfa {l["alfa_clmin"]:.1f})  corrigido = {l["clmin_corrigido"]:.4f}  '
              f'estol capturado: {l["estol_capturado"]}', flush=True)
        linhas.append(l)
    os.makedirs(os.path.dirname(SAIDA), exist_ok=True)
    with open(SAIDA, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        w.writerows(linhas)
    print('->', SAIDA)


if __name__ == '__main__':
    main()
