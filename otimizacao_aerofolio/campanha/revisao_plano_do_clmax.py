'''
REVISAO: o cl_max foi verificado no plano certo?

O perfil otimizado e a secao NORMAL ao enflechamento -- correto para o projeto
transonico, porque e esse o plano onde mora o choque. Mas o clmax_w do
designTool entra na formula de Raymer

    CLmax_clean = 0,9 * clmax_w * cos(Lambda)

e a convencao dessa formula e que clmax_w seja o cl_max do perfil na direcao
da CORRENTE, que e como perfis de asa sao especificados em desenho. O cos(Λ)
ali JA e a correcao de enflechamento.

As duas secoes sao geometricamente diferentes. Normalizando a corrente para
corda unitaria, a secao da corrente e a secao normal com as ordenadas
multiplicadas por cos(Lambda):

    x_corrente = x_normal / cos(L)  ,  y igual   ->  corda 1/cos(L)
    normalizando: x igual, y *= cos(L)

ou seja, (t/c)_corrente = (t/c)_normal * cos(L) -- mais FINA. E perfil mais
fino tende a ter cl_max MENOR.

Se verificamos o cl_max na secao normal e comparamos com um alvo que espera a
secao da corrente, superestimamos. Este script mede as duas.

Rodar com:  python revisao_plano_do_clmax.py
'''

import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import descritores as dsc
import otimiza_secao as ot
import xfoil_runner as xr
from analisa_otimos import carrega
from doe_clmax_xfoil import ALPHA_SEQ, MACH_DECOLAGEM

COS_L = 0.84659            # cos do enflechamento a meia corda (32,158 graus)
VIES_XFOIL = 0.24          # aferido em afere_xfoil.py

# Reynolds de decolagem baseado na corda DA CORRENTE de cada estacao
RE_CORRENTE = {'raiz': 6.06e7, 'meio': 4.23e7, 'ponta': 1.45e7}

ALVO = 1.80


def secao_da_corrente(Au, Al):
    '''
    Converte a secao normal para a secao da corrente, normalizada.
    Basta multiplicar as ordenadas por cos(Lambda): na parametrizacao CST,
    y e linear nos coeficientes, entao escalar os coeficientes escala y.
    '''
    return np.asarray(Au) * COS_L, np.asarray(Al) * COS_L


def main():
    print('Revisao: o cl_max foi medido na secao NORMAL ou na da CORRENTE?\n')
    print(f'cos(Lambda_c/2) = {COS_L:.5f}\n')
    print(f'{"estacao":8s} {"plano":10s} {"t/c":>7s} {"Re":>10s} '
          f'{"clmax":>7s} {"corrigido":>10s} {"vs alvo":>9s}')

    linhas = []
    for e in ot.ESTACOES:
        d = carrega(f'otim_{e}')
        if d is None:
            continue
        Al, Au, _ = ot.desmonta(d['xx_otimo'])
        Re = RE_CORRENTE[e]

        for plano, (au, al) in [('normal', (Au, Al)),
                                ('corrente', secao_da_corrente(Au, Al))]:
            g = dsc.descritores(au, al)
            r = xr.clmax_cst(au, al, Re=Re, Mach=MACH_DECOLAGEM,
                             alpha_seq=ALPHA_SEQ, timeout=900)
            corr = r['clmax'] - VIES_XFOIL
            print(f'{e:8s} {plano:10s} {g["t_max"]:7.4f} {Re:10.2e} '
                  f'{r["clmax"]:7.4f} {corr:10.4f} '
                  f'{corr - ALVO:+9.4f}   ({r["situacao"]})')
            linhas.append({'estacao': e, 'plano': plano,
                           'tc': g['t_max'], 'Re': Re,
                           'clmax_xfoil': r['clmax'], 'clmax_corrigido': corr,
                           'situacao': r['situacao']})
        print()

    cam = os.path.join(AQUI, 'resultados', 'revisao_plano_clmax.csv')
    with open(cam, 'w') as fid:
        cols = list(linhas[0])
        fid.write(','.join(cols) + '\n')
        for l in linhas:
            fid.write(','.join(
                f'{l[c]:.6g}' if isinstance(l[c], float) else str(l[c])
                for c in cols) + '\n')

    print('=' * 72)
    corr_n = {l['estacao']: l['clmax_corrigido']
              for l in linhas if l['plano'] == 'normal'}
    corr_c = {l['estacao']: l['clmax_corrigido']
              for l in linhas if l['plano'] == 'corrente'}
    print('LEITURA')
    print('=' * 72)
    for e in corr_n:
        print(f'  {e:8s}: normal {corr_n[e]:.3f}  ->  corrente '
              f'{corr_c[e]:.3f}   (diferenca {corr_c[e]-corr_n[e]:+.3f})')
    atende_n = sum(v >= ALVO for v in corr_n.values())
    atende_c = sum(v >= ALVO for v in corr_c.values())
    print(f'\n  estacoes que atendem o alvo de {ALVO}:')
    print(f'    medindo na secao normal   : {atende_n} de {len(corr_n)}')
    print(f'    medindo na secao da corrente: {atende_c} de {len(corr_c)}')
    print(f'\n  gravado em {cam}')


if __name__ == '__main__':
    main()
