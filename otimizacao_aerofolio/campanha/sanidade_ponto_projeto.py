'''
O ponto de projeto da asa faz sentido? -- checagem de sanidade.

A otimizacao entregou perfis bons: forma supercritica, plato supersonico no
lugar do choque forte, e c_d que bate com Korn dentro de 34%. Mas o NIVEL do
arrasto de onda esta alto demais para um aviao de transporte.

Um perfil supercritico bem projetado tem arrasto de onda quase NULO no ponto
de projeto -- e o proposito do projeto supercritico. Aqui a MAC otimizada da
0,0071, e integrando na envergadura o arrasto de onda sai em 17% do arrasto
total da aeronave, contra 1 a 3% da pratica de projeto.

Este script mostra que a causa nao e o perfil: e a ESPESSURA DA ASA. Para
cada espessura de raiz, calcula a que distancia da divergencia de arrasto
cada estacao opera. Uma asa bem dimensionada voa NO limite ou logo abaixo
dele, nunca alem.

Rodar com:  python sanidade_ponto_projeto.py
'''

import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import otimiza_secao as ot

K_KORN = 0.95
COS50 = 0.84659          # cos do enflechamento a meia corda
TCT = 0.08               # espessura da ponta (fixa no Lab 02)
MACH_N = ot.MACH_N

# espessura relativa de raiz: a nossa e as tipicas de transporte comercial
TCR_CANDIDATOS = [0.19623, 0.18, 0.16, 0.15, 0.14, 0.12]


def m_dd_secao(tc_n, cl_n):
    '''Korn no plano normal: M_dd = k - (t/c) - cl/10.'''
    return K_KORN - tc_n - cl_n / 10.0


def cd_onda(tc_n, cl_n, mach=MACH_N):
    m_crit = m_dd_secao(tc_n, cl_n) - (0.1 / 80) ** (1 / 3)
    return 20 * max(0.0, mach - m_crit) ** 4


def main():
    print(f'M_n = {MACH_N}   (M = 0,85 com enflechamento de 32,16 graus a c/2)')
    print(f'espessura de ponta fixa em {TCT}\n')
    print('Margem = M_dd - M_n. Positiva significa voar ABAIXO da divergencia')
    print('de arrasto, que e onde uma asa bem dimensionada opera.\n')

    print(f'{"tcr":>7s} | ' + ' | '.join(
        f'{n:>22s}' for n in ('raiz (eta=0,10)', 'meio (MAC)', 'ponta (eta=0,90)'))
        + ' | CDwave secoes')
    print(f'{"":>7s} | ' + ' | '.join(
        f'{"margem":>10s} {"c_d,onda":>10s}' for _ in range(3)) + ' |')
    print('-' * 100)

    b, cr, taper, S_w = 60.1518, 9.8898, 0.24, 368.833
    eta_junc = 0.1011

    for tcr in TCR_CANDIDATOS:
        col, cds_est, etas_est = [], [], []
        for nome, est in ot.ESTACOES.items():
            eta = est['eta']
            cl_n = est['cl_ref']
            # espessura streamwise interpolada linearmente, depois normal
            tc_s = tcr + (TCT - tcr) * eta
            tc_n = tc_s / COS50
            margem = m_dd_secao(tc_n, cl_n) - MACH_N
            cd = cd_onda(tc_n, cl_n)
            col.append(f'{margem:+10.4f} {cd:10.6f}')
            cds_est.append(max(cd, 1e-9))
            etas_est.append(eta)

        # integra na envergadura, como em compara_korn
        eg = np.linspace(eta_junc, 1.0, 400)
        cd_g = np.exp(np.interp(eg, etas_est, np.log(cds_est)))
        c_g = cr * (1 - (1 - taper) * eg)
        CDw = (2 * np.trapezoid(cd_g * c_g, eg) * (b / 2)
               * COS50 ** 3 / S_w)

        marca = '  <-- a nossa' if abs(tcr - 0.19623) < 1e-4 else ''
        print(f'{tcr:7.3f} | ' + ' | '.join(col) + f' | {CDw:12.6f}{marca}')

    print(f'\n{"="*100}')
    print('LEITURA')
    print('=' * 100)
    print('Com tcr = 0,196 (a nossa), a RAIZ opera 0,048 alem da divergencia')
    print('de arrasto e a MAC 0,020 alem. So a ponta, com 8% de espessura,')
    print('tem margem positiva (+0,045) -- e ela e a estacao de menor corda,')
    print('logo a que menos contribui para o arrasto da asa.')
    print('\nNao existe perfil, por melhor que seja, que elimine arrasto de')
    print('onda com margem negativa: ele e imposto pela espessura e pelo')
    print('Mach, nao pela forma. A otimizacao so escolhe COMO gastar o que ja')
    print('esta imposto.')
    print('\nCom tcr entre 0,14 e 0,15 -- a faixa da industria para transporte')
    print('a M = 0,85 -- a raiz passa a ter margem positiva e o arrasto de')
    print('onda da asa cai 72 a 80%, porque a dependencia e de 4a potencia.')
    print('\nOu seja: os c_d que medimos estao CORRETOS, e sao a evidencia de')
    print('que a asa da aeronave B e grossa demais para o Mach de cruzeiro.')
    print('A otimizacao de perfil extraiu o que havia para extrair (-79% na')
    print('MAC); o que sobra e problema de dimensionamento da asa, e so o')
    print('Lab 02 pode resolver.')


if __name__ == '__main__':
    main()
