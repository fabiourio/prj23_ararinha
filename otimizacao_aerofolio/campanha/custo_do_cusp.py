'''
O ganho do cusp sobrevive no nivel da aeronave? -- arrasto de compensacao.

Soltar o batente do ultimo coeficiente do intradorso rendeu 29% de arrasto na
MAC. Mas o momento de arfagem quase dobrou: c_m de -0,097 para -0,165.

A nossa funcao objetivo e 2D e NAO ENXERGA isso. Momento maior significa mais
download na empenagem horizontal, que a asa precisa compensar com mais
sustentacao, que custa arrasto induzido -- na asa e na propria empenagem.

Este script estima esse custo com a geometria da aeronave B e compara com o
ganho de arrasto de onda. Se as duas grandezas forem da mesma ordem, o ganho
do cusp e ilusorio sem uma restricao de c_m.

As estimativas sao grosseiras de proposito (equilibrio rigido, sem
contribuicao da fuselagem, arrasto induzido parabolico). O objetivo e a ORDEM
DE GRANDEZA, que e o que decide se vale investigar a fundo.

Rodar com:  python custo_do_cusp.py
'''

import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import otimiza_secao as ot
from analisa_otimos import carrega
from compara_korn import FATOR_MALHA

# aeronave B (documento de projeto, secao 2) e saidas do designTool
S_W, MAC, B_W = 368.833, 6.899, 60.1518
CHT, LC_H, AR_H = 0.7, 4.0, 4.6
AR_EFF = 11.772
CL_CRUZEIRO = 0.5053
CDIND_CRUZEIRO = 0.008674      # do designTool no peso medio de cruzeiro
E_TAIL = 0.85
COS50 = 0.84659
TAPER, CR = 0.24, 9.8898
ETA_JUNC = 0.1011


def cm_e_cd(nome, sufixo):
    d = carrega(f'otim_{nome}{sufixo}')
    if d is None:
        return None
    i = d['i_otimo']
    f = FATOR_MALHA.get(nome, 1.0)
    return d['hist']['CD'][i] * f


def integra(cds_por_estacao):
    '''c_d de secao -> CD da asa, com a conversao de enflechamento.'''
    etas = np.array([ot.ESTACOES[n]['eta'] for n in ('raiz', 'meio', 'ponta')])
    cds = np.array([cds_por_estacao[n] for n in ('raiz', 'meio', 'ponta')])
    eg = np.linspace(ETA_JUNC, 1.0, 400)
    cd_g = np.exp(np.interp(eg, etas, np.log(cds)))
    c_g = CR * (1 - (1 - TAPER) * eg)
    return 2 * np.trapezoid(cd_g * c_g, eg) * (B_W / 2) * COS50 ** 3 / S_W


def trim(cm_asa):
    '''
    Equilibrio rigido: CM da asa e reagido pela empenagem.
        CM * S_w * MAC = CL_h * S_h * L_h    ->   CL_h = CM / Cht
    A asa compensa o download com sustentacao extra, e as duas pagam arrasto
    induzido.
    '''
    S_h = CHT * S_W * MAC / (LC_H * MAC)
    CL_h = cm_asa / CHT
    dCL_asa = abs(CL_h) * S_h / S_W
    CL_novo = CL_CRUZEIRO + dCL_asa
    dCDi_asa = CDIND_CRUZEIRO * ((CL_novo / CL_CRUZEIRO) ** 2 - 1)
    CDi_emp = CL_h ** 2 * (S_h / S_W) / (np.pi * AR_H * E_TAIL)
    return S_h, CL_h, dCDi_asa, CDi_emp


def main():
    print('Aeronave B: S_w = %.1f m2, MAC = %.3f m, Cht = %.2f, AR_h = %.1f\n'
          % (S_W, MAC, CHT, AR_H))

    cms = {'roteiro': {}, 'cusp': {}}
    cds = {'roteiro': {}, 'cusp': {}}
    # c_m lidos das tabelas (rodadas no nivel 1,0)
    CM = {'raiz': (-0.0513944, -0.149806),
          'meio': (-0.0969922, -0.165102),
          'ponta': (-0.108717, -0.168817)}
    for nome in ot.ESTACOES:
        cds['roteiro'][nome] = cm_e_cd(nome, '')
        cds['cusp'][nome] = cm_e_cd(nome, '_cusp')
        cms['roteiro'][nome], cms['cusp'][nome] = CM[nome]

    if any(v is None for d in cds.values() for v in d.values()):
        print('faltam resultados')
        return 1

    print(f'{"variante":10s} {"CDwave asa":>11s} {"c_m medio":>10s} '
          f'{"CL_emp":>8s} {"dCDi asa":>9s} {"CDi emp":>9s} {"CD trim":>9s}')
    linhas = {}
    for var in ('roteiro', 'cusp'):
        CDw = integra(cds[var])
        cm_med = float(np.mean(list(cms[var].values())))
        S_h, CL_h, dCDi_a, CDi_e = trim(cm_med)
        CD_trim = dCDi_a + CDi_e
        linhas[var] = (CDw, cm_med, CL_h, dCDi_a, CDi_e, CD_trim)
        print(f'{var:10s} {CDw:11.6f} {cm_med:10.4f} {CL_h:8.4f} '
              f'{dCDi_a:9.6f} {CDi_e:9.6f} {CD_trim:9.6f}')

    ganho_onda = linhas['roteiro'][0] - linhas['cusp'][0]
    custo_trim = linhas['cusp'][5] - linhas['roteiro'][5]
    liquido = ganho_onda - custo_trim

    print(f'\n{"="*68}\nBALANCO NO NIVEL DA AERONAVE\n{"="*68}')
    print(f'  ganho de arrasto de onda com o cusp : {ganho_onda:+.6f}')
    print(f'  custo de arrasto de compensacao     : {custo_trim:+.6f}')
    print(f'  saldo liquido                       : {liquido:+.6f}')
    if abs(custo_trim) > 0.5 * abs(ganho_onda):
        print('\n  O custo de compensacao e da MESMA ORDEM do ganho.')
        print('  O ganho de 29% medido em 2D NAO se transfere direto para a')
        print('  aeronave: boa parte volta como arrasto induzido de trimagem.')
        print('  A funcao objetivo 2D nao enxerga isso, entao o resultado do')
        print('  item 9 precisa vir acompanhado de uma restricao de c_m, ou')
        print('  de uma contabilidade no nivel da aeronave.')
    else:
        print('\n  O custo de compensacao e pequeno perto do ganho: o resultado')
        print('  do item 9 se sustenta no nivel da aeronave.')
    print('\n  Ressalvas: equilibrio rigido, sem contribuicao da fuselagem nem')
    print('  do enflechamento no braco, arrasto induzido parabolico e c_m')
    print('  medio das tres estacoes em vez de integrado. E ordem de')
    print('  grandeza, nao numero de projeto.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
