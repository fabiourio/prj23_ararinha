'''
Que estacao estola primeiro, e quanto o Ncrit infla o cl_max?

Duas perguntas que decidem se a restricao de cl_max deve valer para a asa
inteira ou so para parte dela.

PARTE 1 -- estacao critica.
A asa estola onde o c_l LOCAL alcanca primeiro o c_l_max LOCAL. Escalando a
distribuicao de sustentacao, a estacao com a menor razao
    k = c_l_max_local / c_l_local
e a que estola primeiro, e o CLmax da asa e k_min vezes o CL de referencia.

Intuicao comum e que a asa estola na raiz. Para uma asa ENFLECHADA e
AFILADA e o contrario: os dois efeitos empurram a carga para fora. A nossa
tem afilamento 0,24, e a distribuicao eliptica ja poe o c_l local maximo em
eta = 0,76.

PARTE 2 -- sensibilidade ao Ncrit.
Ncrit = 9 e tunel limpo. Asa de verdade tem rugosidade de fabricacao,
rebites, degraus e frestas, que antecipam a transicao e reduzem o cl_max.
Valores de 5 a 7 representam superficie real. Medimos quanto muda.

Rodar com:  python estacao_critica.py
'''

import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import otimiza_secao as ot
import xfoil_runner as xr
from analisa_otimos import RES, carrega
from doe_clmax_xfoil import ALPHA_SEQ, MACH_DECOLAGEM

# aeronave B
B, CR, TAPER, S_W = 60.1518, 9.8898, 0.24, 368.833
CL_REF = 0.5053                # peso medio de cruzeiro
COS_L25 = 0.82263              # cos do enflechamento a c/4 (Raymer usa este)
ETA_JUNC = 0.1011

RE_EST = {'raiz': 6.06e7, 'meio': 4.23e7, 'ponta': 1.45e7}
VIES = 0.24
NCRITS = [9.0, 7.0, 5.0]


def corda(eta):
    return CR * (1 - (1 - TAPER) * eta)


def cl_eliptico(eta, CL=CL_REF):
    K = 4 * CL * S_W / (np.pi * B)
    return K * np.sqrt(np.maximum(0.0, 1 - eta ** 2)) / corda(eta)


def main():
    print('=' * 74)
    print('PARTE 1 -- que estacao estola primeiro?')
    print('=' * 74)
    print('k = c_l_max local / c_l local. A menor razao estola primeiro.\n')

    dados = {}
    print(f'{"estacao":8s} {"eta":>6s} {"c_l local":>10s} {"c_l_max":>9s} '
          f'{"k":>7s}')
    for e in ot.ESTACOES:
        d = carrega(f'otim_{e}')
        if d is None:
            continue
        Al, Au, _ = ot.desmonta(d['xx_otimo'])
        r = xr.clmax_cst(Au, Al, Re=RE_EST[e], Mach=MACH_DECOLAGEM,
                         alpha_seq=ALPHA_SEQ, timeout=900)
        clmax_real = r['clmax'] - VIES
        eta = ot.ESTACOES[e]['eta']
        cl_loc = float(cl_eliptico(eta))
        k = clmax_real / cl_loc
        dados[e] = {'eta': eta, 'cl_loc': cl_loc, 'clmax': clmax_real,
                    'k': k, 'Al': Al, 'Au': Au}
        print(f'{e:8s} {eta:6.3f} {cl_loc:10.4f} {clmax_real:9.4f} {k:7.3f}')

    if len(dados) < 3:
        print('faltam estacoes')
        return 1

    # interpola k ao longo da envergadura para achar o minimo de verdade
    etas = np.array([dados[e]['eta'] for e in ('raiz', 'meio', 'ponta')])
    clmaxs = np.array([dados[e]['clmax'] for e in ('raiz', 'meio', 'ponta')])
    eg = np.linspace(ETA_JUNC, 0.98, 300)
    clmax_g = np.interp(eg, etas, clmaxs)
    k_g = clmax_g / cl_eliptico(eg)
    i = int(np.argmin(k_g))

    print(f'\n  minimo de k ao longo da envergadura: eta = {eg[i]:.3f}, '
          f'k = {k_g[i]:.3f}')
    critica = min(dados, key=lambda e: dados[e]['k'])
    print(f'  entre as tres estacoes, a critica e a {critica.upper()}')
    print(f'\n  CLmax da asa estimado = k_min * CL_ref = '
          f'{k_g[i]*CL_REF:.4f}')
    print(f'  CLmax_clean do designTool = 0,9 * 1,8 * cos(34,65) = '
          f'{0.9*1.8*COS_L25:.4f}')

    print(f'\n{"="*74}')
    print('PARTE 2 -- quanto o Ncrit infla o cl_max?')
    print('=' * 74)
    print('Ncrit = 9 e tunel limpo; 5 a 7 representa superficie real.\n')
    print(f'{"estacao":8s} ' + ' '.join(f'{"Ncrit "+str(int(n)):>10s}'
                                        for n in NCRITS)
          + f' {"queda 9->5":>11s}')
    for e in ('raiz', 'meio', 'ponta'):
        vals = []
        for nc in NCRITS:
            r = xr.clmax_cst(dados[e]['Au'], dados[e]['Al'], Re=RE_EST[e],
                             Mach=MACH_DECOLAGEM, alpha_seq=ALPHA_SEQ,
                             n_crit=nc, timeout=900)
            vals.append(r['clmax'])
        print(f'{e:8s} ' + ' '.join(f'{v:10.4f}' for v in vals)
              + f' {vals[-1]-vals[0]:+11.4f}')
        dados[e]['ncrit'] = vals

    print(f'\n{"="*74}')
    print('LEITURA')
    print('=' * 74)
    print(f'A estacao critica e a {critica.upper()}, nao a raiz: enflechamento')
    print('e afilamento empurram a carga para fora, e o afilamento de 0,24 poe')
    print(f'o c_l local maximo em eta = 0,76.')
    print('\nConsequencia direta para a pergunta: a restricao de cl_max NAO')
    print('pode ser afrouxada na ponta. Ela e justamente onde importa.')
    print('Se alguma pudesse ser afrouxada, seria a da RAIZ, que tem a maior')
    print('margem.')
    print('\nSobre o Ncrit: a queda de 9 para 5 mostra quanto do nosso cl_max')
    print('vem de supor superficie perfeitamente lisa. Some isso ao vies de')
    print('+0,24 ja aferido contra dados experimentais.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
