'''
Lab 04 -- decomposicao da estabilidade direcional (Cn_beta) do aft.avl.

A Tabela 9 (derivadas_estabilidade.py) da Cn_beta < 0 em M 0,85, isto e,
aeronave direcionalmente instavel no AVL. Este script mostra de onde vem
esse valor: recorta o aft.avl final em blocos SURFACE/BODY (como o
ponto_neutro.py) e roda o AVL com componentes removidos:

  completo            asa + winglet + EH + EV + fuselagem + nacele
  sem fuselagem       tira o BODY
  sem nacele          tira o anel da nacele
  sem EV              tira a empenagem vertical
  so fuselagem        so o BODY
  sem fus. e nacele   asa + winglet + EH + EV

Condicao: alfa = 3 graus (proximo do alfa de cruzeiro, 3,9 graus), beta = 0,
controles e it em zero, M 0,2 e M 0,85. Le CY_beta, Cl_beta e Cn_beta do
"st" (eixos de estabilidade, convencao do AVL, sem as inversoes de sinal do
programa de MVO; Cn_beta e Cl_beta nao mudam de sinal nessa conversao).

Saidas em resultados/derivadas/: decomposicao_cnb.csv e decomposicao_cnb.png.
Arquivos temporarios: resultados/_tmp/cnb_*.avl.

Rodar de dentro de avl/:   python estabilidade_direcional.py [--avl aft.avl]
'''

import argparse
import csv
import os
import re
import sys

import numpy as np

from estilo import plt, AZUL, LARANJA, TINTA2
from avl_run import roda
from avl_saida import derivadas
from gera_avl import AQUI

MACHS = (0.2, 0.85)
ALFA = 3.0
SAIDA = os.path.join(AQUI, 'resultados', 'derivadas')

TODOS = ['Wing', 'Winglet', 'Horizontal tail', 'Vertical tail', 'Fuselage', 'Nacelle']
CONFIGS = [  # (slug, rotulo, blocos incluidos)
    ('completo', 'completo', TODOS),
    ('sem_fus', 'sem fuselagem', [b for b in TODOS if b != 'Fuselage']),
    ('sem_nac', 'sem nacele', [b for b in TODOS if b != 'Nacelle']),
    ('sem_ev', 'sem EV', [b for b in TODOS if b != 'Vertical tail']),
    ('so_fus', 'só fuselagem', ['Fuselage']),
    ('sem_fus_nac', 'sem fuselagem e nacele',
     [b for b in TODOS if b not in ('Fuselage', 'Nacelle')]),
]


def blocos(caminho):
    '''Divide o .avl em {'cab': cabecalho, nome: bloco SURFACE/BODY}.'''
    with open(os.path.join(AQUI, caminho), encoding='utf-8') as f:
        txt = f.read()
    out = {}
    for bl in re.split(r'(?=#-{20,}\n)', txt):
        m = re.search(r'(SURFACE|BODY)\n(.+)\n', bl)
        out[m.group(2).strip() if m else 'cab'] = bl
    faltam = [b for b in TODOS if b not in out]
    if faltam:
        raise RuntimeError(f'blocos ausentes em {caminho}: {faltam}')
    return out


def roda_config(arquivo, mach):
    s = roda([f'load {arquivo}', 'oper', 'm', f'mn {mach}', '',
              f'a a {ALFA}', 'x', 'st', '', '', 'quit'])
    d = derivadas(s, 'st')
    return {k: d[k] for k in ('CYb', 'Clb', 'Cnb')}


def figura(linhas, destino):
    rot = [c[1] for c in CONFIGS]
    y = np.arange(len(rot))
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    for k, (M, cor) in enumerate(zip(MACHS, (LARANJA, AZUL))):
        v = [next(l['Cnb'] for l in linhas if l['config'] == r and l['M'] == M) for r in rot]
        ax.barh(y + (k - 0.5)*0.36, v, height=0.34, color=cor,
                label=f'M {M:g}'.replace('.', ','))
        for yi, vi in zip(y + (k - 0.5)*0.36, v):
            ax.text(vi + (0.006 if vi >= 0 else -0.006), yi, f'{vi:+.3f}'.replace('.', ','),
                    va='center', ha='left' if vi >= 0 else 'right', fontsize=7, color=TINTA2)
    ax.axvline(0.0, color=TINTA2, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(rot)
    ax.invert_yaxis()
    lim = max(abs(l['Cnb']) for l in linhas)*1.35
    ax.set_xlim(-lim, lim)
    ax.set_xlabel('Cnβ [1/rad]  (> 0: estável)')
    alfa_txt = f'{ALFA:g}'.replace('.', ',')
    ax.set_title(f'Decomposição de Cnβ por componente (aft.avl, α = {alfa_txt}°)', fontsize=9)
    ax.legend(fontsize=8, loc='upper left')
    fig.tight_layout()
    fig.savefig(destino)
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description='Decomposicao de Cn_beta por componente')
    p.add_argument('--avl', default='aft.avl')
    a = p.parse_args(argv)

    B = blocos(a.avl)
    os.makedirs(os.path.join(AQUI, 'resultados', '_tmp'), exist_ok=True)
    os.makedirs(SAIDA, exist_ok=True)
    linhas = []
    for slug, rotulo, incl in CONFIGS:
        arq = f'resultados/_tmp/cnb_{slug}.avl'
        with open(os.path.join(AQUI, arq), 'w', encoding='utf-8') as f:
            f.write(''.join(B[k] for k in ['cab'] + incl))
        for M in MACHS:
            r = roda_config(arq, M)
            linhas.append({'config': rotulo, 'M': M, 'alfa_graus': ALFA, **r,
                           'blocos': ' + '.join(incl)})
            print(f'  {rotulo:24s} M {M:4.2f}   CYb {r["CYb"]:+.4f}   Clb {r["Clb"]:+.4f}'
                  f'   Cnb {r["Cnb"]:+.4f}')

    with open(os.path.join(SAIDA, 'decomposicao_cnb.csv'), 'w', newline='',
              encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)
    figura(linhas, os.path.join(SAIDA, 'decomposicao_cnb.png'))
    print('salvo resultados/derivadas/decomposicao_cnb.csv e .png')


if __name__ == '__main__':
    main(sys.argv[1:])
