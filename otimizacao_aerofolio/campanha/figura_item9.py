'''
Figura do item 9 do roteiro: o que muda ao alterar a definicao do problema.

Compara, por estacao, o otimo obtido com os batentes sugeridos pelo roteiro
(Al <= -0,05) e o obtido soltando o ultimo coeficiente do intradorso, que e o
que constroi o cusp concavo do bordo de fuga.

A comparacao importa porque, com os batentes do roteiro, os tres otimos
estavam ENCOSTADOS no batente -- nao eram otimos interiores. O c_d estava
limitado por uma escolha de caixa, nao pela fisica.

Rodar com:  python figura_item9.py
'''

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import otimiza_secao as ot
import xfoil_runner as xr
from analisa_otimos import RES, carrega, roda_euler
from estilo import INK, INK2, MUTED, PAL, salvar, style_axes

VARIANTES = [('', 'batentes do roteiro', MUTED),
             ('_cusp', 'com o cusp liberado', PAL[0])]


def figura(nome):
    est = ot.ESTACOES[nome]
    dados = []
    for suf, rot, cor in VARIANTES:
        d = carrega(f'otim_{nome}{suf}')
        if d is None:
            return
        Al, Au, alpha = ot.desmonta(d['xx_otimo'])
        r = roda_euler(Al, Au, alpha, 1.0, cfl=est.get('cfl', 0.20))
        dados.append((rot, cor, Al, Au, r,
                      d['hist']['CD'][d['i_otimo']]))

    fig, axs = plt.subplots(3, 1, figsize=(9, 11),
                            gridspec_kw={'height_ratios': [1.0, 1.3, 1.3]})

    # --- geometria ---
    ax = axs[0]
    style_axes(ax)
    for rot, cor, Al, Au, r, cd in dados:
        f = xr.cst_coords(Au, Al)
        ax.plot(np.real(f['x_coord']), np.real(f['y_coord']), '-',
                color=cor, linewidth=2.0, label=f'{rot}  ($c_d$ = {cd:.5f})')
    ax.axis('equal')
    ax.set_ylabel('$y/c$', color=INK2, fontsize=10)
    leg = ax.legend(fontsize=9, frameon=False, loc='lower center')
    for t in leg.get_texts():
        t.set_color(INK2)
    ax.set_title(f'Estação {nome} — o que o batente do roteiro custava',
                 color=INK, fontsize=12, fontweight='bold', loc='left',
                 pad=26)
    ganho = (dados[1][5] / dados[0][5] - 1) * 100
    ax.text(0.0, 1.02,
            f'soltar o último coeficiente do intradorso vale {ganho:+.1f}% '
            f'de arrasto; o cusp côncavo no bordo de fuga é a diferença',
            transform=ax.transAxes, fontsize=9, color=INK2, va='bottom')

    # --- Cp e Mach ---
    for ax, chave, rot_y in [(axs[1], 'Cp', '$C_p$'),
                             (axs[2], 'Mach', 'Mach local')]:
        style_axes(ax)
        for rot, cor, Al, Au, r, cd in dados:
            ax.plot(r['x'], r[chave], '-', color=cor, linewidth=1.8,
                    label=rot)
        ax.set_ylabel(rot_y, color=INK2, fontsize=10)
        if chave == 'Cp':
            ax.invert_yaxis()
        else:
            ax.axhline(1.0, color=PAL[1], linestyle='--', linewidth=1.1)
            ax.annotate('sônico', (0.98, 1.0), xytext=(-4, 5),
                        textcoords='offset points', ha='right',
                        fontsize=9, color=PAL[1])
            ax.set_xlabel('$x/c$', color=INK2, fontsize=10)
    salvar(fig, os.path.join(RES, f'item9_{nome}.png'))


def main():
    for nome in ot.ESTACOES:
        print(f'--- {nome} ---')
        figura(nome)


if __name__ == '__main__':
    main()
