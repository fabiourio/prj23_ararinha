'''
Lab 02 - Topico 3 (sintese): qual aeronave escolher, e com base em que.

Mostra o processo completo no plano W0 x Wf: v1 (PRJ-22) -> otimizacao
mono-objetivo (v2, min W0) -> frente de Pareto -> escolha final.

Criterio de escolha: a taxa marginal de troca ao longo da frente. Do
extremo de min W0 ate o joelho (B), cada kg de MTOW aceito compra ~1 kg
de combustivel POR MISSAO -- troca vantajosa, ja que o combustivel se
paga a cada voo. Depois do joelho a taxa piora para ~2,3 kg de MTOW por
kg de combustivel, e parar em B se justifica.

Gera Resultados/3_multiobj/equipe_escolha_aeronave.png
'''

import os

import numpy as np
import matplotlib.pyplot as plt

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.constants import gravity

RES = 'Resultados/3_multiobj'
os.makedirs(RES, exist_ok=True)

# v1: baseline do PRJ-22, calculada na hora
ap = standard_airplane('my_airplane')
analyze(ap)
v1 = {'W0': ap['thrust_matching']['W0']/gravity,
      'Wf': ap['thrust_matching']['W_fuel']/gravity}

# v2: otimo mono-objetivo (lab2_opt_equipe_geom.py)
v2 = {'W0': 290117.1, 'Wf': 111036.5}

# frente de Pareto e selecionadas (lab2_opt_equipe_moga.py)
dados = np.loadtxt(RES + '/equipe_moga_frente.csv')
W0f, Wff = dados[:,0], dados[:,1]
A = {'W0': 290758.2, 'Wf': 110994.1}
B = {'W0': 291291.1, 'Wf': 109891.3}
C = {'W0': 292468.2, 'Wf': 109385.4}

taxa_v2B = (B['W0'] - v2['W0'])/(v2['Wf'] - B['Wf'])
taxa_BC = (C['W0'] - B['W0'])/(B['Wf'] - C['Wf'])

print('taxa de troca v2->B: %.2f kg MTOW por kg de combustivel' % taxa_v2B)
print('taxa de troca B->C : %.2f kg MTOW por kg de combustivel' % taxa_BC)

PAL = ['#2a78d6', '#eb6834', '#1baf7a']
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#c3c2b7'
GRID = '#e1e0d9'

def style_axes(ax):
    ax.grid(color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9)

fig, ax = plt.subplots(figsize=(9, 6.5))

# frente
o = np.argsort(W0f)
ax.plot(W0f[o]/1000, Wff[o]/1000, '-', color=PAL[0], linewidth=1.1, alpha=0.5)
ax.plot(W0f/1000, Wff/1000, 'o', color=PAL[0], markersize=4.5,
        markeredgecolor='white', markeredgewidth=0.9,
        label='frente de Pareto')

# v2 e selecionadas
ax.plot(v2['W0']/1000, v2['Wf']/1000, '*', color=INK, markersize=16,
        markeredgecolor='white', markeredgewidth=0.8,
        label='ótimo mono-objetivo (SLSQP)')
ax.annotate('min $W_0$ (SLSQP)', (v2['W0']/1000, v2['Wf']/1000),
            xytext=(-10, 8), textcoords='offset points', ha='right',
            fontsize=9, color=INK2)

for p, nome, dxy in ((A, 'A', (8, 6)), (C, 'C', (8, 6))):
    ax.plot(p['W0']/1000, p['Wf']/1000, 'o', color=PAL[1], markersize=9,
            markeredgecolor='white', markeredgewidth=1.2, zorder=5)
    ax.annotate(nome, (p['W0']/1000, p['Wf']/1000), xytext=dxy,
                textcoords='offset points', fontsize=10, color=INK,
                fontweight='bold')

# escolha: B, com destaque
ax.plot(B['W0']/1000, B['Wf']/1000, 'o', color=PAL[1], markersize=13,
        markeredgecolor=INK, markeredgewidth=1.6, zorder=6,
        label='B: escolha da equipe (joelho)')
ax.annotate('B — escolha da equipe', (B['W0']/1000, B['Wf']/1000),
            xytext=(12, 10), textcoords='offset points', fontsize=10,
            color=INK, fontweight='bold')

# taxas marginais de troca anotadas nos trechos
ax.annotate('até o joelho:\n≈ %.0f kg de MTOW\npor kg de combustível' % taxa_v2B,
            ((v2['W0'] + B['W0'])/2000, (v2['Wf'] + B['Wf'])/2000),
            xytext=(-84, -46), textcoords='offset points',
            fontsize=8.5, color=INK2)
ax.annotate('depois do joelho:\n≈ %.1f kg de MTOW\npor kg de combustível' % taxa_BC,
            ((B['W0'] + C['W0'])/2000, (B['Wf'] + C['Wf'])/2000),
            xytext=(10, 16), textcoords='offset points',
            fontsize=8.5, color=INK2)

ax.set_xlabel('$W_0$ [t]', fontsize=13)
ax.set_ylabel('$W_f$ [t]', fontsize=13)
ax.margins(x=0.10, y=0.14)
ax.legend(fontsize=9, frameon=False, loc='upper right')
style_axes(ax)

# inset: o caminho desde a v1
axi = ax.inset_axes([0.045, 0.06, 0.34, 0.36])
axi.plot(W0f/1000, Wff/1000, 'o', color=PAL[0], markersize=2)
axi.plot(v1['W0']/1000, v1['Wf']/1000, 's', color=INK, markersize=6)
axi.plot(v2['W0']/1000, v2['Wf']/1000, '*', color=INK, markersize=9)
axi.annotate('', xy=(v2['W0']/1000, v2['Wf']/1000),
             xytext=(v1['W0']/1000, v1['Wf']/1000),
             arrowprops=dict(arrowstyle='->', color=INK2, lw=1.1))
axi.annotate('baseline (PRJ-22)', (v1['W0']/1000, v1['Wf']/1000),
             xytext=(-6, -12), textcoords='offset points', ha='right',
             fontsize=7.5, color=INK2)
axi.annotate('$-5{,}4\\%$ $W_0$\n$-9{,}5\\%$ $W_f$',
             ((v1['W0'] + v2['W0'])/2000, (v1['Wf'] + v2['Wf'])/2000),
             xytext=(2, 6), textcoords='offset points',
             fontsize=7.5, color=INK2)
axi.set_title('o caminho desde a baseline', fontsize=8, color=INK2)
axi.tick_params(labelsize=7)
axi.margins(x=0.25, y=0.25)
style_axes(axi)

ax.set_title('Da baseline à escolha final: mono-objetivo aproxima, a frente decide',
             fontsize=12.5, color=INK)

plt.tight_layout()
fig.savefig(RES + '/equipe_escolha_aeronave.png', dpi=150)

plt.show()
