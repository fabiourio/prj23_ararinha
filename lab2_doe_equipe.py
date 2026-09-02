'''
DOE e figuras de apoio da secao 3: cortes no otimo, etapas de W0, delta v1->v2.

Rodar: python lab2_doe_equipe.py
'''

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from lab2_equipe_common import (DV, Xref, Xopt, CON, dv_names,
                                airplane_out, constraints_geq0)
from lab2_plot import PAL, INK, INK2, style_axes

RES = 'Resultados/2_monoobj_equipe'
os.makedirs(RES, exist_ok=True)

W0_opt = 290117.1

Xetapa2 = dict(AR_w=10.7327, xr_w=16.1393, S_w=361.5249,
               sweep_w=0.5949, x_mlg=31.2840, tcr_w=0.2044)


def run_analysis(X):
    x_norm = X / Xref
    return airplane_out(x_norm)


def g_norm(out):
    return constraints_geq0(out)


NPT = 15
TOL_VIAVEL = 1e-3  # folga de 0,1% na g normalizada
unidades = {'AR_w': '', 'xr_w': ' [m]', 'S_w': ' [m$^2$]', 'sweep_w': ' [rad]',
            'x_mlg': ' [m]', 'tcr_w': '', 'z_lg': ' [m]', 'y_mlg': ' [m]'}

cor_con = {'SM_fwd': '#2a78d6', 'alpha_tipback': '#eb6834',
           'alpha_tailstrike': '#1baf7a', 'tank_excess': '#eda100',
           'mlg_track': '#e87ba4', 'mlg_fit': '#008300',
           'SM_aft': '#4a3aa7', 'ground_clearance': '#e34948'}
COR_OUTRAS = '#52514e'
COR_VIAVEL = '#0ca30c'
MUTED = '#c3c2b7'

fig, axs = plt.subplots(2, 4, figsize=(13, 6.5), sharey=True)

neval = 0
usados = set()
for k, (name, _, lb, ub) in enumerate(DV):
    ax = axs.flat[k]
    xs = np.linspace(lb, ub, NPT)
    W0s = np.full(NPT, np.nan)
    classe = [None] * NPT
    for j, xv in enumerate(xs):
        X = Xopt.copy()
        X[k] = xv
        try:
            out = run_analysis(X)
            gg = g_norm(out)
            W0s[j] = out['W0'] / 1000
            if gg.min() >= -TOL_VIAVEL:
                classe[j] = 'viável'
            else:
                nome_con = CON[int(np.argmin(gg))][0]
                classe[j] = nome_con if nome_con in cor_con else 'outras'
            usados.add(classe[j])
            neval += 1
        except Exception:
            pass

    ax.plot(xs, W0s, '-', color=MUTED, linewidth=1.0, zorder=1)
    for cl in set(c for c in classe if c is not None):
        cor = (COR_VIAVEL if cl == 'viável'
               else COR_OUTRAS if cl == 'outras' else cor_con[cl])
        marca = 's' if cl == 'viável' else 'o'
        sel = np.array([c == cl for c in classe])
        ax.plot(xs[sel], W0s[sel], marca, color=cor, markersize=5.5,
                markeredgecolor='white', markeredgewidth=0.8, zorder=3)
    ax.plot([Xopt[k]], [W0_opt / 1000], '*', color=INK, markersize=13,
            markeredgecolor='white', markeredgewidth=0.6, zorder=4)
    ax.set_title(name + unidades[name], fontsize=11, color=INK)
    style_axes(ax)

for ax in axs[:, 0]:
    ax.set_ylabel('$W_0$ [t]', fontsize=11)

ordem_leg = ['viável'] + [c for c in cor_con if c in usados] + \
            (['outras'] if 'outras' in usados else [])
proxies = []
for cl in ordem_leg:
    cor = (COR_VIAVEL if cl == 'viável'
           else COR_OUTRAS if cl == 'outras' else cor_con[cl])
    rotulo = cl if cl == 'viável' else 'bloqueia: ' + cl
    proxies.append(Line2D([], [], marker='s' if cl == 'viável' else 'o',
                          linestyle='', color=cor,
                          markeredgecolor='white', label=rotulo))
proxies.append(Line2D([], [], marker='*', linestyle='', color=INK,
                      markersize=12, label='ótimo (SLSQP)'))
fig.legend(handles=proxies, loc='lower center', ncol=min(5, len(proxies)),
           frameon=False, fontsize=9)

fig.suptitle('Cortes 1-a-1 no ótimo: qual restrição barra cada direção',
             fontsize=13, color=INK)
plt.tight_layout(rect=[0, 0.07, 1, 0.97])
fig.savefig(RES + '/equipe_doe_sensibilidade.png', dpi=150)
print('Cortes 1-a-1: %d analises' % neval)

Xbase = Xref.copy()
out_base = run_analysis(Xbase)
X6 = Xbase.copy()
for name, value in Xetapa2.items():
    X6[dv_names.index(name)] = value
out_6dv = run_analysis(X6)
out_opt = run_analysis(Xopt)

etapas = [('Baseline PRJ-22', out_base, False,
           'trem principal atrás do bordo de fuga (mlg_fit = %.3f)' % out_base['mlg_fit']),
          ('Otimização 6 DVs\n(sem encaixes do trem)', out_6dv, False,
           'ótimo aparente, mas trem ainda atrás da asa (mlg_fit = %.3f)' % out_6dv['mlg_fit']),
          ('Formulação final 8 DVs\n(16 restrições)', out_opt, True,
           'viável: encaixes do trem custam +%.1f t' % ((out_opt['W0'] - out_6dv['W0']) / 1000))]

fig, ax = plt.subplots(figsize=(9, 4.2))

ypos = np.arange(len(etapas))[::-1]
for yi, (nome, out, ok, nota) in zip(ypos, etapas):
    w0 = out['W0'] / 1000
    if ok:
        ax.plot(w0, yi, 'o', color=PAL[0], markersize=11,
                markeredgecolor='white', markeredgewidth=1.2, zorder=3)
    else:
        ax.plot(w0, yi, 'o', markerfacecolor='white', markeredgecolor=PAL[0],
                markersize=11, markeredgewidth=1.6, zorder=3)
    ax.annotate('%.1f t' % w0, (w0, yi), xytext=(0, 12),
                textcoords='offset points', ha='center',
                fontsize=10, color=INK, fontweight='bold')
    ax.annotate(nota, (w0, yi), xytext=(0, -18),
                textcoords='offset points', ha='center',
                fontsize=8.5, color=INK2)

ax.set_yticks(ypos)
ax.set_yticklabels([e[0] for e in etapas], fontsize=10)
ax.set_xlabel('$W_0$ [t]', fontsize=12)
ax.set_xlim(285, 310)
ax.set_ylim(-0.7, 2.7)
ax.set_title('MTOW nas etapas da otimização (marcador vazado = inviável)',
             fontsize=12, color=INK)
style_axes(ax)

plt.tight_layout()
fig.savefig(RES + '/equipe_etapas_w0.png', dpi=150)
print('Etapas: baseline %.1f t | 6 DVs %.1f t | final %.1f t'
      % (out_base['W0'] / 1000, out_6dv['W0'] / 1000, out_opt['W0'] / 1000))

COR_MAIS = '#e34948'
COR_MENOS = '#2a78d6'


def pct(v1, v2):
    return 100.0 * (v2 - v1) / abs(v1)


g_base = out_base
g_opt2 = out_opt

linhas_dv = [
    ('$y_{mlg}$',        '%.2f m',  5.50,   6.95),
    ('$AR_w$',           '%.2f',    8.00,   9.80),
    ('$(t/c)_{r,w}$',    '%.3f',    0.180,  0.2093),
    ('$\\Lambda_w$',     '%.3f rad', 0.580, 0.6108),
    ('trem (comprim. $-z_{lg}$)', '%.2f m', 5.75, 6.025),
    ('$x_{r,w}$',        '%.2f m', 17.00,  16.011),
    ('$x_{mlg}$',        '%.2f m', 31.00,  29.263),
    ('$S_w$',            '%.1f m$^2$', 390.0, 353.67),
]
linhas_res = [
    ('$b_w$',      '%.1f m', g_base['b_w'],          g_opt2['b_w']),
    ('$W_0$',      '%.1f t', g_base['W0'] / 1000,      g_opt2['W0'] / 1000),
    ('$W_{vazio}$','%.1f t', g_base['W_empty'] / 1000, g_opt2['W_empty'] / 1000),
    ('$W_f$',      '%.1f t', g_base['Wf'] / 1000,      g_opt2['Wf'] / 1000),
]

fig, (axd, axr) = plt.subplots(2, 1, figsize=(9, 7.2), sharex=True,
                               gridspec_kw={'height_ratios': [2, 1.05]})

for ax, linhas, titulo in ((axd, linhas_dv, 'Variáveis de projeto'),
                           (axr, linhas_res, 'Resultados')):
    nomes = [l[0] for l in linhas]
    deltas = [pct(l[2], l[3]) for l in linhas]
    ypos = np.arange(len(linhas))[::-1]
    cores = [COR_MAIS if d > 0 else COR_MENOS for d in deltas]
    ax.barh(ypos, deltas, height=0.55, color=cores)
    for yi, d, (nome, fmt, v1, v2) in zip(ypos, deltas, linhas):
        lado = 1 if d > 0 else -1
        ax.annotate('%+.1f%%' % d, (d, yi), xytext=(6 * lado, 0),
                    textcoords='offset points', va='center',
                    ha='left' if d > 0 else 'right',
                    fontsize=9, color=INK, fontweight='bold')
        ax.annotate(('v1 ' + fmt + ' → v2 ' + fmt) % (v1, v2), (0, yi),
                    xytext=(-6 * lado, 0), textcoords='offset points',
                    va='center', ha='right' if d > 0 else 'left',
                    fontsize=8, color=INK2)
    ax.axvline(0, color=INK2, linewidth=0.9)
    ax.set_yticks(ypos)
    ax.set_yticklabels(nomes, fontsize=10)
    ax.set_title(titulo, fontsize=11, color=INK, loc='left')
    ax.set_xlim(-32, 32)
    style_axes(ax)

axr.set_xlabel('variação v1 → v2 [%]', fontsize=12)
fig.suptitle('Aeronave da equipe: v1 (PRJ-22) → v2 (otimizada, $-5{,}42\\%$ de MTOW)',
             fontsize=13, color=INK)

plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(RES + '/equipe_v1v2_delta.png', dpi=150)
print('v1 -> v2: W0 %+.2f%% | Wf %+.2f%% | W_empty %+.2f%% | b_w %+.2f%%'
      % (pct(g_base['W0'], g_opt2['W0']), pct(g_base['Wf'], g_opt2['Wf']),
         pct(g_base['W_empty'], g_opt2['W_empty']), pct(g_base['b_w'], g_opt2['b_w'])))

plt.show()
