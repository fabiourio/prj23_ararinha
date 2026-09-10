'''
Lab 02 - Topico 2 (storytelling): exploracao do espaco de projeto da
aeronave da equipe, seguindo o passo-a-passo das aulas (DOE antes de
confiar no otimo).

Gera duas figuras em Resultados/2_monoobj_equipe/:

1. equipe_doe_sensibilidade.png -- cortes 1-a-1 em torno do otimo: cada
   painel varre uma DV com as demais fixas no otimo, mostrando a
   sensibilidade de W0 e a janela viavel que as 16 restricoes recortam.
   Nas DVs "de posicao" (xr_w, x_mlg, z_lg, y_mlg) a curva e' plana --
   elas nao movem o objetivo, sao posicionadas pelas restricoes.

2. equipe_etapas_w0.png -- W0 nas etapas da otimizacao: baseline PRJ-22,
   otimizacao com 6 DVs (sem as restricoes de encaixe do trem) e
   formulacao final com 8 DVs. A versao de 6 DVs parece melhor, mas deixa
   o trem principal atras do bordo de fuga da asa; as restricoes de
   realismo custam ~1,9 t.
'''

import os

import numpy as np
import matplotlib.pyplot as plt

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.constants import gravity

rad2deg = 180/np.pi

RES = 'Resultados/2_monoobj_equipe'
os.makedirs(RES, exist_ok=True)

DV = [('AR_w',      8.0,   7.0,  12.0),
      ('xr_w',     17.0,  14.0,  20.0),
      ('S_w',     390.0, 330.0, 460.0),
      ('sweep_w',  0.58,  0.40,  0.70),
      ('x_mlg',    31.0,  28.0,  35.0),
      ('tcr_w',     0.18,  0.12,  0.24),
      ('z_lg',     -5.75, -7.50, -4.50),
      ('y_mlg',     5.50,  4.00,   9.00)]

dv_names = [d[0] for d in DV]
Xbase = np.array([d[1] for d in DV])

CON = [('deltaS_wlan',      '>=',  0.0),
       ('SM_fwd',           '<=',  0.30),
       ('SM_aft',           '>=',  0.05),
       ('frac_nlg_fwd',     '<=',  0.18),
       ('frac_nlg_aft',     '>=',  0.03),
       ('alpha_tipback',    '>=', 15.0),
       ('alpha_tailstrike', '>=', 10.0),
       ('phi_overturn',     '<=', 63.0),
       ('ground_clearance', '>=',  0.5),
       ('tank_excess',      '>=',  0.0),
       ('b_w',              '<=', 64.9),
       ('mlg_track',        '<=', 13.9),
       ('h_tail',           '<=', 20.0),
       ('T0req',            '<=',  1.0),
       ('vt_fit',           '<=',  1.0),
       ('mlg_fit',          '<=',  1.0)]

# otimo da formulacao final (saida do lab2_opt_equipe_geom.py)
Xopt = np.array([9.8001, 16.0113, 353.6712, 0.6108,
                 29.2629, 0.2093, -6.0247, 6.9500])
W0_opt = 290117.1

# otimo da etapa de 6 DVs (lab2_opt_equipe.py, sem encaixes do trem)
Xetapa2 = dict(AR_w=10.7327, xr_w=16.1393, S_w=361.5249,
               sweep_w=0.5949, x_mlg=31.2840, tcr_w=0.2044)

def run_analysis(X):

    airplane = standard_airplane('my_airplane_v1')
    for name, value in zip(dv_names, X):
        airplane['inputs'][name] = value

    analyze(airplane)

    g = airplane['geometry']
    b = airplane['balance']
    lg = airplane['landing_gear']
    tm = airplane['thrust_matching']
    i = airplane['inputs']

    eta = i['y_mlg']/(g['b_w']/2)
    c_mlg = g['cr_w'] - (g['cr_w'] - g['ct_w'])*eta
    te_mlg = i['xr_w'] + (g['xt_w'] - i['xr_w'])*eta + c_mlg

    out = {'W0'               : tm['W0']/gravity,
           'Wf'               : tm['W_fuel']/gravity,
           'W_empty'          : tm['W_empty']/gravity,
           'deltaS_wlan'      : tm['deltaS_wlan']/i['S_w'],
           'SM_fwd'           : b['SM_fwd'],
           'SM_aft'           : b['SM_aft'],
           'frac_nlg_fwd'     : lg['frac_nlg_fwd'],
           'frac_nlg_aft'     : lg['frac_nlg_aft'],
           'alpha_tipback'    : lg['alpha_tipback']*rad2deg,
           'alpha_tailstrike' : lg['alpha_tailstrike']*rad2deg,
           'phi_overturn'     : lg['phi_overturn']*rad2deg,
           'ground_clearance' : lg['ground_clearance'],
           'tank_excess'      : b['tank_excess'],
           'b_w'              : g['b_w'],
           'mlg_track'        : lg['mlg_track'],
           'h_tail'           : (i['zr_v'] + g['b_v']) - i['z_lg'],
           'T0req'            : max(tm['T0req'].values())/tm['T0'],
           'vt_fit'           : (g['xr_v'] + g['cr_v'])/i['L_f'],
           'mlg_fit'          : i['x_mlg']/te_mlg}

    return out

def g_norm(out):
    gg = []
    for name, sense, lim in CON:
        v = out[name]
        if lim == 0.0:
            gg.append(v if sense == '>=' else -v)
        elif sense == '>=':
            gg.append(v/lim - 1)
        else:
            gg.append(1 - v/lim)
    return np.array(gg)

# paleta e tintas (mesma convencao dos demais scripts)
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
    ax.tick_params(colors=INK2, labelsize=8)

### FIGURA 1: cortes 1-a-1 em torno do otimo
# cada ponto e' colorido pela restricao que bloqueia aquela direcao (a mais
# violada); com 6 restricoes ativas no otimo, qualquer passo em uma DV
# isolada sai do conjunto viavel -- o otimo esta "encurralado", e o corte
# mostra QUEM barra cada direcao. Cores iguais as do historico.

NPT = 15
TOL_VIAVEL = 1e-3   # 0,1% de folga normalizada
unidades = {'AR_w': '', 'xr_w': ' [m]', 'S_w': ' [m$^2$]', 'sweep_w': ' [rad]',
            'x_mlg': ' [m]', 'tcr_w': '', 'z_lg': ' [m]', 'y_mlg': ' [m]'}

# mesma atribuicao de cores do lab2_opt_equipe_geom (ativas na ordem da CON)
cor_con = {'SM_fwd': '#2a78d6', 'alpha_tipback': '#eb6834',
           'alpha_tailstrike': '#1baf7a', 'tank_excess': '#eda100',
           'mlg_track': '#e87ba4', 'mlg_fit': '#008300',
           'SM_aft': '#4a3aa7', 'ground_clearance': '#e34948'}
COR_OUTRAS = '#52514e'
COR_VIAVEL = '#0ca30c'

fig, axs = plt.subplots(2, 4, figsize=(13, 6.5), sharey=True)

neval = 0
usados = set()
for k, (name, _, lb, ub) in enumerate(DV):
    ax = axs.flat[k]
    xs = np.linspace(lb, ub, NPT)
    W0s = np.full(NPT, np.nan)
    classe = [None]*NPT
    for j, xv in enumerate(xs):
        X = Xopt.copy()
        X[k] = xv
        try:
            out = run_analysis(X)
            gg = g_norm(out)
            W0s[j] = out['W0']/1000
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
        # viavel ganha forma propria (quadrado) alem da cor de status
        marca = 's' if cl == 'viável' else 'o'
        sel = np.array([c == cl for c in classe])
        ax.plot(xs[sel], W0s[sel], marca, color=cor, markersize=5.5,
                markeredgecolor='white', markeredgewidth=0.8, zorder=3)
    ax.plot([Xopt[k]], [W0_opt/1000], '*', color=INK, markersize=13,
            markeredgecolor='white', markeredgewidth=0.6, zorder=4)
    ax.set_title(name + unidades[name], fontsize=11, color=INK)
    style_axes(ax)

for ax in axs[:, 0]:
    ax.set_ylabel('$W_0$ [t]', fontsize=11)

from matplotlib.lines import Line2D
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

### FIGURA 2: W0 por etapa da otimizacao

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
           'viável: encaixes do trem custam +%.1f t' % ((out_opt['W0'] - out_6dv['W0'])/1000))]

fig, ax = plt.subplots(figsize=(9, 4.2))

ypos = np.arange(len(etapas))[::-1]
for yi, (nome, out, ok, nota) in zip(ypos, etapas):
    w0 = out['W0']/1000
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
      % (out_base['W0']/1000, out_6dv['W0']/1000, out_opt['W0']/1000))

### FIGURA 3: comparacao v1 (baseline PRJ-22) x v2 (otimizada)
# variacao percentual das DVs e dos resultados; dupla divergente azul/vermelho
# (aumento = vermelho, reducao = azul), sem juizo de valor -- e' so o sinal.
# z_lg vira comprimento do trem (-z_lg) para o sinal ter sentido fisico.

COR_MAIS = '#e34948'
COR_MENOS = '#2a78d6'

def pct(v1, v2):
    return 100.0*(v2 - v1)/abs(v1)

g_base = run_analysis(Xbase)
g_opt2 = out_opt   # ja calculado acima

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
    ('$W_0$',      '%.1f t', g_base['W0']/1000,      g_opt2['W0']/1000),
    ('$W_{vazio}$','%.1f t', g_base['W_empty']/1000, g_opt2['W_empty']/1000),
    ('$W_f$',      '%.1f t', g_base['Wf']/1000,      g_opt2['Wf']/1000),
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
        ax.annotate('%+.1f%%' % d, (d, yi), xytext=(6*lado, 0),
                    textcoords='offset points', va='center',
                    ha='left' if d > 0 else 'right',
                    fontsize=9, color=INK, fontweight='bold')
        ax.annotate(('v1 ' + fmt + ' → v2 ' + fmt) % (v1, v2), (0, yi),
                    xytext=(-6*lado, 0), textcoords='offset points',
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
