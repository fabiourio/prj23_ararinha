'''
Lab 02 - Topico 3 (storytelling): evidencia de convergencia do NSGA-II.

Roda uma versao "ingenua" do problema multiobjetivo -- populacao inicial
aleatoria, sem sementes SLSQP, com metade do orcamento (pop 60 x 120
geracoes) -- e compara com a frente final do lab2_opt_equipe_moga.py
(lida do CSV) e com os otimos mono-objetivo (ancoras).

A frente ingenua estagna dominada pelas ancoras; a frente final, semeada
e com orcamento maior, encosta nelas. E' a resposta da pergunta 3 da
Secao 4 do roteiro em forma de figura.

Gera Resultados/3_multiobj/equipe_moga_convergencia.png
'''

import os
import time
import warnings

import numpy as np
import matplotlib.pyplot as plt
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.constants import gravity

rad2deg = 180/np.pi

np.seterr(all='ignore')
warnings.filterwarnings('ignore')

RES = 'Resultados/3_multiobj'
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
Xref = np.array([d[1] for d in DV])
bounds_phys = np.array([[d[2], d[3]] for d in DV])
bounds_norm = np.sort(bounds_phys/Xref[:,None], axis=1)

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

anchor_W0 = {'W0': 290117.1, 'Wf': 111036.5}
anchor_Wf = {'W0': 291961.0, 'Wf': 109117.9}

def run_analysis(x):

    X = np.asarray(x)*Xref

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

def constraints_geq0(out):
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

out_ref = run_analysis(np.ones(len(DV)))
W0_ref = out_ref['W0']
Wf_ref = out_ref['Wf']

class AirplaneProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(n_var=len(DV), n_obj=2, n_ieq_constr=len(CON),
                         xl=bounds_norm[:,0], xu=bounds_norm[:,1])

    def _evaluate(self, x, out_pymoo, *args, **kwargs):
        try:
            out = run_analysis(x)
            F = [out['W0']/W0_ref, out['Wf']/Wf_ref]
            G = list(-constraints_geq0(out))
            if not (np.all(np.isfinite(F)) and np.all(np.isfinite(G))):
                raise ValueError('analise divergiu')
        except Exception:
            F = [10.0, 10.0]
            G = [1e3]*len(CON)
        out_pymoo['F'] = F
        out_pymoo['G'] = G

# rodada "ingenua": populacao inicial aleatoria (sem sementes SLSQP) e
# metade do orcamento; RNG proprio para manter a reprodutibilidade
POP, NGEN = 60, 120
rng = np.random.default_rng(2)
X0_pop = bounds_norm[:,0] + rng.random((POP, len(DV)))*(bounds_norm[:,1] - bounds_norm[:,0])

problem = AirplaneProblem()
algorithm = NSGA2(pop_size=POP, sampling=X0_pop, eliminate_duplicates=True)

t0 = time.time()
res = minimize(problem, algorithm, ('n_gen', NGEN), seed=1, verbose=False)
print('Rodada ingenua: %d pontos, %.0f s' % (len(res.F), time.time() - t0))

W0_naive = res.F[:,0]*W0_ref
Wf_naive = res.F[:,1]*Wf_ref

# frente final (script principal)
dados = np.loadtxt(RES + '/equipe_moga_frente.csv')
W0_final, Wf_final = dados[:,0], dados[:,1]

print('min W0: ingenua %.1f | final %.1f | ancora %.1f'
      % (W0_naive.min(), W0_final.min(), anchor_W0['W0']))
print('min Wf: ingenua %.1f | final %.1f | ancora %.1f'
      % (Wf_naive.min(), Wf_final.min(), anchor_Wf['Wf']))

### FIGURA

PAL = ['#2a78d6', '#eb6834']
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#c3c2b7'
GRID = '#e1e0d9'

fig, ax = plt.subplots(figsize=(8.5, 6))

o = np.argsort(W0_naive)
ax.plot(W0_naive[o]/1000, Wf_naive[o]/1000, 'o-', color=MUTED, markersize=7,
        linewidth=1.0, markeredgecolor='white', markeredgewidth=0.8,
        label='NSGA-II ingênuo (pop 60 × 120 ger., sem sementes)')
ax.annotate('estagnou dominado:\n%d ponto(s), %+.1f%% da âncora'
            % (len(W0_naive), 100*(W0_naive.min()-anchor_W0['W0'])/anchor_W0['W0']),
            (W0_naive.min()/1000, Wf_naive[np.argmin(W0_naive)]/1000),
            xytext=(-12, -26), textcoords='offset points',
            fontsize=9, color=INK2, ha='right')

o = np.argsort(W0_final)
ax.plot(W0_final[o]/1000, Wf_final[o]/1000, 'o-', color=PAL[0], markersize=5,
        linewidth=1.0, markeredgecolor='white', markeredgewidth=0.8,
        label='NSGA-II final (pop 80 × 200 ger., semeado)')

for anc, lab, dxy in ((anchor_W0, 'SLSQP min $W_0$', (-10, -16)),
                      (anchor_Wf, 'SLSQP min $W_f$', (6, -16))):
    ax.plot(anc['W0']/1000, anc['Wf']/1000, '*', color=INK, markersize=15,
            markeredgecolor='white', markeredgewidth=0.8,
            label='ótimos SLSQP (âncoras)' if anc is anchor_W0 else None)
    ax.annotate(lab, (anc['W0']/1000, anc['Wf']/1000), xytext=dxy,
                textcoords='offset points', fontsize=9, color=INK2)

ax.set_xlabel('$W_0$ [t]', fontsize=13)
ax.set_ylabel('$W_f$ [t]', fontsize=13)
ax.margins(x=0.08, y=0.12)
ax.legend(fontsize=9, frameon=False, loc='lower left')
ax.grid(color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
for side in ('left', 'bottom'):
    ax.spines[side].set_color(MUTED)
ax.tick_params(colors=INK2, labelsize=9)

ax.set_title('Convergência do NSGA-II verificada pelos ótimos mono-objetivo',
             fontsize=12, color=INK)

plt.tight_layout()
fig.savefig(RES + '/equipe_moga_convergencia.png', dpi=150)

plt.show()
