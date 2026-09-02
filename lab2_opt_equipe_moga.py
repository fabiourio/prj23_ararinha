'''
Lab 02 secao 4: min {W0, Wf} com NSGA-II (pymoo), mesmo problema do geom.

Populacao semeada com SLSQP interior (ancoras exatas ficam em restricao ativa
e somem por arredondamento). Frente salva em equipe_moga_frente.csv.

Rodar: python lab2_opt_equipe_moga.py  (~1-2 min)
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

from lab2_equipe_common import (DV, Xref, bounds_norm, dv_names, CON,
                                anchor_W0, anchor_Wf,
                                airplane_out, constraints_geq0)
from lab2_plot import PAL, INK, INK2, style_axes, planform

np.seterr(all='ignore')
warnings.filterwarnings('ignore')

RES = 'Resultados/3_multiobj'
os.makedirs(RES, exist_ok=True)

neval = [0]
out_ref = airplane_out(np.ones(len(DV)))
W0_ref = out_ref['W0']
Wf_ref = out_ref['Wf']


def run_analysis(x):
    out = airplane_out(x)
    neval[0] += 1
    return out


class AirplaneProblem(ElementwiseProblem):
    '''Problema pymoo: 2 objetivos, 16 desigualdades (G <= 0).'''

    def __init__(self):
        super().__init__(n_var=len(DV),
                         n_obj=2,
                         n_ieq_constr=len(CON),
                         xl=bounds_norm[:, 0],
                         xu=bounds_norm[:, 1])

    def _evaluate(self, x, out_pymoo, *args, **kwargs):
        try:
            out = run_analysis(x)
            F = [out['W0'] / W0_ref, out['Wf'] / Wf_ref]
            G = list(-constraints_geq0(out))  # pymoo: G <= 0
            if not (np.all(np.isfinite(F)) and np.all(np.isfinite(G))):
                raise ValueError('analise divergiu')
        except Exception:
            # geometria extrema do LHS pode divergir o ponto fixo do W0
            F = [10.0, 10.0]
            G = [1e3] * len(CON)

        out_pymoo['F'] = F
        out_pymoo['G'] = G


POP_SIZE = 80
N_GEN = 200


def slsqp_interior(objetivo, margem):
    '''Ponto estritamente viavel (g - margem >= 0) para semear o NSGA-II.'''

    from scipy.optimize import minimize as sp_minimize

    ref = W0_ref if objetivo == 'W0' else Wf_ref
    r = sp_minimize(lambda x: run_analysis(x)[objetivo] / ref,
                    np.ones(len(DV)),
                    constraints=[{'type': 'ineq',
                                  'fun': lambda x: constraints_geq0(run_analysis(x)) - margem}],
                    bounds=bounds_norm, method='slsqp',
                    options={'maxiter': 200, 'ftol': 1e-8})
    return r.x


problem = AirplaneProblem()
seed_W0 = slsqp_interior('W0', 0.005)
seed_Wf = slsqp_interior('Wf', 0.005)

# LHS manual: pymoo.LHS nao respeita np.random.seed(1)
rng = np.random.default_rng(1)
nv = len(DV)
u = (np.stack([rng.permutation(POP_SIZE) for _ in range(nv)], axis=1)
     + rng.random((POP_SIZE, nv))) / POP_SIZE
X0_pop = bounds_norm[:, 0] + u * (bounds_norm[:, 1] - bounds_norm[:, 0])
X0_pop[0] = np.ones(nv)
X0_pop[1] = seed_W0
X0_pop[2] = seed_Wf

algorithm = NSGA2(pop_size=POP_SIZE, sampling=X0_pop, eliminate_duplicates=True)

neval[0] = 0
t_start = time.time()
res = minimize(problem, algorithm, ('n_gen', N_GEN), seed=1, verbose=True)
t_elapsed = time.time() - t_start

order = np.argsort(res.F[:, 0])
F = res.F[order]
X = res.X[order]

W0_front = F[:, 0] * W0_ref
Wf_front = F[:, 1] * Wf_ref

print('')
print('=' * 70)
print('OTIMIZACAO MULTIOBJETIVO (NSGA-II)')
print('=' * 70)
print('Individuos por geracao (pop_size): %d' % POP_SIZE)
print('Numero de geracoes (n_gen):        %d' % N_GEN)
print('Execucoes do analyze:              %d' % neval[0])
print('Tempo de otimizacao:               %.1f s' % t_elapsed)
print('Pontos na frente de Pareto:        %d' % len(F))
print('')
print('Verificacao de convergencia com os otimos mono-objetivo (SLSQP):')
print('  min W0 da frente:  %10.1f kgf | ancora SLSQP: %10.1f kgf (%+.2f %%)'
      % (W0_front[0], anchor_W0['W0'], 100 * (W0_front[0] - anchor_W0['W0']) / anchor_W0['W0']))
print('  min Wf da frente:  %10.1f kgf | ancora SLSQP: %10.1f kgf (%+.2f %%)'
      % (Wf_front[-1], anchor_Wf['Wf'], 100 * (Wf_front[-1] - anchor_Wf['Wf']) / anchor_Wf['Wf']))

f1n = (F[:, 0] - F[:, 0].min()) / max(F[:, 0].max() - F[:, 0].min(), 1e-12)
f2n = (F[:, 1] - F[:, 1].min()) / max(F[:, 1].max() - F[:, 1].min(), 1e-12)
idxA, idxC = 0, len(F) - 1
idxB = int(np.argmin(f1n ** 2 + f2n ** 2))  # joelho na frente normalizada

sel_idx = [idxA, idxB, idxC]
sel_names = ['A (min W0)', 'B (joelho)', 'C (min Wf)']

np.savetxt(RES + '/equipe_moga_frente.csv',
           np.column_stack([W0_front, Wf_front, X * Xref]),
           header='W0_kgf Wf_kgf ' + ' '.join(dv_names), fmt='%.6g')

print('')
print('%-12s %10s %10s' % ('aeronave', 'W0 [kgf]', 'Wf [kgf]'))
for name, k in zip(sel_names, sel_idx):
    print('%-12s %10.1f %10.1f' % (name, W0_front[k], Wf_front[k]))

print('')
print('%-10s' % 'DV' + ''.join('%12s' % n for n in sel_names))
for j, dvn in enumerate(dv_names):
    print('%-10s' % dvn + ''.join('%12.4f' % (X[k, j] * Xref[j]) for k in sel_idx))

fig, ax = plt.subplots(figsize=(8.5, 6))

ax.plot(W0_front / 1000, Wf_front / 1000, '-', color=PAL[0], linewidth=1.1, alpha=0.55)
ax.plot(W0_front / 1000, Wf_front / 1000, 'o', color=PAL[0], markersize=5,
        markeredgecolor='white', markeredgewidth=1.0,
        label='frente de Pareto (NSGA-II)')

for anc, lab, dxy in ((anchor_W0, 'SLSQP min $W_0$', (-10, -16)),
                      (anchor_Wf, 'SLSQP min $W_f$', (10, -4))):
    ax.plot(anc['W0'] / 1000, anc['Wf'] / 1000, '*', color=INK, markersize=15,
            markeredgecolor='white', markeredgewidth=0.8,
            label='ótimos SLSQP (âncoras)' if anc is anchor_W0 else None)
    ax.annotate(lab, (anc['W0'] / 1000, anc['Wf'] / 1000), xytext=dxy,
                textcoords='offset points', fontsize=9, color=INK2)

sel_marks = ['o', 's', '^']
for name, k, m in zip(sel_names, sel_idx, sel_marks):
    ax.plot(W0_front[k] / 1000, Wf_front[k] / 1000, m, color=PAL[1], markersize=10,
            markeredgecolor='white', markeredgewidth=1.2, zorder=5,
            label='selecionadas (A, B, C)' if m == 'o' else None)
    ax.annotate(name.split()[0], (W0_front[k] / 1000, Wf_front[k] / 1000),
                xytext=(8, 7), textcoords='offset points',
                fontsize=10, color=INK, fontweight='bold')

ax.set_xlabel('$W_0$ [t]', fontsize=13)
ax.set_ylabel('$W_f$ [t]', fontsize=13)
ax.margins(x=0.10, y=0.12)
ax.legend(fontsize=9, frameon=False, loc='lower left')
style_axes(ax)

axi = ax.inset_axes([0.58, 0.58, 0.39, 0.38])
axi.plot(W0_front / 1000, Wf_front / 1000, 'o', color=PAL[0], markersize=2.5)
axi.plot(out_ref['W0'] / 1000, out_ref['Wf'] / 1000, 's', color=INK, markersize=6)
axi.annotate('baseline PRJ-22', (out_ref['W0'] / 1000, out_ref['Wf'] / 1000),
             xytext=(-8, -14), textcoords='offset points',
             fontsize=8, color=INK2, ha='right')
axi.set_title('contexto completo', fontsize=8, color=INK2)
axi.tick_params(labelsize=7)
style_axes(axi)

plt.tight_layout()
fig.savefig(RES + '/equipe_moga_pareto.png', dpi=150)

fig = plt.figure(figsize=(12, 6.5))
gs = fig.add_gridspec(2, 2, width_ratios=[2.4, 1], hspace=0.4, wspace=0.15)
axp = fig.add_subplot(gs[:, 0])
axb1 = fig.add_subplot(gs[0, 1])
axb2 = fig.add_subplot(gs[1, 1])

for name, k, c in zip(sel_names, sel_idx, PAL):
    airplane = standard_airplane('my_airplane')
    for j, dvn in enumerate(dv_names):
        airplane['inputs'][dvn] = X[k, j] * Xref[j]
    analyze(airplane)
    planform(axp, airplane, c, name)

axp.set_title('Planformas das aeronaves da frente de Pareto', fontsize=13)
axp.set_xlabel('x [m]', fontsize=12)
axp.set_ylabel('y [m]', fontsize=12)
axp.set_aspect('equal')
axp.legend(fontsize=9, frameon=False)
style_axes(axp)

letras = ['A', 'B', 'C']
jS = dv_names.index('S_w')
jT = dv_names.index('tcr_w')
Svals = [X[k, jS] * Xref[jS] for k in sel_idx]
Tvals = [X[k, jT] * Xref[jT] for k in sel_idx]

for axb, vals, titulo, fmt in ((axb1, Svals, '$S_w$ [m$^2$]', '%.0f'),
                               (axb2, Tvals, '$(t/c)_{r,w}$', '%.3f')):
    ypos = np.arange(len(letras))[::-1]
    for yi, v, c in zip(ypos, vals, PAL):
        axb.plot(v, yi, 'o', color=c, markersize=9,
                 markeredgecolor='white', markeredgewidth=1.0)
        axb.annotate(fmt % v, (v, yi), xytext=(0, 9),
                     textcoords='offset points', ha='center',
                     fontsize=9, color=INK2)
    axb.set_yticks(ypos)
    axb.set_yticklabels(letras, fontsize=10)
    vmin, vmax = min(vals), max(vals)
    pad = 0.35 * (vmax - vmin) + 1e-9
    axb.set_xlim(vmin - pad, vmax + pad)
    axb.set_ylim(-0.6, 2.9)
    axb.set_title(titulo, fontsize=11, color=INK)
    style_axes(axb)

plt.tight_layout()
fig.savefig(RES + '/equipe_moga_planformas.png', dpi=150)

plt.show()
