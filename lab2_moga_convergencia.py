'''
Compara NSGA-II ingenuo (sem sementes) com a frente final e ancoras SLSQP.

Precisa de equipe_moga_frente.csv gerado por lab2_opt_equipe_moga.py.
Rodar: python lab2_moga_convergencia.py
'''

import os
import time
import warnings

import numpy as np
import matplotlib.pyplot as plt
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem

from lab2_equipe_common import (DV, Xref, bounds_norm, CON,
                                anchor_W0, anchor_Wf,
                                airplane_out, constraints_geq0)
from lab2_plot import PAL, INK, INK2, MUTED, GRID, style_axes

np.seterr(all='ignore')
warnings.filterwarnings('ignore')

RES = 'Resultados/3_multiobj'
os.makedirs(RES, exist_ok=True)

out_ref = airplane_out(np.ones(len(DV)))
W0_ref = out_ref['W0']
Wf_ref = out_ref['Wf']


class AirplaneProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(n_var=len(DV), n_obj=2, n_ieq_constr=len(CON),
                         xl=bounds_norm[:, 0], xu=bounds_norm[:, 1])

    def _evaluate(self, x, out_pymoo, *args, **kwargs):
        try:
            out = airplane_out(x)
            F = [out['W0'] / W0_ref, out['Wf'] / Wf_ref]
            G = list(-constraints_geq0(out))
            if not (np.all(np.isfinite(F)) and np.all(np.isfinite(G))):
                raise ValueError('analise divergiu')
        except Exception:
            F = [10.0, 10.0]
            G = [1e3] * len(CON)
        out_pymoo['F'] = F
        out_pymoo['G'] = G


POP, NGEN = 60, 120
rng = np.random.default_rng(2)
X0_pop = bounds_norm[:, 0] + rng.random((POP, len(DV))) * (bounds_norm[:, 1] - bounds_norm[:, 0])

problem = AirplaneProblem()
algorithm = NSGA2(pop_size=POP, sampling=X0_pop, eliminate_duplicates=True)

t0 = time.time()
res = minimize(problem, algorithm, ('n_gen', NGEN), seed=1, verbose=False)
print('Rodada ingenua: %d pontos, %.0f s' % (len(res.F), time.time() - t0))

W0_naive = res.F[:, 0] * W0_ref
Wf_naive = res.F[:, 1] * Wf_ref

dados = np.loadtxt(RES + '/equipe_moga_frente.csv')
W0_final, Wf_final = dados[:, 0], dados[:, 1]

print('min W0: ingenua %.1f | final %.1f | ancora %.1f'
      % (W0_naive.min(), W0_final.min(), anchor_W0['W0']))
print('min Wf: ingenua %.1f | final %.1f | ancora %.1f'
      % (Wf_naive.min(), Wf_final.min(), anchor_Wf['Wf']))

fig, ax = plt.subplots(figsize=(8.5, 6))

o = np.argsort(W0_naive)
ax.plot(W0_naive[o] / 1000, Wf_naive[o] / 1000, 'o-', color=MUTED, markersize=7,
        linewidth=1.0, markeredgecolor='white', markeredgewidth=0.8,
        label='NSGA-II ingênuo (pop 60 × 120 ger., sem sementes)')
ax.annotate('estagnou dominado:\n%d ponto(s), %+.1f%% da âncora'
            % (len(W0_naive), 100 * (W0_naive.min() - anchor_W0['W0']) / anchor_W0['W0']),
            (W0_naive.min() / 1000, Wf_naive[np.argmin(W0_naive)] / 1000),
            xytext=(-12, -26), textcoords='offset points',
            fontsize=9, color=INK2, ha='right')

o = np.argsort(W0_final)
ax.plot(W0_final[o] / 1000, Wf_final[o] / 1000, 'o-', color=PAL[0], markersize=5,
        linewidth=1.0, markeredgecolor='white', markeredgewidth=0.8,
        label='NSGA-II final (pop 80 × 200 ger., semeado)')

for anc, lab, dxy in ((anchor_W0, 'SLSQP min $W_0$', (-10, -16)),
                      (anchor_Wf, 'SLSQP min $W_f$', (6, -16))):
    ax.plot(anc['W0'] / 1000, anc['Wf'] / 1000, '*', color=INK, markersize=15,
            markeredgecolor='white', markeredgewidth=0.8,
            label='ótimos SLSQP (âncoras)' if anc is anchor_W0 else None)
    ax.annotate(lab, (anc['W0'] / 1000, anc['Wf'] / 1000), xytext=dxy,
                textcoords='offset points', fontsize=9, color=INK2)

ax.set_xlabel('$W_0$ [t]', fontsize=13)
ax.set_ylabel('$W_f$ [t]', fontsize=13)
ax.margins(x=0.08, y=0.12)
ax.legend(fontsize=9, frameon=False, loc='lower left')
style_axes(ax)
ax.set_title('Convergência do NSGA-II verificada pelos ótimos mono-objetivo',
             fontsize=12, color=INK)

plt.tight_layout()
fig.savefig(RES + '/equipe_moga_convergencia.png', dpi=150)

plt.show()
