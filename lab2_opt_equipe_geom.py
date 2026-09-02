'''
Lab 02 secao 3: min W0 da aeronave da equipe (8 DVs, 16 restricoes), SLSQP.

Saida: tabelas no terminal + Resultados/2_monoobj_equipe/*.png
Rodar: python lab2_opt_equipe_geom.py
'''

import os
import time
import pprint

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.plots import plot_geometry

from lab2_equipe_common import (DV, Xref, bounds_norm, con_names, CON,
                                dv_names, airplane_out, constraints_geq0)
from lab2_plot import PAL, INK2, style_axes, planform, sideview

RES = 'Resultados/2_monoobj_equipe'
os.makedirs(RES, exist_ok=True)

Xlist = []
flist = []
glist = []
neval = [0]
cache = {'x': None, 'out': None}


def run_analysis(x):
    '''Envolve airplane_out com cache e historico para os graficos.'''

    if cache['x'] is not None and np.array_equal(x, cache['x']):
        return cache['out']

    out = airplane_out(x)
    neval[0] += 1

    cache['x'] = np.array(x)
    cache['out'] = out

    Xlist.append(np.asarray(x) * Xref)
    flist.append(out['W0'])
    glist.append(constraints_geq0(out))

    return out


def objfun(x):
    return run_analysis(x)['W0'] / W0_ref


def constraints(x):
    return constraints_geq0(run_analysis(x))


cons = [{'type': 'ineq', 'fun': constraints}]

X0 = Xref.copy()
x0 = X0 / Xref

out_ref = run_analysis(x0)
W0_ref = out_ref['W0']
Wf_ref = out_ref['Wf']

# descarta avaliacoes do ponto de partida antes de contar desempenho
Xlist.clear()
flist.clear()
glist.clear()
neval[0] = 0
cache['x'] = None

options = {'maxiter': 200, 'ftol': 1e-8}
xk = [x0.copy()]

t_start = time.time()
result = minimize(objfun, x0,
                  constraints=cons, bounds=bounds_norm,
                  method='slsqp', options=options,
                  callback=lambda x: xk.append(x.copy()))
t_elapsed = time.time() - t_start

print(result)

xopt = result.x
Xopt = xopt * Xref
out_opt = run_analysis(xopt)
g_opt = constraints(xopt)

print('')
print('=' * 70)
print('VARIAVEIS DE PROJETO')
print('=' * 70)
print('%-18s %12s %12s %12s' % ('', 'inicial', 'otimizado', 'variacao'))
for k, name in enumerate(dv_names):
    print('%-18s %12.4f %12.4f %11.2f%%' % (name, X0[k], Xopt[k],
                                              100 * (Xopt[k] - X0[k]) / X0[k]))

print('')
print('=' * 70)
print('OBJETIVO')
print('=' * 70)
print('%-18s %12.1f %12.1f %11.2f%%' % ('W0 [kgf]', W0_ref, out_opt['W0'],
                                          100 * (out_opt['W0'] - W0_ref) / W0_ref))
print('%-18s %12.1f %12.1f %11.2f%%' % ('Wf [kgf]', Wf_ref, out_opt['Wf'],
                                          100 * (out_opt['Wf'] - Wf_ref) / Wf_ref))

print('')
print('=' * 70)
print('RESTRICOES  (g normalizada >= 0; ATIVA se |g| < 1e-4)')
print('=' * 70)
print('%-18s %10s %10s %10s %10s  %s' % ('', 'inicial', 'otimizado', 'limite', 'g_norm', 'estado'))

tol_ativo = 1e-4
n_ativas = 0
for k, (name, sense, lim) in enumerate(CON):
    if abs(g_opt[k]) < tol_ativo:
        estado = 'ATIVA'
        n_ativas += 1
    elif g_opt[k] < 0:
        estado = 'VIOLADA'
    else:
        estado = ''
    print('%-18s %10.4f %10.4f %4s %5.2f %10.2e  %s' % (name, out_ref[name], out_opt[name],
                                                          sense, lim, g_opt[k], estado))

print('')
for k, name in enumerate(dv_names):
    for b in bounds_norm[k]:
        if abs(xopt[k] - b) < tol_ativo:
            print('bound ATIVO: %s = %.4f' % (name, b * Xref[k]))
            n_ativas += 1

print('')
print('=' * 70)
print('DESEMPENHO')
print('=' * 70)
print('Iteracoes do SLSQP (nit):      %d' % result.nit)
print('Chamadas de objfun (nfev):     %d' % result.nfev)
print('Execucoes do analyze:          %d' % neval[0])
print('Tempo de otimizacao:           %.2f s' % t_elapsed)
print('Restricoes ativas:             %d' % n_ativas)
print('=' * 70)

airplane_opt = standard_airplane('my_airplane')
for name, value in zip(dv_names, Xopt):
    airplane_opt['inputs'][name] = value
analyze(airplane_opt)

airplane_ref = standard_airplane('my_airplane')
analyze(airplane_ref)

Xhist = np.array(Xlist)
ghist = np.array(glist)

fig, axs = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
n_ev = len(flist)

for k, name in enumerate(dv_names):
    axs[0].plot(Xhist[:, k] / Xref[k], '-', linewidth=1.6, color=PAL[k], label=name)
axs[0].set_ylabel('DV / valor inicial', fontsize=12)
axs[0].legend(fontsize=8, ncol=4, frameon=False)
axs[0].set_xlim(-2, n_ev * 1.09)

axs[1].plot(flist, '-', linewidth=1.8, color=PAL[0])
axs[1].annotate('%.0f kgf\n(%.2f%%)' % (flist[-1], 100 * (flist[-1] - W0_ref) / W0_ref),
                (n_ev - 1, flist[-1]), xytext=(5, 0), textcoords='offset points',
                fontsize=8, color=INK2, va='center')
axs[1].set_ylabel('$W_0$ [kgf]', fontsize=13)

# 16 curvas juntas ficam ilegiveis: colorir so as ativas no otimo
ativas_idx = [k for k in range(len(CON)) if abs(g_opt[k]) < tol_ativo]
ic = 0
tem_inativa = False
for k, name in enumerate(con_names):
    if k in ativas_idx:
        axs[2].plot(ghist[:, k], '-', linewidth=1.6, color=PAL[ic], label=name)
        ic += 1
    else:
        axs[2].plot(ghist[:, k], '-', linewidth=0.9, color='#c3c2b7',
                    label=None if tem_inativa else 'inativas')
        tem_inativa = True
axs[2].axhline(0, color=INK2, linewidth=0.8)
axs[2].set_ylabel('$g$ normalizada', fontsize=13)
axs[2].set_xlabel('avaliações', fontsize=13)
axs[2].set_ylim(-0.2, 2.0)
axs[2].legend(fontsize=8, ncol=4, frameon=False, loc='upper right')

for ax in axs:
    style_axes(ax)

plt.tight_layout()
fig.savefig(RES + '/equipe_geom_historico.png', dpi=150)

fig, axs = plt.subplots(2, 1, figsize=(11, 9))

planform(axs[0], airplane_ref, '#898781', 'baseline')
planform(axs[0], airplane_opt, '#e34948', 'otimizado')
axs[0].set_title('Planta', fontsize=13)
axs[0].set_ylabel('y [m]', fontsize=12)
axs[0].legend()

sideview(axs[1], airplane_ref, '#898781', 'baseline')
sideview(axs[1], airplane_opt, '#e34948', 'otimizado')
axs[1].set_title('Vista lateral', fontsize=13)
axs[1].set_xlabel('x [m]', fontsize=12)
axs[1].set_ylabel('z [m]', fontsize=12)

for ax in axs:
    ax.set_aspect('equal')
    style_axes(ax)

plt.tight_layout()
fig.savefig(RES + '/equipe_geom_planformas.png', dpi=150)

print('')
print('=' * 70)
print('DICIONARIO DA AERONAVE OTIMIZADA')
print('=' * 70)
print('inputs = ' + pprint.pformat(airplane_opt['inputs'], sort_dicts=False, width=100))

plot_geometry(airplane_opt, figname=RES + '/equipe_geom_3dview.png')

plt.show()
