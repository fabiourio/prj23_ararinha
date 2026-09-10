import time

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.plots import plot_geometry
from designTool.constants import gravity

rad2deg = 180/np.pi

# nome, valor de referencia, limite inferior, limite superior
DV = [('AR_w',      8.0,   7.0,  12.0),
      ('xr_w',     17.0,  14.0,  20.0),
      ('S_w',     390.0, 330.0, 460.0),
      ('sweep_w',  0.58,  0.40,  0.70),
      ('x_mlg',    31.0,  28.0,  35.0),
      ('tcr_w',     0.18,  0.12,  0.24)]

dv_names = [d[0] for d in DV]
Xref = np.array([d[1] for d in DV])
bounds_phys = np.array([[d[2], d[3]] for d in DV])

# nome, sentido, limite. Todas convertidas para g >= 0 em constraints().
CON = [('deltaS_wlan',      '>=',  0.0),
       ('SM_fwd',           '<=',  0.30),
       ('SM_aft',           '>=',  0.05),
       ('frac_nlg_fwd',     '<=',  0.18),
       ('frac_nlg_aft',     '>=',  0.03),
       ('alpha_tipback',    '>=', 15.0),
       ('alpha_tailstrike', '>=', 10.0),
       ('phi_overturn',     '<=', 63.0),
       ('tank_excess',      '>=',  0.0),
       ('b_w',              '<=', 65.0),
       ('mlg_track',        '<=', 14.0),
       ('h_tail',           '<=', 20.0),
       ('T0req',            '<=',  1.0)]

con_names = [c[0] for c in CON]

Xlist = []
flist = []
glist = []

neval = [0]

cache = {'x': None, 'out': None}

def run_analysis(x):

    if cache['x'] is not None and np.array_equal(x, cache['x']):
        return cache['out']

    X = np.asarray(x)*Xref

    airplane = standard_airplane('my_airplane_v1')
    for name, value in zip(dv_names, X):
        airplane['inputs'][name] = value

    analyze(airplane)
    neval[0] = neval[0] + 1

    g = airplane['geometry']
    b = airplane['balance']
    lg = airplane['landing_gear']
    tm = airplane['thrust_matching']
    i = airplane['inputs']

    out = {'W0'               : tm['W0']/gravity,
           'Wf'               : tm['W_fuel']/gravity,
           # fracao da area da asa, e nao m2: o valor bruto vale ~130 e
           # dominaria a jacobiana, ja que o limite zero impede normalizar
           # pelo limite como nas demais restricoes
           'deltaS_wlan'      : tm['deltaS_wlan']/i['S_w'],
           'SM_fwd'           : b['SM_fwd'],
           'SM_aft'           : b['SM_aft'],
           'frac_nlg_fwd'     : lg['frac_nlg_fwd'],
           'frac_nlg_aft'     : lg['frac_nlg_aft'],
           'alpha_tipback'    : lg['alpha_tipback']*rad2deg,
           'alpha_tailstrike' : lg['alpha_tailstrike']*rad2deg,
           'phi_overturn'     : lg['phi_overturn']*rad2deg,
           'tank_excess'      : b['tank_excess'],
           'b_w'              : g['b_w'],
           'mlg_track'        : lg['mlg_track'],
           'h_tail'           : (i['zr_v'] + g['b_v']) - i['z_lg'],
           'T0req'            : max(tm['T0req'].values())/tm['T0']}

    cache['x'] = np.array(x)
    cache['out'] = out

    Xlist.append(X)
    flist.append(out['W0'])
    glist.append(constraints(x))

    return out

def objfun(x):

    out = run_analysis(x)

    f = out['W0']/W0_ref

    return f

def constraints(x):

    out = run_analysis(x)

    g = []
    for name, sense, lim in CON:
        v = out[name]
        # normalizacao pelo limite, e nao pelo valor corrente: tank_excess
        # vale 0.012 no baseline e dominaria a jacobiana se dividido por si mesmo
        if lim == 0.0:
            g.append(v if sense == '>=' else -v)
        elif sense == '>=':
            g.append(v/lim - 1)
        else:
            g.append(1 - v/lim)

    return np.array(g)

cons = [{'type': 'ineq', 'fun': constraints}]

X0 = Xref.copy()
x0 = X0/Xref
# DVs de referencia negativa (z_lg) invertem a ordem dos limites ao normalizar
bounds = np.sort(bounds_phys/Xref[:,None], axis=1)

out_ref = run_analysis(x0)
W0_ref = out_ref['W0']
Wf_ref = out_ref['Wf']
g_ref = constraints(x0)

Xlist.clear()
flist.clear()
glist.clear()
neval[0] = 0
cache['x'] = None

options = {'maxiter': 200, 'ftol': 1e-8}

xk = [x0.copy()]

t_start = time.time()
result = minimize(objfun, x0,
                  constraints = cons, bounds = bounds,
                  method = 'slsqp', options = options,
                  callback = lambda x: xk.append(x.copy()))
t_elapsed = time.time() - t_start

print(result)

xopt = result.x
Xopt = xopt*Xref
out_opt = run_analysis(xopt)
g_opt = constraints(xopt)

print('')
print('='*70)
print('VARIAVEIS DE PROJETO')
print('='*70)
print('%-18s %12s %12s %12s'%('', 'inicial', 'otimizado', 'variacao'))
for k, name in enumerate(dv_names):
    print('%-18s %12.4f %12.4f %11.2f%%'%(name, X0[k], Xopt[k],
                                          100*(Xopt[k]-X0[k])/X0[k]))

print('')
print('='*70)
print('OBJETIVO')
print('='*70)
print('%-18s %12.1f %12.1f %11.2f%%'%('W0 [kgf]', W0_ref, out_opt['W0'],
                                      100*(out_opt['W0']-W0_ref)/W0_ref))
print('%-18s %12.1f %12.1f %11.2f%%'%('Wf [kgf]', Wf_ref, out_opt['Wf'],
                                      100*(out_opt['Wf']-Wf_ref)/Wf_ref))

print('')
print('='*70)
print('RESTRICOES  (g normalizada >= 0; ATIVA se |g| < 1e-4)')
print('='*70)
print('%-18s %10s %10s %10s %10s  %s'%('', 'inicial', 'otimizado', 'limite', 'g_norm', 'estado'))

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
    print('%-18s %10.4f %10.4f %4s %5.2f %10.2e  %s'%(name, out_ref[name], out_opt[name],
                                                      sense, lim, g_opt[k], estado))

print('')
for k, name in enumerate(dv_names):
    for b in bounds[k]:
        if abs(xopt[k]-b) < tol_ativo:
            print('bound ATIVO: %s = %.4f'%(name, b*Xref[k]))
            n_ativas += 1

print('')
print('='*70)
print('DESEMPENHO')
print('='*70)
print('Iteracoes do SLSQP (nit):      %d'%result.nit)
print('Chamadas de objfun (nfev):     %d'%result.nfev)
print('Execucoes do analyze:          %d'%neval[0])
print('Tempo de otimizacao:           %.2f s'%t_elapsed)
print('Restricoes ativas:             %d'%n_ativas)
print('='*70)

Xhist = np.array(Xlist)
ghist = np.array(glist)

fig, axs = plt.subplots(3, 1, figsize=(8,9), sharex=True)

for k, name in enumerate(dv_names):
    axs[0].plot(Xhist[:,k]/Xref[k], 'o-', markersize=3, label=name)
axs[0].set_ylabel('DV / valor inicial', fontsize=12)
axs[0].legend(fontsize=9, ncol=3)

axs[1].plot(flist, 'o-', markersize=3, color='k')
axs[1].set_ylabel('$W_0$ [kgf]', fontsize=13)

for k, name in enumerate(con_names):
    axs[2].plot(ghist[:,k], 'o-', markersize=2.5, linewidth=1, label=name)
axs[2].axhline(0, color='gray', linewidth=0.8)
axs[2].set_ylabel('$g$ normalizada', fontsize=13)
axs[2].set_xlabel('avaliações', fontsize=13)
axs[2].set_ylim(-0.2, 3.0)
axs[2].legend(fontsize=7, ncol=3, loc='upper right')

plt.tight_layout()
fig.savefig('Resultados/2_monoobj_equipe/equipe_etapa2_historico.png', dpi=150)

airplane_opt = standard_airplane('my_airplane_v1')
for name, value in zip(dv_names, Xopt):
    airplane_opt['inputs'][name] = value
analyze(airplane_opt)

airplane_ref = standard_airplane('my_airplane_v1')
analyze(airplane_ref)

def planform(ax, ap, color, label):

    i = ap['inputs']
    g = ap['geometry']

    for s in (1, -1):
        ax.plot([i['xr_w'], g['xt_w'], g['xt_w']+g['ct_w'], i['xr_w']+g['cr_w'], i['xr_w']],
                s*np.array([0, g['yt_w'], g['yt_w'], 0, 0]), color=color, lw=1.5,
                label=label if s == 1 else None)
        ax.plot([g['xr_h'], g['xt_h'], g['xt_h']+g['ct_h'], g['xr_h']+g['cr_h'], g['xr_h']],
                s*np.array([0, g['yt_h'], g['yt_h'], 0, 0]), color=color, lw=1.5)
        ax.plot([0, i['L_f']], [s*i['D_f']/2, s*i['D_f']/2], color=color, lw=1.0)

    ax.plot([i['x_mlg']]*2, [-i['y_mlg'], i['y_mlg']], 'o', color=color, ms=4)
    ax.plot([i['x_nlg']], [0], 'o', color=color, ms=4)

def sideview(ax, ap, color, label):

    i = ap['inputs']
    g = ap['geometry']

    th = np.linspace(0, 2*np.pi, 120)
    ax.plot(i['L_f']/2 + i['L_f']/2*np.cos(th), i['D_f']/2*np.sin(th),
            color=color, lw=1.0, label=label)

    ax.plot([g['xr_v'], g['xt_v'], g['xt_v']+g['ct_v'], g['xr_v']+g['cr_v'], g['xr_v']],
            [i['zr_v'], g['zt_v'], g['zt_v'], i['zr_v'], i['zr_v']], color=color, lw=1.5)
    ax.plot([i['xr_w'], i['xr_w']+g['cr_w']], [i['zr_w'], i['zr_w']], color=color, lw=2.5)
    ax.plot([g['xr_h'], g['xr_h']+g['cr_h']], [i['zr_h'], i['zr_h']], color=color, lw=2.0)

    ax.plot([i['x_mlg'], i['x_nlg']], [i['z_lg'], i['z_lg']], 'o', color=color, ms=4)
    ax.plot([i['x_tailstrike']], [i['z_tailstrike']], 's', color=color, ms=4)
    ax.axhline(i['z_lg'], color=color, lw=0.6, ls=':')

fig, axs = plt.subplots(2, 1, figsize=(11,9))

planform(axs[0], airplane_ref, 'tab:gray', 'baseline')
planform(axs[0], airplane_opt, 'tab:red', 'otimizado')
axs[0].set_title('Planta', fontsize=13)
axs[0].set_ylabel('y [m]', fontsize=12)
axs[0].legend()

sideview(axs[1], airplane_ref, 'tab:gray', 'baseline')
sideview(axs[1], airplane_opt, 'tab:red', 'otimizado')
axs[1].set_title('Vista lateral', fontsize=13)
axs[1].set_xlabel('x [m]', fontsize=12)
axs[1].set_ylabel('z [m]', fontsize=12)

for ax in axs:
    ax.set_aspect('equal')
    ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig('Resultados/2_monoobj_equipe/equipe_etapa2_planformas.png', dpi=150)

plt.show()
