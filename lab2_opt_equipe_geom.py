import os
import time
import pprint

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.plots import plot_geometry
from designTool.constants import gravity

rad2deg = 180/np.pi

# figuras organizadas por topico do roteiro
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
Xref = np.array([d[1] for d in DV])
bounds_phys = np.array([[d[2], d[3]] for d in DV])

CON = [('deltaS_wlan',      '>=',  0.0),
       ('SM_fwd',           '<=',  0.30),
       ('SM_aft',           '>=',  0.05),
       ('frac_nlg_fwd',     '<=',  0.18),
       ('frac_nlg_aft',     '>=',  0.03),
       ('alpha_tipback',    '>=', 15.0),
       ('alpha_tailstrike', '>=', 10.0),
       ('phi_overturn',     '<=', 63.0),
       # fundo da nacele (z_n - D_n/2) nao pode raspar no solo (z = z_lg)
       ('ground_clearance', '>=',  0.5),
       ('tank_excess',      '>=',  0.0),
       # a Tab. 1-1 do Anexo 14 define a categoria E como "up to but not
       # including" 65 m e 14 m, entao os limites sao estritos
       ('b_w',              '<=', 64.9),
       ('mlg_track',        '<=', 13.9),
       ('h_tail',           '<=', 20.0),
       ('T0req',            '<=',  1.0),
       ('vt_fit',           '<=',  1.0),
       ('mlg_fit',          '<=',  1.0)]

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

    # corda e bordo de fuga da asa na estacao lateral do trem principal
    # (LE interpolado entre raiz e ponta: xt_w ja carrega o termo de 1/4 de
    # corda, pois sweep_w e' medido a 1/4 de corda em geometry.py)
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
           # raiz da empenagem vertical deve terminar dentro da fuselagem
           'vt_fit'           : (g['xr_v'] + g['cr_v'])/i['L_f'],
           # trem principal deve ficar sob a asa, nao atras dela
           'mlg_fit'          : i['x_mlg']/te_mlg}

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

airplane_opt = standard_airplane('my_airplane_v1')
for name, value in zip(dv_names, Xopt):
    airplane_opt['inputs'][name] = value
analyze(airplane_opt)

airplane_ref = standard_airplane('my_airplane_v1')
analyze(airplane_ref)

Xhist = np.array(Xlist)
ghist = np.array(glist)

# paleta categorica validada para daltonismo (ordem fixa) + tintas de apoio
PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100',
       '#e87ba4', '#008300', '#4a3aa7', '#e34948']
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

fig, axs = plt.subplots(3, 1, figsize=(9,10), sharex=True)

n_ev = len(flist)

# DVs: rotulos diretos so nos maiores movimentos (o resto fica na legenda)
destaque_dv = ('y_mlg', 'AR_w', 'tcr_w', 'S_w')
for k, name in enumerate(dv_names):
    axs[0].plot(Xhist[:,k]/Xref[k], '-', linewidth=1.6, color=PAL[k], label=name)
    if name in destaque_dv:
        axs[0].annotate(name, (n_ev-1, Xhist[-1,k]/Xref[k]),
                        xytext=(5, 0), textcoords='offset points',
                        fontsize=8, color=INK2, va='center')
axs[0].set_ylabel('DV / valor inicial', fontsize=12)
axs[0].legend(fontsize=8, ncol=4, frameon=False)
axs[0].set_xlim(-2, n_ev*1.09)

axs[1].plot(flist, '-', linewidth=1.8, color=PAL[0])
axs[1].annotate('%.0f kgf\n(%.2f%%)'%(flist[-1], 100*(flist[-1]-W0_ref)/W0_ref),
                (n_ev-1, flist[-1]), xytext=(5, 0), textcoords='offset points',
                fontsize=8, color=INK2, va='center')
axs[1].set_ylabel('$W_0$ [kgf]', fontsize=13)

# restricoes: 16 series estoura o limite de cores com significado, entao
# so as ativas no otimo ganham cor; as demais viram contexto em cinza
ativas_idx = [k for k in range(len(CON)) if abs(g_opt[k]) < tol_ativo]
ic = 0
tem_label_inativa = False
for k, name in enumerate(con_names):
    if k in ativas_idx:
        axs[2].plot(ghist[:,k], '-', linewidth=1.6, color=PAL[ic], label=name)
        ic += 1
    else:
        axs[2].plot(ghist[:,k], '-', linewidth=0.9, color=MUTED,
                    label=None if tem_label_inativa else 'inativas')
        tem_label_inativa = True
axs[2].axhline(0, color=INK2, linewidth=0.8)
axs[2].set_ylabel('$g$ normalizada', fontsize=13)
axs[2].set_xlabel('avaliações', fontsize=13)
axs[2].set_ylim(-0.2, 2.0)
axs[2].legend(fontsize=8, ncol=4, frameon=False, loc='upper right')

for ax in axs:
    style_axes(ax)

plt.tight_layout()
fig.savefig(RES + '/equipe_geom_historico.png', dpi=150)

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

    ax.plot([i['x_mlg']]*2, [-i['y_mlg'], i['y_mlg']], 'o', color=color, ms=5)
    ax.plot([i['x_nlg']], [0], 'o', color=color, ms=5)

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

    ax.plot([i['x_mlg'], i['x_nlg']], [i['z_lg'], i['z_lg']], 'o', color=color, ms=5)
    ax.plot([i['x_tailstrike']], [i['z_tailstrike']], 's', color=color, ms=5)
    ax.axhline(i['z_lg'], color=color, lw=0.6, ls=':')

fig, axs = plt.subplots(2, 1, figsize=(11,9))

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
print('='*70)
print('DICIONARIO DA AERONAVE OTIMIZADA')
print('='*70)
print('inputs = ' + pprint.pformat(airplane_opt['inputs'], sort_dicts=False, width=100))

plot_geometry(airplane_opt, figname=RES + '/equipe_geom_3dview.png')

plt.show()
