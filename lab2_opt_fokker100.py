import os

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.constants import gravity

Xlist = []
flist = []
g1list = []

neval = [0]

b_w_max = 30.0

# figuras organizadas por topico do roteiro
RES = 'Resultados/1_monoobj_fokker100'
os.makedirs(RES, exist_ok=True)

cache = {'x': None, 'out': None}

Xref = np.array([7.5, 90.0])

def run_analysis(x):

    if cache['x'] is not None and np.array_equal(x, cache['x']):
        return cache['out']

    X = np.asarray(x)*Xref

    AR_w = X[0]
    S_w = X[1]

    airplane = standard_airplane('fokker100')
    airplane['inputs']['AR_w'] = AR_w
    airplane['inputs']['S_w'] = S_w

    analyze(airplane)
    neval[0] = neval[0] + 1

    MTOW = airplane['thrust_matching']['W0']/gravity # [kgf]
    b_w = airplane['geometry']['b_w'] # [m]

    cache['x'] = np.array(x)
    cache['out'] = (MTOW, b_w)

    Xlist.append(X)
    flist.append(MTOW)
    g1list.append(b_w/b_w_max - 1)

    return MTOW, b_w

def objfun(x):

    MTOW, b_w = run_analysis(x)

    f = MTOW/MTOW_ref

    return f

def confun(x):

    MTOW, b_w = run_analysis(x)

    # >= 0
    g1 = 1 - b_w/b_w_max

    return g1

con1 = {'type': 'ineq',
        'fun': confun}

cons = [con1]

X0 = np.array([7.5, 90.0])
bounds_phys = np.array([[7.0, 12.0], [80.0, 120.0]])

x0 = X0/Xref
bounds = bounds_phys/Xref[:,None]

MTOW_ref, b_w_ref = run_analysis(x0)

Xlist.clear()
flist.clear()
g1list.clear()
neval[0] = 0
cache['x'] = None

options = {'maxiter': 200, 'ftol': 1e-8}

xk = [x0.copy()]

result = minimize(objfun, x0,
                  constraints = cons, bounds = bounds,
                  method = 'slsqp', options = options,
                  callback = lambda x: xk.append(x.copy()))

print(result)
xopt = result.x
Xopt = xopt*Xref

MTOW_opt, b_w_opt = run_analysis(xopt)

print('')
print('='*58)
print('RESULTADOS DA OTIMIZACAO')
print('='*58)
print('                 AR_w      S_w [m2]   MTOW [kgf]   b_w [m]')
print('Starting Point   %6.3f    %7.2f    %8.1f    %6.2f'%(X0[0], X0[1], MTOW_ref, b_w_ref))
print('Optimized Point  %6.3f    %7.2f    %8.1f    %6.2f'%(Xopt[0], Xopt[1], MTOW_opt, b_w_opt))
print('-'*58)
print('Melhoria relativa do objetivo: %.3f %%'%(100*(MTOW_ref - MTOW_opt)/MTOW_ref))
print('Chamadas da funcao de analise: %d'%neval[0])
print('Chamadas de objfun (nfev):     %d'%result.nfev)
print('Iteracoes do SLSQP (nit):      %d'%result.nit)
print('-'*58)
print('ATIVIDADE DAS RESTRICOES (tolerancia de 1e-4)')

tol_ativo = 1e-4

g1_opt = b_w_opt/b_w_max - 1
ativo_g1 = abs(g1_opt) < tol_ativo
ativo_AR = min(abs(xopt[0]-bounds[0][0]), abs(xopt[0]-bounds[0][1])) < tol_ativo
ativo_S = min(abs(xopt[1]-bounds[1][0]), abs(xopt[1]-bounds[1][1])) < tol_ativo

print('  g1 = b_w/30 - 1 = %+.4e  (b_w = %.3f m) -> %s'%(g1_opt, b_w_opt, 'ATIVA' if ativo_g1 else 'inativa'))
print('  bound AR_w em [%.0f, %.0f]    -> %s'%(bounds_phys[0][0], bounds_phys[0][1], 'ATIVO' if ativo_AR else 'inativo'))
print('  bound S_w  em [%.0f, %.0f]  -> %s'%(bounds_phys[1][0], bounds_phys[1][1], 'ATIVO' if ativo_S else 'inativo'))
print('')
print('  Norma do gradiente do objetivo no otimo: |jac| = %.3e'%np.linalg.norm(result.jac[:2]))

if not (ativo_g1 or ativo_AR or ativo_S):
    print('  => OTIMO IRRESTRITO (interior). Nenhuma restricao ativa:')
    print('     o gradiente do objetivo se anula, e a condicao de otimalidade')
    print('     e simplesmente grad(f) = 0, sem contribuicao dos multiplicadores.')
else:
    print('  => OTIMO RESTRINGIDO.')
print('='*58)

Xhist = np.array(Xlist)

# paleta e tintas (mesma convencao dos demais scripts do lab)
PAL = ['#2a78d6', '#eb6834']
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

fig, axs = plt.subplots(4, 1, figsize=(7,9), sharex=True)

# cada DV no seu painel, em unidades fisicas (eixo y duplo num painel so
# induz correlacao falsa entre escalas arbitrarias)
axs[0].plot(Xhist[:,0], '-', linewidth=1.6, color=PAL[0])
axs[0].set_ylabel(r'$AR_w$', fontsize=12)

axs[1].plot(Xhist[:,1], '-', linewidth=1.6, color=PAL[1])
axs[1].set_ylabel(r'$S_w$ [m$^2$]', fontsize=12)

axs[2].plot(flist, '-', linewidth=1.8, color=PAL[0])
axs[2].set_ylabel('MTOW [kgf]', fontsize=12)

axs[3].plot(g1list, '-', linewidth=1.6, color=PAL[0], label=r'$g_1 = b_w/30 - 1$')
axs[3].axhline(0, color=INK2, linewidth=0.8)
axs[3].set_ylabel('$g$', fontsize=12)
axs[3].set_xlabel('avaliações', fontsize=12)
axs[3].legend(fontsize=9, frameon=False)

for ax in axs:
    style_axes(ax)

plt.tight_layout()
fig.savefig(RES + '/fokker100_historico.png', dpi=150)

PLOT_CONTOUR = True

if PLOT_CONTOUR:

    ngrid = 60
    AR_grid = np.linspace(bounds_phys[0][0], bounds_phys[0][1], ngrid)
    S_grid = np.linspace(bounds_phys[1][0], bounds_phys[1][1], ngrid)
    AR_mesh, S_mesh = np.meshgrid(AR_grid, S_grid)
    MTOW_mesh = np.zeros_like(AR_mesh)
    bw_mesh = np.zeros_like(AR_mesh)

    print('')
    print('Gerando contorno do espaco de projeto (%d pontos)...'%(ngrid*ngrid))

    for ii in range(ngrid):
        for jj in range(ngrid):
            airplane = standard_airplane('fokker100')
            airplane['inputs']['AR_w'] = AR_mesh[ii,jj]
            airplane['inputs']['S_w'] = S_mesh[ii,jj]
            analyze(airplane)
            MTOW_mesh[ii,jj] = airplane['thrust_matching']['W0']/gravity
            bw_mesh[ii,jj] = airplane['geometry']['b_w']

    fig = plt.figure(figsize=(7.5,6))
    ax = plt.gca()

    pcm = ax.pcolormesh(AR_mesh, S_mesh, MTOW_mesh,
                        cmap='viridis', shading='gouraud', zorder=0)

    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label('MTOW [kgf]', fontsize=13)

    cs = ax.contour(AR_mesh, S_mesh, MTOW_mesh, levels=15,
                    colors='white', linewidths=0.5, alpha=0.5, zorder=1)
    ax.clabel(cs, inline=True, fontsize=6, fmt='%.0f', colors='white')

    ax.contour(AR_mesh, S_mesh, bw_mesh - b_w_max, levels=[0],
               colors='r', linewidths=2.5, zorder=3)

    xk_arr = np.array(xk)*Xref
    ax.plot(xk_arr[:,0], xk_arr[:,1], '--o', color='#e34948', markersize=4,
            markeredgecolor='w', markeredgewidth=0.5, zorder=4)
    ax.plot(X0[0], X0[1], 'o', color='lime', markersize=10,
            markeredgecolor='k', label='Ponto de partida', zorder=5)
    ax.plot(Xopt[0], Xopt[1], 'o', color='deepskyblue', markersize=10,
            markeredgecolor='k', label='Otimo', zorder=5)

    ax.plot([], [], '-r', linewidth=2.5, label=r'$b_w = 30$ m')

    ax.set_xlabel(r'$AR_w$', fontsize=15)
    ax.set_ylabel(r'$S_w$ [m$^2$]', fontsize=15)
    ax.legend(loc='lower left', framealpha=0.9)
    plt.tight_layout()
    fig.savefig(RES + '/fokker100_espaco_projeto.png', dpi=150)

plt.show()
