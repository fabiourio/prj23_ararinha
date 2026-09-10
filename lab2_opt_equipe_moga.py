'''
Lab 02 - Secao 4: Otimizacao multiobjetivo da aeronave da equipe
Minimizar W0 e Wf simultaneamente com NSGA-II (pymoo), sujeito as mesmas
restricoes da otimizacao mono-objetivo (lab2_opt_equipe_geom.py).

A convergencia da frente de Pareto e verificada comparando o extremo de
minimo W0 com o otimo do SLSQP da Secao 3.
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

# a populacao aleatoria do NSGA-II gera geometrias extremas em que o ponto
# fixo do W0 diverge (overflow/NaN); os avisos sao esperados e inofensivos
np.seterr(all='ignore')
warnings.filterwarnings('ignore')

# figuras e dados organizados por topico do roteiro
RES = 'Resultados/3_multiobj'
os.makedirs(RES, exist_ok=True)

# mesmas DVs e limites da otimizacao mono-objetivo
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
       ('ground_clearance', '>=',  0.5),
       ('tank_excess',      '>=',  0.0),
       ('b_w',              '<=', 64.9),
       ('mlg_track',        '<=', 13.9),
       ('h_tail',           '<=', 20.0),
       ('T0req',            '<=',  1.0),
       ('vt_fit',           '<=',  1.0),
       ('mlg_fit',          '<=',  1.0)]

con_names = [c[0] for c in CON]

neval = [0]

def run_analysis(x):

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
           'vt_fit'           : (g['xr_v'] + g['cr_v'])/i['L_f'],
           'mlg_fit'          : i['x_mlg']/te_mlg}

    return out

def constraints_geq0(out):
    # mesma normalizacao do script mono-objetivo: viavel se g >= 0
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

# referencia (aeronave do PRJ-22) para normalizar os objetivos
out_ref = run_analysis(np.ones(len(DV)))
W0_ref = out_ref['W0']
Wf_ref = out_ref['Wf']

# otimos mono-objetivo (SLSQP) do mesmo problema, usados como ancoras para
# verificar a convergencia dos extremos da frente: min W0 vem do
# lab2_opt_equipe_geom.py; min Wf vem de uma rodada identica trocando o
# objetivo. Sem as sementes abaixo, o NSGA-II puro (pop 60, 120 geracoes)
# estagnou ~1% acima dessas ancoras -- evidencia de nao-convergencia.
anchor_W0 = {'W0': 290117.1, 'Wf': 111036.5,
             'x': np.array([9.8001, 16.0113, 353.6712, 0.6108,
                            29.2629, 0.2093, -6.0247, 6.9500])/Xref}
anchor_Wf = {'W0': 291961.0, 'Wf': 109117.9,
             'x': np.array([9.8416, 15.8699, 374.4731, 0.6066,
                            29.3427, 0.1872, -6.0106, 6.9500])/Xref}

bounds_norm = np.sort(bounds_phys/Xref[:,None], axis=1)

class AirplaneProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(n_var=len(DV),
                         n_obj=2,
                         n_ieq_constr=len(CON),
                         xl=bounds_norm[:,0],
                         xu=bounds_norm[:,1])

    def _evaluate(self, x, out_pymoo, *args, **kwargs):

        # designs extremos podem divergir a analise; devolvemos um
        # individuo fortemente inviavel para o NSGA-II descarta-lo
        try:
            out = run_analysis(x)
            F = [out['W0']/W0_ref, out['Wf']/Wf_ref]
            # ATENCAO: pymoo considera viavel G <= 0, o oposto do scipy
            G = list(-constraints_geq0(out))
            if not (np.all(np.isfinite(F)) and np.all(np.isfinite(G))):
                raise ValueError('analise divergiu')
        except Exception:
            F = [10.0, 10.0]
            G = [1e3]*len(CON)

        out_pymoo['F'] = F
        out_pymoo['G'] = G

POP_SIZE = 80
N_GEN = 200

problem = AirplaneProblem()

# sementes estritamente viaveis: otimos SLSQP do mesmo problema com as
# restricoes apertadas em `margem` -- ficam no interior do conjunto viavel,
# entao sobrevivem a reavaliacao do NSGA-II (as ancoras exatas ficam em cima
# de restricoes ativas e sao descartadas por violacao de arredondamento)
def slsqp_interior(objetivo, margem):
    from scipy.optimize import minimize as sp_minimize
    ref = W0_ref if objetivo == 'W0' else Wf_ref
    r = sp_minimize(lambda x: run_analysis(x)[objetivo]/ref,
                    np.ones(len(DV)),
                    constraints=[{'type': 'ineq',
                                  'fun': lambda x: constraints_geq0(run_analysis(x)) - margem}],
                    bounds=bounds_norm, method='slsqp',
                    options={'maxiter': 200, 'ftol': 1e-8})
    return r.x

seed_W0 = slsqp_interior('W0', 0.005)
seed_Wf = slsqp_interior('Wf', 0.005)

# populacao inicial: LHS semeado manualmente (o LHS do pymoo nao usa o RNG
# global do numpy e quebraria a reprodutibilidade) + baseline + sementes
rng = np.random.default_rng(1)
nv = len(DV)
u = (np.stack([rng.permutation(POP_SIZE) for _ in range(nv)], axis=1)
     + rng.random((POP_SIZE, nv)))/POP_SIZE
X0_pop = bounds_norm[:,0] + u*(bounds_norm[:,1] - bounds_norm[:,0])
X0_pop[0] = np.ones(nv)
X0_pop[1] = seed_W0
X0_pop[2] = seed_Wf

algorithm = NSGA2(pop_size=POP_SIZE, sampling=X0_pop, eliminate_duplicates=True)

neval[0] = 0
t_start = time.time()
res = minimize(problem, algorithm, ('n_gen', N_GEN), seed=1, verbose=True)
t_elapsed = time.time() - t_start

# ordena a frente por W0
order = np.argsort(res.F[:,0])
F = res.F[order]
X = res.X[order]

W0_front = F[:,0]*W0_ref
Wf_front = F[:,1]*Wf_ref

print('')
print('='*70)
print('OTIMIZACAO MULTIOBJETIVO (NSGA-II)')
print('='*70)
print('Individuos por geracao (pop_size): %d'%POP_SIZE)
print('Numero de geracoes (n_gen):        %d'%N_GEN)
print('Execucoes do analyze:              %d'%neval[0])
print('Tempo de otimizacao:               %.1f s'%t_elapsed)
print('Pontos na frente de Pareto:        %d'%len(F))
print('')
print('Verificacao de convergencia com os otimos mono-objetivo (SLSQP):')
print('  min W0 da frente:  %10.1f kgf | ancora SLSQP: %10.1f kgf (%+.2f %%)'
      %(W0_front[0], anchor_W0['W0'], 100*(W0_front[0]-anchor_W0['W0'])/anchor_W0['W0']))
print('  min Wf da frente:  %10.1f kgf | ancora SLSQP: %10.1f kgf (%+.2f %%)'
      %(Wf_front[-1], anchor_Wf['Wf'], 100*(Wf_front[-1]-anchor_Wf['Wf'])/anchor_Wf['Wf']))

# selecao de 3 aeronaves de regioes distintas da frente:
# A = extremo de minimo W0, C = extremo de minimo Wf,
# B = "joelho" (ponto mais proximo da utopia na frente normalizada)
f1n = (F[:,0] - F[:,0].min())/max(F[:,0].max() - F[:,0].min(), 1e-12)
f2n = (F[:,1] - F[:,1].min())/max(F[:,1].max() - F[:,1].min(), 1e-12)
idxA = 0
idxC = len(F) - 1
idxB = int(np.argmin(f1n**2 + f2n**2))

sel_idx = [idxA, idxB, idxC]
sel_names = ['A (min W0)', 'B (joelho)', 'C (min Wf)']

# paleta categorica validada para daltonismo (3 primeiros slots valem para
# scatter) + tintas de apoio, mesma convencao do lab2_opt_equipe_geom
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

# frente salva em CSV para replotar sem precisar re-rodar o NSGA-II
np.savetxt(RES + '/equipe_moga_frente.csv',
           np.column_stack([W0_front, Wf_front, X*Xref]),
           header='W0_kgf Wf_kgf ' + ' '.join(dv_names), fmt='%.6g')

print('')
print('%-12s %10s %10s'%('aeronave', 'W0 [kgf]', 'Wf [kgf]'))
for name, k in zip(sel_names, sel_idx):
    print('%-12s %10.1f %10.1f'%(name, W0_front[k], Wf_front[k]))

print('')
print('%-10s'%'DV' + ''.join('%12s'%n for n in sel_names))
for j, dvn in enumerate(dv_names):
    print('%-10s'%dvn + ''.join('%12.4f'%(X[k,j]*Xref[j]) for k in sel_idx))

### FRENTE DE PARETO
# eixo principal ampliado na regiao da frente (a baseline, 14 t acima,
# esmagaria os pontos num canto); contexto completo fica no inset

fig, ax = plt.subplots(figsize=(8.5, 6))

ax.plot(W0_front/1000, Wf_front/1000, '-', color=PAL[0], linewidth=1.1, alpha=0.55)
ax.plot(W0_front/1000, Wf_front/1000, 'o', color=PAL[0], markersize=5,
        markeredgecolor='white', markeredgewidth=1.0,
        label='frente de Pareto (NSGA-II)')

for anc, lab, dxy in ((anchor_W0, 'SLSQP min $W_0$', (-10, -16)),
                      (anchor_Wf, 'SLSQP min $W_f$', (10, -4))):
    ax.plot(anc['W0']/1000, anc['Wf']/1000, '*', color=INK, markersize=15,
            markeredgecolor='white', markeredgewidth=0.8,
            label='ótimos SLSQP (âncoras)' if anc is anchor_W0 else None)
    ax.annotate(lab, (anc['W0']/1000, anc['Wf']/1000), xytext=dxy,
                textcoords='offset points', fontsize=9, color=INK2)

sel_marks = ['o', 's', '^']
for name, k, m in zip(sel_names, sel_idx, sel_marks):
    ax.plot(W0_front[k]/1000, Wf_front[k]/1000, m, color=PAL[1], markersize=10,
            markeredgecolor='white', markeredgewidth=1.2, zorder=5,
            label='selecionadas (A, B, C)' if m == 'o' else None)
    ax.annotate(name.split()[0], (W0_front[k]/1000, Wf_front[k]/1000),
                xytext=(8, 7), textcoords='offset points',
                fontsize=10, color=INK, fontweight='bold')

ax.set_xlabel('$W_0$ [t]', fontsize=13)
ax.set_ylabel('$W_f$ [t]', fontsize=13)
ax.margins(x=0.10, y=0.12)
ax.legend(fontsize=9, frameon=False, loc='lower left')
style_axes(ax)

axi = ax.inset_axes([0.58, 0.58, 0.39, 0.38])
axi.plot(W0_front/1000, Wf_front/1000, 'o', color=PAL[0], markersize=2.5)
axi.plot(out_ref['W0']/1000, out_ref['Wf']/1000, 's', color=INK, markersize=6)
axi.annotate('baseline PRJ-22', (out_ref['W0']/1000, out_ref['Wf']/1000),
             xytext=(-8, -14), textcoords='offset points',
             fontsize=8, color=INK2, ha='right')
axi.annotate('frente', (W0_front[0]/1000, Wf_front[0]/1000),
             xytext=(6, 8), textcoords='offset points', fontsize=8, color=INK2)
axi.set_title('contexto completo', fontsize=8, color=INK2)
axi.tick_params(labelsize=7)
style_axes(axi)

plt.tight_layout()
fig.savefig(RES + '/equipe_moga_pareto.png', dpi=150)

### PLANFORMAS DAS AERONAVES SELECIONADAS

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

# as tres planformas sao quase identicas -- o trade-off real esta em S_w e
# t/c, entao um painel lateral mostra essas diferencas em escala legivel
fig = plt.figure(figsize=(12, 6.5))
gs = fig.add_gridspec(2, 2, width_ratios=[2.4, 1], hspace=0.4, wspace=0.15)
axp = fig.add_subplot(gs[:, 0])
axb1 = fig.add_subplot(gs[0, 1])
axb2 = fig.add_subplot(gs[1, 1])

for name, k, c in zip(sel_names, sel_idx, PAL):
    airplane = standard_airplane('my_airplane_v1')
    for j, dvn in enumerate(dv_names):
        airplane['inputs'][dvn] = X[k,j]*Xref[j]
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
Svals = [X[k, jS]*Xref[jS] for k in sel_idx]
Tvals = [X[k, jT]*Xref[jT] for k in sel_idx]

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
    pad = 0.35*(vmax - vmin) + 1e-9
    axb.set_xlim(vmin - pad, vmax + pad)
    axb.set_ylim(-0.6, 2.9)
    axb.set_title(titulo, fontsize=11, color=INK)
    style_axes(axb)

plt.tight_layout()
fig.savefig(RES + '/equipe_moga_planformas.png', dpi=150)

plt.show()
