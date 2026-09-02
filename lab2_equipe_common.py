'''
Definicao unica do problema da equipe (Lab 02, secoes 3 e 4).

Todos os scripts de otimizacao importam daqui para nao divergir em DVs,
restrições ou normalizacao. A analise fisica continua no designTool.

Convencoes
----------
DV : tuplas (nome, valor_ref, lim_inf, lim_sup). O otimizador trabalha com
     x = X / Xref; bounds_norm vem de bounds_phys/Xref com sort por linha
     (necessario porque z_lg < 0 inverte lo/hi apos dividir por Xref).

CON : tuplas (nome, sentido, limite). constraints_geq0() devolve g com
      g >= 0 viavel. Formato do roteiro: (v/lim - 1) para >= e (1 - v/lim)
      para <=. Limite zero usa v ou -v direto (deltaS_wlan, tank_excess).

Sinais : SciPy minimize usa g >= 0. pymoo usa G <= 0, entao negamos g ao
         passar para o NSGA-II.

Referencias fixas
-----------------
anchor_W0 / anchor_Wf : otimos SLSQP usados para checar convergencia da
                        frente de Pareto (pergunta 3 da secao 4).
Xopt : mesmo ponto que anchor_W0, em unidades fisicas (DOE, figuras).
'''

import numpy as np

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.constants import gravity

rad2deg = 180 / np.pi

# nome, Xref, lim_inf, lim_sup  (baseline PRJ-22 = Xref)
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
bounds_norm = np.sort(bounds_phys / Xref[:, None], axis=1)

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
       ('b_w',              '<=', 64.9),   # Anexo 14 cat. E: strictly < 65 m
       ('mlg_track',        '<=', 13.9),   # idem, strictly < 14 m
       ('h_tail',           '<=', 20.0),
       ('T0req',            '<=',  1.0),
       ('vt_fit',           '<=',  1.0),
       ('mlg_fit',          '<=',  1.0)]

con_names = [c[0] for c in CON]

anchor_W0 = {'W0': 290117.1, 'Wf': 111036.5,
             'x': np.array([9.8001, 16.0113, 353.6712, 0.6108,
                            29.2629, 0.2093, -6.0247, 6.9500]) / Xref}
anchor_Wf = {'W0': 291961.0, 'Wf': 109117.9,
             'x': np.array([9.8416, 15.8699, 374.4731, 0.6066,
                            29.3427, 0.1872, -6.0106, 6.9500]) / Xref}

Xopt = anchor_W0['x'] * Xref


def airplane_out(x_norm):
    '''
    Avalia um design da equipe.

    Parametros
    ----------
    x_norm : array (8,)  DVs divididas por Xref (partida = ones(8)).

    Retorna
    -------
    dict com W0, Wf em kgf e demais metricas ja nas unidades das CON.
    '''

    X = np.asarray(x_norm) * Xref

    airplane = standard_airplane('my_airplane')
    for name, value in zip(dv_names, X):
        airplane['inputs'][name] = value

    analyze(airplane)

    g = airplane['geometry']
    b = airplane['balance']
    lg = airplane['landing_gear']
    tm = airplane['thrust_matching']
    i = airplane['inputs']

    # TE da asa na estacao y_mlg (sweep_w medido a 1/4 de corda no designTool)
    eta = i['y_mlg'] / (g['b_w'] / 2)
    c_mlg = g['cr_w'] - (g['cr_w'] - g['ct_w']) * eta
    te_mlg = i['xr_w'] + (g['xt_w'] - i['xr_w']) * eta + c_mlg

    return {
        'W0': tm['W0'] / gravity,
        'Wf': tm['W_fuel'] / gravity,
        'W_empty': tm['W_empty'] / gravity,
        # fracao da area de pouso sobre S_w (nao m2); evita escala ~130 na jacobiana
        'deltaS_wlan': tm['deltaS_wlan'] / i['S_w'],
        'SM_fwd': b['SM_fwd'],
        'SM_aft': b['SM_aft'],
        'frac_nlg_fwd': lg['frac_nlg_fwd'],
        'frac_nlg_aft': lg['frac_nlg_aft'],
        'alpha_tipback': lg['alpha_tipback'] * rad2deg,
        'alpha_tailstrike': lg['alpha_tailstrike'] * rad2deg,
        'phi_overturn': lg['phi_overturn'] * rad2deg,
        'ground_clearance': lg['ground_clearance'],
        'tank_excess': b['tank_excess'],
        'b_w': g['b_w'],
        'mlg_track': lg['mlg_track'],
        'h_tail': (i['zr_v'] + g['b_v']) - i['z_lg'],
        'T0req': max(tm['T0req'].values()) / tm['T0'],
        'vt_fit': (g['xr_v'] + g['cr_v']) / i['L_f'],   # raiz da EV dentro da fuselagem
        'mlg_fit': i['x_mlg'] / te_mlg,                 # MLG sob a asa, nao atras do TE
    }


def constraints_geq0(out):
    '''Restricoes no formato SciPy: viavel quando todas >= 0.'''

    gg = []
    for name, sense, lim in CON:
        v = out[name]
        if lim == 0.0:
            # tank_excess ~0.012 no baseline: nao dividir pelo valor corrente
            gg.append(v if sense == '>=' else -v)
        elif sense == '>=':
            gg.append(v / lim - 1)
        else:
            gg.append(1 - v / lim)
    return np.array(gg)
