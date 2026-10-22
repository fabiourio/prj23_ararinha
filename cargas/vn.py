'''
Envelope de voo (diagrama V-n) pelo 14 CFR Part 25, velocidades em EAS.

- n limite de manobra: 25.337(b) (positivo) e 25.337(c) (-1 ate V_C,
  variando linearmente ate 0 em V_D);
- V_S1 com CLmax e V_A = V_S1*sqrt(n_max): 25.335(c);
- V_C pelo Mach de cruzeiro na altitude considerada; V_D pela margem minima
  de 0,07 Mach de 25.335(b)(2) (altitude limitada por compressibilidade);
- V_B: 25.335(d) (curva de estol cruzando a rajada de V_C, limitada a V_C);
- rajada: U_ref de 25.341(a)(5) (56 ft/s ao nivel do mar, 44 ft/s a
  15000 ft, 20,86 ft/s a 60000 ft; metade em V_D), formula quase-estatica de
  Pratt com o fator de alivio K_g = 0,88 mu/(5,3 + mu).
'''

import numpy as np

from designTool.auxiliary import atmosphere
from designTool.constants import ft2m, gravity, kt2ms, lb2N

RHO0 = atmosphere(0.0)['density']
MARGEM_MD = 0.07


def n_max_manobra(W):
    '''25.337(b): 2,1 + 24000/(W + 10000), W em lb, entre 2,5 e 3,8.'''
    W_lb = W/lb2N
    return float(np.clip(2.1 + 24000.0/(W_lb + 10000.0), 2.5, 3.8))


def u_ref(h):
    '''Rajada de referencia em V_C [m/s EAS], 25.341(a)(5)(i).'''
    h_ft = h/ft2m
    if h_ft <= 15000.0:
        u = 56.0 + (44.0 - 56.0)*h_ft/15000.0
    else:
        u = 44.0 + (20.86 - 44.0)*(h_ft - 15000.0)/45000.0
    return u*ft2m


def razao_massa(W_S, rho, c, CLa):
    return 2.0*W_S/(rho*c*CLa*gravity)


def fator_alivio(mu):
    return 0.88*mu/(5.3 + mu)


def dn_rajada(V, U, W_S, CLa, Kg):
    '''Incremento de fator de carga da rajada (Pratt), V e U em EAS.'''
    return Kg*RHO0*U*V*CLa/(2.0*W_S)


def v_estol(W, S, CL):
    return np.sqrt(2.0*W/(RHO0*S*abs(CL)))


def envelope(W, S, c, h, CLmax, CLmin, CLa, M_C, margem_MD=MARGEM_MD):
    '''
    Velocidades caracteristicas [m/s EAS] e fatores de carga do V-n para o
    peso W [N] na altitude h [m]. CLa [1/rad] da aeronave.
    '''
    atm = atmosphere(h)
    rho, a = atm['density'], atm['speed_of_sound']
    raiz_sig = np.sqrt(rho/RHO0)
    W_S = W/S

    n_max, n_min = n_max_manobra(W), -1.0
    V_S1 = v_estol(W, S, CLmax)
    V_Sneg = v_estol(W, S, CLmin)
    V_A = V_S1*np.sqrt(n_max)
    V_C = M_C*a*raiz_sig
    M_D = M_C + margem_MD
    V_D = M_D*a*raiz_sig

    mu = razao_massa(W_S, rho, c, CLa)
    Kg = fator_alivio(mu)
    U_C = u_ref(h)
    U = {'V_B': U_C, 'V_C': U_C, 'V_D': 0.5*U_C}
    n_g_C = 1.0 + dn_rajada(V_C, U_C, W_S, CLa, Kg)
    V_B = min(V_S1*np.sqrt(n_g_C), V_C)

    V = {'V_S1': V_S1, 'V_Sneg': V_Sneg, 'V_A': V_A, 'V_B': V_B, 'V_C': V_C, 'V_D': V_D}
    dn = {k: dn_rajada(V[k], U[k], W_S, CLa, Kg) for k in U}
    return {
        'W': W, 'S': S, 'c': c, 'h': h, 'rho': rho, 'a': a,
        'CLmax': CLmax, 'CLmin': CLmin, 'CLa': CLa, 'M_C': M_C, 'M_D': M_D,
        'n_max': n_max, 'n_min': n_min, **V,
        'mu': mu, 'Kg': Kg, 'U_ref': U_C, 'U': U,
        'n_rajada_pos': {k: 1.0 + dn[k] for k in dn},
        'n_rajada_neg': {k: 1.0 - dn[k] for k in dn},
        'TAS': {k: v/raiz_sig for k, v in V.items()},
        'M': {k: v/raiz_sig/a for k, v in V.items()},
        # 25.335(a)(1): V_C >= V_B + 1,32 U_ref (nos, com U_ref em ft/s)
        'folga_VC_VB': V_C - (V_B + 1.32*(U_C/ft2m)*kt2ms),
        # n de estol em V_C e em V_D (limite do envelope quando V_A > V_D)
        'n_estol_VC': (V_C/V_S1)**2, 'n_estol_VD': (V_D/V_S1)**2,
    }


def contorno_manobra(e, n_pts=200):
    '''
    Contorno fechado do envelope de manobra: estol positivo ate
    min(V_A, V_D), n_max ate V_D, descida em V_D, -1 de V_C ate a curva de
    estol negativo e volta pela curva de estol negativo.
    '''
    V1 = np.linspace(0.0, min(e['V_A'], e['V_D']), n_pts)
    n1 = np.minimum((V1/e['V_S1'])**2, e['n_max'])
    Vs = [V1, [e['V_D'], e['V_D']]]
    ns = [n1, [n1[-1] if e['V_A'] >= e['V_D'] else e['n_max'], 0.0]]
    if e['V_A'] < e['V_D']:
        Vs.insert(1, [e['V_A'], e['V_D']])
        ns.insert(1, [e['n_max'], e['n_max']])
    Vs += [[e['V_C'], e['V_Sneg']]]
    ns += [[-1.0, -1.0]]
    V2 = np.linspace(e['V_Sneg'], 0.0, n_pts)
    Vs.append(V2)
    ns.append(-(V2/e['V_Sneg'])**2)
    return np.concatenate([np.asarray(v, float) for v in Vs]), np.concatenate([np.asarray(n, float) for n in ns])


def contorno_rajada(e):
    '''Poligono das rajadas: 1 +- dn em V_B, V_C e V_D.'''
    ks = ('V_B', 'V_C', 'V_D')
    V = [0.0] + [e[k] for k in ks] + [e[k] for k in ks[::-1]] + [0.0]
    n = ([1.0] + [e['n_rajada_pos'][k] for k in ks]
         + [e['n_rajada_neg'][k] for k in ks[::-1]] + [1.0])
    return np.array(V), np.array(n)
