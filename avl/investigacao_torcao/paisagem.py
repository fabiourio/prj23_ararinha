'''
B1/B2: paisagem do objetivo em torno do otimo armazenado e variantes
suavizadas, avaliadas direto (CDff no it fixo de producao e restricao de
estol original). Saida: runs/paisagem.json
'''

import json
import os

import numpy as np

from comum import Aval, aeronave, LC, ETAS_A, RUNS, AVL_DIR, param_para_torcao, ajusta_param

av = aeronave(LC)
a = Aval(av, 'paisagem')
with open(os.path.join(AVL_DIR, 'resultados', 'torcao_otimizada.json')) as f:
    xs = np.array(json.load(f)['twist'])
tw_s = np.concatenate([[0.0], xs])
r0 = a.avalia(xs)
f0, g0 = r0['CDff'], r0['g'].min()
print(f'otimo armazenado: CDff {f0:.3f}  gmin {g0:.4f}', flush=True)

# ---------------------------------------------------------------- B1: gradiente e curvatura por variavel
H = 0.5
grad = []
for i in range(10):
    lo_ok = xs[i] - H >= -8.0
    hi_ok = xs[i] + H <= 3.0 + 1e-9
    if lo_ok and hi_ok:
        pts = [-H, +H]
    else:
        pts = [-H, -2*H]          # no limite superior: unilateral para dentro
    vals = []
    for d in pts:
        x = xs.copy(); x[i] += d
        r = a.avalia(x)
        vals.append((d, r['CDff'], float(r['g'].min()), r['alfa_on'], r['eta_on']))
    (d1, f1, g1, *_), (d2, f2, g2, *_) = vals
    # ajuste quadratico por (0,f0), (d1,f1), (d2,f2)
    A = np.array([[d1, d1**2], [d2, d2**2]])
    b_, c_ = np.linalg.solve(A, [f1 - f0, f2 - f0])
    gb, gc = np.linalg.solve(A, [g1 - g0, g2 - g0])
    grad.append({'i': i, 'eta': float(ETAS_A[i + 1]), 'x': float(xs[i]),
                 'dCDff_dx': float(b_), 'd2CDff_dx2': float(2*c_),
                 'dg_dx': float(gb), 'pontos': vals})
    print(f'  eta {ETAS_A[i+1]:.3f}: dCD/dx {b_:+.3f} count/°  d2CD/dx2 {2*c_:+.3f} count/°²  '
          f'dgmin/dx {gb:+.3f}  pts {[(round(v[0],2), round(v[1],3), round(v[2],3)) for v in vals]}',
          flush=True)

# ---------------------------------------------------------------- B2: variantes suavizadas
def interp_sobre(tw, i_de, i_ate):
    t = tw.copy()
    e = ETAS_A
    t[i_de + 1:i_ate] = np.interp(e[i_de + 1:i_ate], [e[i_de], e[i_ate]], [t[i_de], t[i_ate]])
    return t

var = {
    'armazenado': tw_s,
    'interior_linear(0.10-0.56)': interp_sobre(tw_s, 1, 6),
    'sem_pico_0.48': interp_sobre(tw_s, 4, 6),
    'vale_preenchido(0.10-0.48)': interp_sobre(tw_s, 1, 5),
    'sem_ondas_0.22-0.40(0.10->0.48 lin, mantem pico)': interp_sobre(tw_s, 1, 5),
    'ajuste_pchip4': param_para_torcao('pchip4', ajusta_param('pchip4', tw_s)),
    'ajuste_pl3': param_para_torcao('pl3', ajusta_param('pl3', tw_s)),
}
# 'sem ondas' e 'vale preenchido' coincidem; mantem so um
var.pop('sem_ondas_0.22-0.40(0.10->0.48 lin, mantem pico)')

res_var = {}
for nome, tw in var.items():
    r = a.avalia(tw[1:])
    ent = {'twist': tw.tolist(), 'CDff': r['CDff'], 'gmin': float(r['g'].min()),
           'alfa_on': r['alfa_on'], 'eta_on': r['eta_on'],
           'alfa_min_all': r['alfa_min_all'], 'eta_min_all': r['eta_min_all']}
    # deslocamento uniforme da parte externa (eta >= 0.56) que recupera g >= 0
    if ent['gmin'] < 0:
        ext = ETAS_A[1:] >= 0.56 - 1e-9
        d_a, g_a = 0.0, ent['gmin']
        d_b = -max(0.5, -g_a)
        for _ in range(6):
            x = tw[1:].copy(); x[ext] += d_b
            rb = a.avalia(x); g_b = float(rb['g'].min())
            if abs(g_b) < 0.02:
                break
            d_c = d_b - g_b*(d_b - d_a)/(g_b - g_a)
            d_a, g_a, d_b = d_b, g_b, d_c
        ent['reparo'] = {'delta_ext': float(d_b), 'CDff': rb['CDff'], 'gmin': g_b,
                         'alfa_on': rb['alfa_on'], 'eta_on': rb['eta_on'],
                         'twist': np.concatenate([[0.0], x]).tolist()}
    res_var[nome] = ent
    rep = ent.get('reparo')
    print(f'  {nome:32s} CDff {ent["CDff"]:.3f} (Δ {ent["CDff"]-f0:+.3f})  gmin {ent["gmin"]:+.3f}  '
          f'onset α {ent["alfa_on"]:.2f} η {ent["eta_on"]:.3f}  first-anywhere η {ent["eta_min_all"]:.3f}'
          + (f'  | reparo Δext {rep["delta_ext"]:+.2f}° -> CDff {rep["CDff"]:.3f} '
             f'(Δ {rep["CDff"]-f0:+.3f}) gmin {rep["gmin"]:+.3f}' if rep else ''), flush=True)

with open(os.path.join(RUNS, 'paisagem.json'), 'w') as f:
    json.dump({'f0': f0, 'g0': float(g0), 'H': H, 'grad': grad, 'variantes': res_var},
              f, indent=1, default=float)
print('ok')
