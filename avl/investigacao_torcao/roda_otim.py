'''
Uma otimizacao registrada (cada avaliacao e cada iteracao do SLSQP) em
runs/<tag>.json. Uso (de qualquer pasta):
  python roda_otim.py TAG [--param full] [--margem 0.2] [--cg aft]
        [--modelo completo] [--eps 0.2] [--x0 zeros|washout|armazenado|v1,v2,...]
        [--alfa-min A]
'''

import argparse
import json
import time

import numpy as np

from comum import (Aval, otimiza, verifica, param_para_torcao, ajusta_param,
                   aeronave, LC, N_PARAM, ETAS_A, AVL_DIR)
import os

ap = argparse.ArgumentParser()
ap.add_argument('tag')
ap.add_argument('--param', default='full')
ap.add_argument('--margem', type=float, default=0.2)
ap.add_argument('--cg', default='aft')
ap.add_argument('--modelo', default='completo')
ap.add_argument('--eps', type=float, default=0.2)
ap.add_argument('--x0', default='zeros')
ap.add_argument('--alfa-min', type=float, default=None)
ap.add_argument('--sem-estol', action='store_true')
ap.add_argument('--sem-verif', action='store_true')
ap.add_argument('--maxiter', type=int, default=100)
a_ = ap.parse_args()

av = aeronave(LC)
n = N_PARAM[a_.param]
if a_.x0 == 'zeros':
    x0 = np.zeros(n)
elif a_.x0 == 'washout':
    x0 = ajusta_param(a_.param, -4.0*ETAS_A)
elif a_.x0 == 'armazenado':
    with open(os.path.join(AVL_DIR, 'resultados', 'torcao_otimizada.json')) as f:
        x0 = ajusta_param(a_.param, np.concatenate([[0.0], json.load(f)['twist']]))
else:
    x0 = np.array([float(v) for v in a_.x0.split(',')])
    if len(x0) == 11 and n != 11:
        x0 = ajusta_param(a_.param, x0)

t0 = time.time()
a = Aval(av, a_.tag, param=a_.param, margem=a_.margem, cg=a_.cg, modelo=a_.modelo,
         com_estol=not a_.sem_estol, alfa_min=a_.alfa_min)
print(f'[{a_.tag}] it fixo {a.it:.4f}; x0 {np.round(x0, 3).tolist()}', flush=True)
base = a.avalia(np.zeros(n))
r = otimiza(a, x0, eps=a_.eps, maxiter=a_.maxiter)
fim = a.avalia(r.x)
tw = param_para_torcao(a_.param, r.x)
extra = {'x0': x0.tolist(), 'eps': a_.eps, 'p_opt': r.x.tolist(), 'twist_opt': tw.tolist(),
         'CDff_opt': fim['CDff'], 'e_opt': fim['e'], 'CDff_base': base['CDff'],
         'gmin_opt': float(np.min(fim['g'])) if a.com_estol else None,
         'alfa_on_opt': fim.get('alfa_on'), 'eta_on_opt': fim.get('eta_on'),
         'margem_ext_opt': fim.get('margem_ext'),
         'sucesso': bool(r.success), 'mensagem': str(r.message), 'nit': int(r.nit),
         'nfev': len(a.evals), 'tempo_min': (time.time() - t0)/60}
a.salva(extra)
print(f'[{a_.tag}] FIM {r.message} nit {r.nit} nev {len(a.evals)} CDff '
      f'{base["CDff"]:.3f} -> {fim["CDff"]:.3f} twist {np.round(tw, 2).tolist()} '
      f'gmin {extra["gmin_opt"]}  {extra["tempo_min"]:.1f} min', flush=True)
if a_.modelo == 'completo' and not a_.sem_verif:
    extra['verif'] = verifica(tw, a_.tag)
    a.salva(extra)
    v = extra['verif']
    print(f'[{a_.tag}] VERIF CDff(aft,retrim) {v["aft"]["CDff"]:.3f} e {v["aft"]["e"]:.4f} '
          f'alfa_st {v["aft"]["alfa"]:.2f} eta {v["aft"]["eta_crit"]:.3f} CLmax {v["aft"]["CL"]:.3f} '
          f'marg_ail {v["aft"]["margem_aileron"]:.2f} folga_sombra {v["sombra"]["folga"]:.2f}',
          flush=True)
