import os, sys
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import otimiza_secao as ot
from analisa_otimos import RES, carrega
from poco_de_arrasto import _um_ponto

d = carrega('otim_meio')
Al, Au, a0 = ot.desmonta(d['xx_otimo'])
novos = []
for da in (-0.15, -0.05):
    for cfl in (0.10, 0.05, 0.025):
        s = _um_ponto((Al, Au, a0 + np.radians(da), cfl))
        ok = s.get('convergiu') and np.isfinite(s.get('CD', np.nan))
        print(f"d_alpha={da:+.2f} cfl={cfl}: cd={s.get('CD')} convergiu={s.get('convergiu')}", flush=True)
        if ok:
            novos.append((da, s)); break
cam = os.path.join(RES, 'poco_de_arrasto_meio.csv')
cols = ['alpha', 'CL', 'CD', 'CM', 'mach_max', 'convergiu', 'n_iter', 'residuo']
with open(cam, 'a') as fid:
    for da, s in novos:
        fid.write(f"{da:.2f}," + ','.join(
            f"{s[c]:.8g}" if isinstance(s[c], float) else str(s[c]) for c in cols) + "\n")
print(f"{len(novos)} ponto(s) acrescentado(s) a {cam}")
