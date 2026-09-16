import os, sys, tempfile, shutil
import numpy as np
AQUI = os.path.dirname(os.path.abspath(__file__))
PACOTE = os.path.normpath(os.path.join(AQUI, '..', 'eulerblock', 'package'))
sys.path.insert(0, PACOTE); sys.path.insert(0, AQUI)
import otimiza_secao as ot
from analisa_otimos import carrega
from eulerblock import euler_mod as eb

d = carrega('otim_meio')
Ali, Aui, _ = ot.desmonta(d['xx_ini'])   # NACA 1411 reescalado p/ t/c=0.1772
alvo = 0.7319
alpha = np.radians(3.0)
tmp = tempfile.mkdtemp(prefix='trim_'); cwd = os.getcwd(); os.chdir(tmp)
try:
    for it in range(6):
        r = eb.run_cst(Ali, Aui, ot.NCHORD, alpha, ot.MACH_N, gamma=1.4,
                       order=2, iter=ot.ITER, dt=ot.DT, CFL=0.2,
                       use_local_dt=1, res_NK=ot.RES_NK, res_tol=ot.RES_TOL,
                       reinitialize=0, adj_funcs=['cl_jlow'], plot=False,
                       NJ=ot.NJ, s0=ot.S0)
        print(f'iter {it}: alpha = {np.degrees(alpha):.4f} deg  '
              f'cl = {r["CL"]:.5f}  cd = {r["CD"]:.6f}  cm = {r["CM"]:.5f}',
              flush=True)
        err = r['CL'] - alvo
        if abs(err) < 1e-3:
            print(f'\nTRIMADO: alpha = {np.degrees(alpha):.4f} deg, '
                  f'cl = {r["CL"]:.5f}, cd = {r["CD"]:.6f}, cm = {r["CM"]:.5f}')
            break
        alpha -= err / r['grads']['cl_jlow']['alpha']
finally:
    os.chdir(cwd); shutil.rmtree(tmp, ignore_errors=True)
