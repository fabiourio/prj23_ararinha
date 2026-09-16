import json
import os
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
PACOTE = os.path.normpath(os.path.join(AQUI, '..', 'eulerblock', 'package'))
sys.path.insert(0, PACOTE)
sys.path.insert(0, AQUI)

from eulerblock import euler_mod as eb                          # noqa: E402
import otimiza_secao as ot                                      # noqa: E402

TOL_CL = 1e-4          # casamento do c_l; o ruido do solver e ~1e-5 no c_d
MAX_AVAL = 8


def acha_alpha(estacao='meio'):
    est = ot.ESTACOES[estacao]
    cl_ref, tc_ref = est['cl_ref'], est['tc_ref']
    cfl = est.get('cfl', 0.20)
    Al, Au = ot.ponto_de_partida(tc_ref)

    pasta = os.path.join(AQUI, 'resultados', f'ativ2_{estacao}')
    os.makedirs(pasta, exist_ok=True)
    cwd = os.getcwd()

    def roda(alpha_deg):
        try:
            os.chdir(pasta)
            r = eb.run_cst(Al, Au, ot.NCHORD, alpha_deg * np.pi / 180,
                           ot.MACH_N, gamma=1.4, order=2, iter=ot.ITER,
                           dt=ot.DT, CFL=cfl, use_local_dt=1,
                           res_NK=ot.RES_NK, res_tol=ot.RES_TOL,
                           reinitialize=0, plot=False,
                           adj_funcs=['cl_jlow', 'cd_jlow'],
                           NJ=ot.NJ, s0=ot.S0)
        finally:
            os.chdir(cwd)
        return (float(r['CL']), float(r['CD']),
                float(r['grads']['cl_jlow']['alpha']),
                float(np.real(r['maxt'])))

    print(f'=== item 2 -- estacao {estacao} ===')
    print(f'  M_n = {ot.MACH_N}   alvo c_l = {cl_ref}   '
          f't/c da partida = {tc_ref}')

    hist = []
    alpha, t0 = 3.0, time.time()
    for i in range(MAX_AVAL):
        cl, cd, dcl_da, maxt = roda(alpha)
        hist.append(dict(alpha_deg=alpha, cl=cl, cd=cd, maxt=maxt))
        print(f'  [{i + 1}] alpha = {alpha:8.4f} deg   c_l = {cl:.6f}   '
              f'c_d = {cd:.6f}   dc_l/dalpha = {dcl_da:.4f}/rad', flush=True)
        if abs(cl - cl_ref) < TOL_CL:
            break
        # passo de Newton; o gradiente vem por radiano
        alpha += (cl_ref - cl) / (dcl_da * np.pi / 180)

    print(f'\n  RESULTADO: alpha = {hist[-1]["alpha_deg"]:.6f} graus, '
          f'c_l = {hist[-1]["cl"]:.6f}, c_d = {hist[-1]["cd"]:.6f}')
    print(f'  ({len(hist)} avaliacoes, {time.time() - t0:.0f} s)')

    saida = dict(estacao=estacao, mach_n=ot.MACH_N, cl_ref=cl_ref,
                 tc_ref=tc_ref, hist=hist,
                 alpha_final_deg=hist[-1]['alpha_deg'],
                 cl_final=hist[-1]['cl'], cd_final=hist[-1]['cd'],
                 Al=list(Al), Au=list(Au))
    nome = ('ativ2_alpha_naca.json' if estacao == 'meio'
            else f'ativ2_alpha_naca_{estacao}.json')
    cam = os.path.join(AQUI, 'resultados', nome)
    with open(cam, 'w') as fid:
        json.dump(saida, fid, indent=2)
    print(f'  gravado em {os.path.relpath(cam, AQUI)}')
    return saida


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    est = args[0] if args else 'meio'
    if est not in ot.ESTACOES:
        print(f'estacao invalida: {est}. Use: {list(ot.ESTACOES)}')
        sys.exit(1)
    acha_alpha(est)
