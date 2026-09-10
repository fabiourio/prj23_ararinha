'''
O gradiente de cl_max do XFoil e utilizavel num otimizador?

Esta e a pergunta que decide se da para rodar o XFoil DENTRO do laco de
otimizacao, em vez de usar uma substituta geometrica.

O teste e o mesmo que aplicamos ao adjunto do eulerblock: calcular a derivada
por diferenca central em VARIOS passos. Se o valor for estavel entre passos,
a funcao e suave o bastante e o gradiente serve. Se variar muito, o que
estamos medindo e ruido, e o SLSQP andaria para o lado errado.

Contexto: o cl_max vem de um max sobre uma grade discreta de alpha, entao ele
e quantizado -- pula de degrau conforme o pico se desloca. A magnitude desse
degrau e o que este teste mede na pratica.

Rodar com:  python testa_gradiente_clmax.py
'''

import os
import sys
import time
from multiprocessing import Pool

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import otimiza_secao as ot
import xfoil_runner as xr
from doe_clmax_xfoil import ALPHA_SEQ, MACH_DECOLAGEM, RE_DECOLAGEM

RES = os.path.join(AQUI, 'resultados')

# passos de diferenca central, cobrindo tres ordens de grandeza
PASSOS = [3e-2, 1e-2, 3e-3, 1e-3, 3e-4]
IDX = [0, 1]              # coeficientes do extradorso a auditar


def _avalia(args):
    Al, Au, rotulo = args
    r = xr.clmax_cst(Au, Al, Re=RE_DECOLAGEM, Mach=MACH_DECOLAGEM,
                     alpha_seq=ALPHA_SEQ, timeout=600)
    return rotulo, r['clmax'], r['situacao'], r['passo_final']


def main():
    os.makedirs(RES, exist_ok=True)
    Al0, Au0 = ot.ponto_de_partida(ot.ESTACOES['meio']['tc_ref'])

    tarefas = [(Al0, Au0, ('centro', 0, 0))]
    for i in IDX:
        for h in PASSOS:
            for s in (+1, -1):
                au = np.asarray(Au0, dtype=float).copy()
                au[i] += s * h
                tarefas.append((Al0, au, (f'Au{i}', h, s)))

    print(f'{len(tarefas)} chamadas do XFoil, {len(PASSOS)} passos por '
          f'coeficiente\n', flush=True)
    t0 = time.time()
    with Pool(6) as pool:
        saidas = pool.map(_avalia, tarefas)
    print(f'{time.time()-t0:.0f} s\n')

    res = {rot: (cl, sit, passo) for rot, cl, sit, passo in saidas}
    cl0, sit0, _ = res[('centro', 0, 0)]
    print(f'ponto central: cl_max = {cl0:.4f}  ({sit0})\n')

    linhas = []
    for i in IDX:
        print(f'--- d(cl_max)/d(Au{i}) ---')
        print(f'{"passo":>9s} {"cl(+h)":>9s} {"cl(-h)":>9s} {"derivada":>11s} '
              f'  situacoes')
        derivadas = []
        for h in PASSOS:
            clp, sp, _ = res[(f'Au{i}', h, +1)]
            clm, sm, _ = res[(f'Au{i}', h, -1)]
            if not (np.isfinite(clp) and np.isfinite(clm)):
                print(f'{h:9.0e} {"--":>9s} {"--":>9s} {"nao convergiu":>11s}')
                continue
            der = (clp - clm) / (2 * h)
            derivadas.append(der)
            print(f'{h:9.0e} {clp:9.4f} {clm:9.4f} {der:11.3f}   {sp} / {sm}')
            linhas.append({'coef': f'Au{i}', 'passo': h, 'derivada': der,
                           'cl_mais': clp, 'cl_menos': clm})
        if len(derivadas) >= 2:
            d = np.array(derivadas)
            espalhamento = (d.max() - d.min()) / max(abs(d.mean()), 1e-12)
            print(f'  espalhamento entre passos: {espalhamento:.1%}', end='')
            if espalhamento < 0.10:
                print('   -> ESTAVEL, gradiente utilizavel')
            elif espalhamento < 0.40:
                print('   -> limitrofe')
            else:
                print('   -> RUIDO, gradiente inutilizavel')
        print()

    caminho = os.path.join(RES, 'gradiente_clmax.csv')
    with open(caminho, 'w') as fid:
        fid.write('coef,passo,derivada,cl_mais,cl_menos\n')
        for l in linhas:
            fid.write(f"{l['coef']},{l['passo']:.3e},{l['derivada']:.6f},"
                      f"{l['cl_mais']:.6f},{l['cl_menos']:.6f}\n")
    print(f'gravado em {caminho}')
    print('\nCriterio: para servir num otimizador de gradiente, a derivada')
    print('precisa ser estavel entre passos. Compare com o adjunto do')
    print('eulerblock, cujo pior componente variou 0,04% entre passos.')


if __name__ == '__main__':
    main()
