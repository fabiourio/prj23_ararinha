'''
Fecha a investigacao do poco de arrasto: duas perguntas que sobraram.

A varredura fina (poco_de_arrasto.py) mostrou que o poco EXISTE e e estreito
-- os dois pontos que o cercam convergiram direito --, mas deixou duas
pontas soltas:

  1. SEIS de 21 pontos estouraram o limite de 20000 iteracoes sem atingir a
     tolerancia. O roda_euler() do pos-processamento aceita qualquer valor
     FINITO, sem olhar o residuo, entao esses pontos entraram na polar da
     entrega com qualidade desconhecida. Aqui repetimos so eles, recuando o
     CFL como o otimizador faz, e vemos se o valor muda.

  2. O poco sobrevive ao refinamento de malha? Esta e a pergunta que decide
     se ele e fisica ou artefato do nivel 1,0. Se o poco sumir na malha fina,
     o otimizador explorou uma peculiaridade do discreto e o perfil nao tem o
     desempenho que a polar promete. Reavaliamos cinco pontos em torno do
     otimo no nivel 1,5.

Rodar com:  python poco_confirmacao.py [estacao]
'''

import os
import sys
from multiprocessing import Pool

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
PACOTE = os.path.normpath(os.path.join(AQUI, '..', 'eulerblock', 'package'))
sys.path.insert(0, PACOTE)
sys.path.insert(0, AQUI)

import otimiza_secao as ot
from analisa_otimos import RES, carrega
from poco_de_arrasto import _um_ponto

RUINS = [-0.7, -0.1, 0.6, 0.7, 0.8, 0.9]      # os que nao convergiram
MALHA = [-0.4, -0.2, 0.0, 0.2, 0.4]           # vizinhanca para o teste de malha
CFL_RECUO = [0.10, 0.05, 0.025]
NIVEL_FINO = 1.5


def _ponto_malha(args):
    '''Um ponto no nivel de malha pedido, com recuo de CFL.'''
    import shutil
    import tempfile

    from eulerblock import euler_mod as eb

    Al, Au, alpha, nivel = args
    nchord, nj = int(30 * nivel) + 1, int(48 * nivel) + 1
    s0 = 48 / (nj - 1) * 0.5e-2
    tmp = tempfile.mkdtemp(prefix='pocom_')
    cwd = os.getcwd()
    try:
        os.chdir(tmp)
        for cfl in [0.20] + CFL_RECUO:
            r = eb.run_cst(Al, Au, nchord, alpha, ot.MACH_N, gamma=1.4,
                           order=2, iter=ot.ITER, dt=ot.DT, CFL=cfl,
                           use_local_dt=1, res_NK=ot.RES_NK,
                           res_tol=ot.RES_TOL, reinitialize=0, adj_funcs=[],
                           plot=False, NJ=nj, s0=s0)
            if np.isfinite(r['CL']) and np.isfinite(r['CD']):
                return dict(alpha=alpha, nivel=nivel, cfl=cfl,
                            CL=float(r['CL']), CD=float(r['CD']))
        return dict(alpha=alpha, nivel=nivel, cfl=np.nan, CL=np.nan,
                    CD=np.nan)
    except Exception:                                          # noqa: BLE001
        return dict(alpha=alpha, nivel=nivel, cfl=np.nan, CL=np.nan,
                    CD=np.nan)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    estacao = sys.argv[1] if len(sys.argv) > 1 else 'meio'
    d = carrega(f'otim_{estacao}')
    Al, Au, a0 = ot.desmonta(d['xx_otimo'])

    # ---- pergunta 1: os pontos nao-convergidos mudam com CFL menor? ----
    print('=' * 74)
    print('1. Os pontos que nao convergiram mudam ao recuar o CFL?')
    print('=' * 74, flush=True)
    antes = {}
    cam_fino = os.path.join(RES, f'poco_de_arrasto_{estacao}.csv')
    with open(cam_fino) as fid:
        cab = fid.readline().strip().split(',')
        for linha in fid:
            v = linha.strip().split(',')
            reg = dict(zip(cab, v))
            antes[round(float(reg['d_alpha_deg']), 1)] = float(reg['CD'])

    tarefas = [(Al, Au, a0 + np.radians(da), c)
               for da in RUINS for c in CFL_RECUO]
    with Pool(4) as pool:
        saidas = pool.map(_um_ponto, tarefas)

    melhor = {}
    for s, (_, _, _, cfl) in zip(saidas, tarefas):
        da = round(np.degrees(s['alpha'] - a0), 1)
        if s['convergiu'] and np.isfinite(s['CD']):
            melhor.setdefault(da, (cfl, s['CD'], s['CL']))

    print(f'{"d_alpha":>8s} {"c_d antes":>10s} {"c_d convergido":>15s} '
          f'{"variacao":>10s}  CFL')
    for da in RUINS:
        if da in melhor:
            cfl, cd, cl = melhor[da]
            print(f'{da:+8.1f} {antes[da]:10.6f} {cd:15.6f} '
                  f'{(cd/antes[da]-1)*100:+9.2f}%  {cfl}')
        else:
            print(f'{da:+8.1f} {antes[da]:10.6f} {"nao convergiu":>15s}')

    # ---- pergunta 2: o poco sobrevive ao refinamento? ----
    print(f'\n{"="*74}')
    print(f'2. O poco sobrevive na malha nivel {NIVEL_FINO}?')
    print('=' * 74, flush=True)
    tarefas = [(Al, Au, a0 + np.radians(da), NIVEL_FINO) for da in MALHA]
    with Pool(3) as pool:
        finos = pool.map(_ponto_malha, tarefas)

    print(f'{"d_alpha":>8s} {"c_d nivel 1,0":>14s} {"c_d nivel 1,5":>14s} '
          f'{"c_l 1,5":>9s}')
    linhas = []
    for da, s in zip(MALHA, finos):
        print(f'{da:+8.1f} {antes[da]:14.6f} {s["CD"]:14.6f} '
              f'{s["CL"]:9.5f}')
        linhas.append({'d_alpha_deg': da, 'CD_nivel_1.0': antes[da],
                       'CD_nivel_1.5': s['CD'], 'CL_nivel_1.5': s['CL']})

    cds = np.array([s['CD'] for s in finos])
    if np.all(np.isfinite(cds)):
        i = int(np.argmin(cds))
        prof_10 = max(antes[d] for d in MALHA) / min(antes[d] for d in MALHA)
        prof_15 = np.nanmax(cds) / np.nanmin(cds)
        print(f'\n  profundidade do poco (c_d max / c_d min na faixa):')
        print(f'    nivel 1,0: {prof_10:.2f}x')
        print(f'    nivel 1,5: {prof_15:.2f}x')
        print(f'  minimo no nivel 1,5 em d_alpha = {MALHA[i]:+.1f} deg')
        if abs(MALHA[i]) > 1e-9:
            print('  ATENCAO: o minimo na malha fina NAO esta no ponto de')
            print('  projeto. O otimizador mirou onde a malha grossa marcava.')

    cam = os.path.join(RES, f'poco_malha_{estacao}.csv')
    with open(cam, 'w') as fid:
        cols = list(linhas[0])
        fid.write(','.join(cols) + '\n')
        for l in linhas:
            fid.write(','.join(f'{l[c]:.8g}' for c in cols) + '\n')
    print(f'\n  gravado em {cam}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
