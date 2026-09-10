'''
Verificacao de malha do codigo de Euler -- Lab 03, PRJ-23.

Ninguem deve gastar horas otimizando com um solver nao verificado. Este
script mede como c_d e c_l variam com o refinamento, na CONDICAO DE PROJETO
da nossa asa (M_n = 0,7196, normal ao enflechamento de 50% da corda), e nao
na condicao de exemplo do professor.

Ja e sabido que a malha pela metade e inaceitavel: no exemplo do professor
(M = 0,75) o c_d vai de 0,01577 para 0,03557, +126%, embora o caso rode em
6 s em vez de 68 s. O objetivo aqui e achar onde o c_d estabiliza e fixar o
nivel com justificativa numerica, nao por gosto.

Cada nivel roda em diretorio proprio: o eulerblock.exe escreve grid.xyz,
settings.txt, wall.dat, solution.vtk e derivatives.dat com nomes fixos no
diretorio corrente.

Rodar com:  python verificacao_malha.py
Saida:      resultados/verificacao_malha.csv
'''

import os
import shutil
import sys
import tempfile
import time
from multiprocessing import Pool

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
PACOTE = os.path.normpath(os.path.join(AQUI, '..', 'eulerblock', 'package'))
sys.path.insert(0, PACOTE)
sys.path.insert(0, AQUI)

RES = os.path.join(AQUI, 'resultados')

# --- condicao de projeto (ver documento de projeto, secoes 3.2 e 4) ---------
MACH_N = 0.7196          # M * cos(Lambda_c/2), Lambda_c/2 = 32,158 graus
ALPHA_DEG = 2.0          # angulo fixo para a comparacao entre malhas

# NACA 1411 -- ponto de partida do roteiro
AU = [0.16146332, 0.18349204, 0.14126241, 0.18194397]
AL = [-0.1489439, -0.10330027, -0.10305128, -0.10514982]

# malha de referencia do professor (single_run.py)
NCHORD0, NJ0, S00 = 30, 48, 0.5e-2

NIVEIS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]


def roda_nivel(nivel):
    '''Roda um nivel de malha em diretorio isolado.'''
    from eulerblock import euler_mod as eb

    nchord = int(NCHORD0 * nivel) + 1
    nj = int(NJ0 * nivel) + 1
    s0 = NJ0 / (nj - 1) * S00

    tmp = tempfile.mkdtemp(prefix=f'malha_{nivel:.2f}_')
    cwd = os.getcwd()
    try:
        os.chdir(tmp)
        t0 = time.time()
        # ATENCAO a ordem dos argumentos: run_cst(Al, Au, ...) recebe o
        # INTRADORSO primeiro, ao contrario de cstfoil(Au, Al, ...), que
        # recebe o extradorso primeiro. Inverter gera um perfil de dentro
        # para fora, o gerador de malha produz celulas de area negativa e o
        # solver morre sem escrever wall.dat.
        r = eb.run_cst(AL, AU, nchord, ALPHA_DEG * np.pi / 180.0, MACH_N,
                       gamma=1.4, order=2,
                       iter=20000, dt=0.001, CFL=0.2, use_local_dt=1,
                       res_NK=1e-5, res_tol=1e-8,
                       reinitialize=0, adj_funcs=[], plot=False,
                       NJ=nj, s0=s0)
        dt = time.time() - t0
        maxt = float(np.real(r['maxt']))
        if not 0.05 < maxt < 0.25:
            raise ValueError(
                f'espessura maxima {maxt:.5f} fora do esperado para o '
                'NACA 1411 -- verifique a ordem de Al/Au')
        return {'nivel': nivel, 'nchord': nchord, 'nj': nj,
                'celulas': (nchord - 1) * (nj - 1),
                'CL': r['CL'], 'CD': r['CD'], 'CM': r['CM'],
                'tempo_s': dt}
    except Exception as exc:                                   # noqa: BLE001
        return {'nivel': nivel, 'nchord': nchord, 'nj': nj,
                'celulas': (nchord - 1) * (nj - 1),
                'CL': np.nan, 'CD': np.nan, 'CM': np.nan,
                'tempo_s': np.nan, 'erro': f'{type(exc).__name__}: {exc}'}
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    os.makedirs(RES, exist_ok=True)
    print(f'Verificacao de malha: NACA 1411, M_n = {MACH_N}, '
          f'alpha = {ALPHA_DEG} graus')
    print(f'{len(NIVEIS)} niveis, 2 processos '
          f'(o DOE do XFoil pode estar usando o resto da maquina)\n')

    t0 = time.time()
    with Pool(2) as pool:
        linhas = list(pool.imap(roda_nivel, NIVEIS))
    linhas.sort(key=lambda l: l['nivel'])

    print(f'{"nivel":>6s} {"nchord":>7s} {"NJ":>5s} {"celulas":>9s} '
          f'{"CL":>9s} {"CD":>10s} {"CM":>9s} {"tempo":>8s} {"dCD vs ref":>11s}')
    ref = next((l['CD'] for l in linhas if l['nivel'] == max(NIVEIS)), np.nan)
    for l in linhas:
        d = (l['CD'] - ref) / ref * 100 if np.isfinite(l['CD']) and ref else np.nan
        print(f"{l['nivel']:6.2f} {l['nchord']:7d} {l['nj']:5d} "
              f"{l['celulas']:9d} {l['CL']:9.5f} {l['CD']:10.6f} "
              f"{l['CM']:9.5f} {l['tempo_s']:7.1f}s {d:10.2f}%"
              + (f"   {l['erro']}" if 'erro' in l else ''))

    caminho = os.path.join(RES, 'verificacao_malha.csv')
    cols = ['nivel', 'nchord', 'nj', 'celulas', 'CL', 'CD', 'CM', 'tempo_s']
    with open(caminho, 'w') as fid:
        fid.write(','.join(cols) + '\n')
        for l in linhas:
            fid.write(','.join(f'{l[c]:.8g}' for c in cols) + '\n')
    print(f'\ngravado em {caminho}')
    print(f'tempo total: {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
