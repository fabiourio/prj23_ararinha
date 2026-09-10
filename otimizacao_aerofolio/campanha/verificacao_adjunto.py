'''
Verificacao dos gradientes adjuntos do eulerblock -- Lab 03, PRJ-23.

O SLSQP vai confiar inteiramente nos gradientes que o metodo adjunto entrega.
Se estiverem errados, a otimizacao produz lixo depois de horas de maquina.
Conferir custa poucos minutos, entao nao ha desculpa.

Metodo: diferenca central em torno do ponto de partida, com res_tol apertado
(1e-8) para que o ruido do solver nao domine a subtracao. Comparamos
dc_d/dalpha, dc_l/dalpha e dc_d/dA_u[i] para alguns i.

E a aplicacao direta da Aula 03, que apresenta diferencas finitas, passo
complexo, diferenciacao automatica e metodo adjunto como formas de obter
gradientes: aqui usamos a primeira para auditar a ultima.

ATENCAO a duas armadilhas ja documentadas:
  - run_cst(Al, Au, ...) recebe o INTRADORSO primeiro, ao contrario de
    cstfoil(Au, Al, ...);
  - cada avaliacao precisa de diretorio proprio, porque o eulerblock.exe
    escreve com nomes fixos no diretorio corrente.

Rodar com:  python verificacao_adjunto.py
Saida:      resultados/verificacao_adjunto.csv
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

MACH_N = 0.7196
ALPHA0 = 2.0 * np.pi / 180.0

AU = np.array([0.16146332, 0.18349204, 0.14126241, 0.18194397])
AL = np.array([-0.1489439, -0.10330027, -0.10305128, -0.10514982])

NCHORD, NJ, S0 = 31, 49, 0.5e-2      # nivel 1,0 -- o mesmo da otimizacao

# Passos da diferenca central. Grandes demais captam a nao-linearidade;
# pequenos demais deixam o ruido do solver dominar. Testamos dois de cada
# para poder distinguir os dois regimes.
PASSOS_ALPHA = [1e-3, 1e-4]          # rad
PASSOS_CST = [1e-3, 1e-4]

# quais coeficientes do extradorso auditar (0 = o que controla o raio do nariz)
IDX_AU = [0, 2]


def _avalia(args):
    '''Roda um ponto e devolve (CL, CD) -- e, se pedido, os gradientes.'''
    al, au, alpha, com_adjunto = args
    from eulerblock import euler_mod as eb

    tmp = tempfile.mkdtemp(prefix='adj_')
    cwd = os.getcwd()
    try:
        os.chdir(tmp)
        adj = ['cl_jlow', 'cd_jlow'] if com_adjunto else []
        r = eb.run_cst(al, au, NCHORD, alpha, MACH_N,
                       gamma=1.4, order=2,
                       iter=20000, dt=0.001, CFL=0.2, use_local_dt=1,
                       res_NK=1e-5, res_tol=1e-8,
                       reinitialize=0, adj_funcs=adj, plot=False,
                       NJ=NJ, s0=S0)
        maxt = float(np.real(r['maxt']))
        if not 0.05 < maxt < 0.25:
            raise ValueError(f'espessura {maxt:.4f} suspeita -- ordem de Al/Au?')
        saida = {'CL': float(r['CL']), 'CD': float(r['CD'])}
        if com_adjunto:
            g = r['grads']
            saida['dCL_dalpha'] = float(g['cl_jlow']['alpha'])
            saida['dCD_dalpha'] = float(g['cd_jlow']['alpha'])
            saida['dCL_dAu'] = np.asarray(g['cl_jlow']['dAu'], dtype=float)
            saida['dCD_dAu'] = np.asarray(g['cd_jlow']['dAu'], dtype=float)
        return saida
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    os.makedirs(RES, exist_ok=True)
    t0 = time.time()

    # monta a lista de avaliacoes: o ponto central (com adjunto) e os pares
    # perturbados de cada diferenca central
    tarefas = [(AL, AU, ALPHA0, True)]
    rotulos = [('centro', None, None)]

    for h in PASSOS_ALPHA:
        tarefas.append((AL, AU, ALPHA0 + h, False))
        rotulos.append(('alpha', h, +1))
        tarefas.append((AL, AU, ALPHA0 - h, False))
        rotulos.append(('alpha', h, -1))

    for i in IDX_AU:
        for h in PASSOS_CST:
            for s in (+1, -1):
                au = AU.copy()
                au[i] += s * h
                tarefas.append((AL, au, ALPHA0, False))
                rotulos.append((f'Au{i}', h, s))

    print(f'{len(tarefas)} avaliacoes do Euler '
          f'(1 com adjunto + {len(tarefas)-1} para diferencas centrais)')
    print('malha nivel 1,0; res_tol = 1e-8 para o ruido nao dominar\n', flush=True)

    with Pool(4) as pool:
        saidas = pool.map(_avalia, tarefas)

    centro = saidas[0]
    print(f"ponto central: CL = {centro['CL']:.6f}  CD = {centro['CD']:.6f}")
    print(f"  adjunto: dCL/dalpha = {centro['dCL_dalpha']:+.6e}")
    print(f"           dCD/dalpha = {centro['dCD_dalpha']:+.6e}")
    print(f"           dCD/dAu    = {np.array2string(centro['dCD_dAu'], precision=4)}\n")

    # agrupa os pares (+h, -h)
    pares = {}
    for rot, s in zip(rotulos[1:], saidas[1:]):
        nome, h, sinal = rot
        pares.setdefault((nome, h), {})[sinal] = s

    linhas = []
    print(f'{"quantidade":>12s} {"passo":>8s} {"adjunto":>14s} '
          f'{"dif. central":>14s} {"erro rel.":>11s}')
    for (nome, h), par in pares.items():
        if +1 not in par or -1 not in par:
            continue
        for fun in ('CL', 'CD'):
            df = (par[+1][fun] - par[-1][fun]) / (2 * h)
            if nome == 'alpha':
                adj = centro[f'd{fun}_dalpha']
            else:
                i = int(nome[2:])
                adj = centro[f'd{fun}_dAu'][i]
            erro = abs(df - adj) / max(abs(adj), 1e-14)
            marca = '  OK' if erro < 0.05 else ('  ATENCAO' if erro < 0.25
                                                else '  DIVERGE')
            print(f'd{fun}/d{nome:<7s} {h:8.0e} {adj:+14.6e} {df:+14.6e} '
                  f'{erro:10.2%}{marca}')
            linhas.append({'funcao': fun, 'variavel': nome, 'passo': h,
                           'adjunto': adj, 'dif_central': df, 'erro_rel': erro})

    caminho = os.path.join(RES, 'verificacao_adjunto.csv')
    cols = ['funcao', 'variavel', 'passo', 'adjunto', 'dif_central', 'erro_rel']
    with open(caminho, 'w') as fid:
        fid.write(','.join(cols) + '\n')
        for l in linhas:
            fid.write(','.join(
                f'{l[c]:.10g}' if isinstance(l[c], float) else str(l[c])
                for c in cols) + '\n')
    print(f'\ngravado em {caminho}')
    print(f'tempo total: {time.time() - t0:.0f} s')

    piores = [l for l in linhas if l['erro_rel'] >= 0.25]
    if piores:
        print('\nATENCAO: ha gradientes divergindo. NAO rode a otimizacao')
        print('antes de entender o motivo -- o SLSQP confia neles cegamente.')
    else:
        print('\nGradientes conferem. A otimizacao pode usar o adjunto.')


if __name__ == '__main__':
    main()
