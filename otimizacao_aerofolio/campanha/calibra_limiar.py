'''
Calibracao final do limiar da restricao de sustentacao maxima.

Junta os dois DOEs (LHS geral + extensao para perfis grossos, que cobre a
espessura da estacao da raiz) e produz o limiar de bluntez que cada estacao
deve respeitar.

DOIS CRITERIOS, e a escolha entre eles importa:

  conservador  -- limiar com 95% de precisao: dos perfis aceitos, 95% de fato
                  atingem cl_max >= 1,80. Seguro, mas rejeita ~40% dos perfis
                  bons, e cada perfil bom rejeitado e arrasto pago a toa.
                  Chega a rejeitar o proprio NACA 1411, cujo cl_max medido e
                  1,98.

  nao-viesado  -- limiar onde o valor ESPERADO de cl_max e exatamente 1,80,
                  obtido da regressao cl_max ~ f(bluntez, t/c). Erra para os
                  dois lados em vez de sempre para o lado caro.

Adotamos o NAO-VIESADO, porque o cl_max do otimo e verificado no XFoil ao
final: se ficar abaixo de 1,80, aperta-se o limiar e re-otimiza com partida
quente. Um limiar conservador so faria sentido sem essa verificacao.

Rodar com:  python calibra_limiar.py
'''

import itertools
import os
import sys

import numpy as np
from scipy import optimize, stats

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import otimiza_secao as ot
from analise_doe_clmax import (PRECISAO_MIN, RES, le_csv, pontos_confiaveis,
                               qualidade_limiar)
from doe_clmax_xfoil import CLMAX_ALVO


def carrega_tudo():
    '''Junta os DOEs disponiveis num conjunto so.'''
    partes = []
    for nome in ('lhs', 'grossos'):
        cam = os.path.join(RES, f'doe_clmax_{nome}.csv')
        if not os.path.isfile(cam):
            print(f'  (ausente: {nome})')
            continue
        d = le_csv(cam)
        m = pontos_confiaveis(d)
        partes.append({'nome': nome,
                       'cl': d['clmax'][m],
                       't01': d['t_01'][m],
                       'tmax': d['t_max'][m]})
        print(f'  {nome}: {m.sum()} perfis, '
              f't/c de {d["t_max"][m].min():.4f} a {d["t_max"][m].max():.4f}')
    cl = np.concatenate([p['cl'] for p in partes])
    t01 = np.concatenate([p['t01'] for p in partes])
    tmax = np.concatenate([p['tmax'] for p in partes])
    return cl, t01, tmax


def superficie(bl, tmax, cl):
    '''Quadratica completa em (bluntez, t/c). Devolve preditor e R2_cv.'''
    X = np.column_stack([bl, tmax])
    mu, sd = X.mean(0), X.std(0)
    Z = (X - mu) / sd

    def base(Z):
        cols = [np.ones(len(Z)), Z[:, 0], Z[:, 1],
                Z[:, 0] ** 2, Z[:, 0] * Z[:, 1], Z[:, 1] ** 2]
        return np.column_stack(cols)

    A = base(Z)
    coef, *_ = np.linalg.lstsq(A, cl, rcond=None)
    r2 = 1 - np.sum((cl - A @ coef) ** 2) / np.sum((cl - cl.mean()) ** 2)

    rng = np.random.default_rng(23)
    ordem = rng.permutation(len(cl))
    folds = np.array_split(ordem, 5)
    e, t = [], []
    for i in range(5):
        te = folds[i]
        tr = np.concatenate([folds[j] for j in range(5) if j != i])
        c, *_ = np.linalg.lstsq(A[tr], cl[tr], rcond=None)
        e.append(np.sum((cl[te] - A[te] @ c) ** 2))
        t.append(np.sum((cl[te] - cl[tr].mean()) ** 2))
    r2cv = 1 - np.sum(e) / np.sum(t)

    def preve(bl_q, tmax_q):
        Zq = (np.array([[bl_q, tmax_q]]) - mu) / sd
        return float(base(Zq) @ coef)

    return preve, r2, r2cv


def main():
    print('=== dados ===')
    cl, t01, tmax = carrega_tudo()
    bl = t01 / np.sqrt(tmax)
    print(f'  total: {len(cl)} perfis, t/c de {tmax.min():.4f} a {tmax.max():.4f}')

    rho, _ = stats.spearmanr(bl, cl)
    print(f'  correlacao bluntez x cl_max: rho = {rho:.3f}')

    preve, r2, r2cv = superficie(bl, tmax, cl)
    print(f'  superficie cl_max ~ f(bluntez, t/c): R2 = {r2:.3f}, '
          f'R2_cv = {r2cv:.3f}')

    print(f'\n=== criterio CONSERVADOR (precisao >= {PRECISAO_MIN:.0%}) ===')
    q = qualidade_limiar(bl, cl >= CLMAX_ALVO, maior_melhor=True)
    if q:
        print(f'  bluntez >= {q[0]:.4f}   precisao {q[1]:.1%}   '
              f'cobertura {q[2]:.1%}')
        conservador = q[0]
    else:
        print('  nenhum limiar atinge a precisao exigida')
        conservador = np.nan

    print(f'\n=== criterio NAO-VIESADO (E[cl_max] = {CLMAX_ALVO}) ===')
    print('  limiar por estacao, resolvido na espessura de cada uma:')
    print(f'{"estacao":>8s} {"t/c":>8s} {"bluntez*":>10s} '
          f'{"NACA 1411 escalado":>20s} {"situacao":>12s}')

    limiares = {}
    for nome, est in ot.ESTACOES.items():
        tc = est['tc_ref']
        if not (tmax.min() <= tc <= tmax.max()):
            print(f'{nome:>8s} {tc:8.4f}   FORA da faixa amostrada '
                  f'[{tmax.min():.3f}, {tmax.max():.3f}] -- extrapolacao')
        try:
            b_star = optimize.brentq(lambda b: preve(b, tc) - CLMAX_ALVO,
                                     bl.min(), bl.max())
        except ValueError:
            b_star = np.nan
        limiares[nome] = b_star

        Al0, Au0 = ot.ponto_de_partida(tc)
        import descritores as dsc
        d0 = dsc.descritores(Au0, Al0)
        b0 = ot.bluntez(ot.t01_de(Al0, Au0), d0['t_max'])
        sit = 'atende' if b0 >= b_star else 'abaixo'
        print(f'{nome:>8s} {tc:8.4f} {b_star:10.4f} {b0:20.4f} {sit:>12s}')

    print(f'\n  (o conservador seria {conservador:.4f} para todas)')
    print('\n=== resumo para otimiza_secao.py ===')
    print('BLUNTEZ_MIN = {')
    for nome, b in limiares.items():
        print(f"    '{nome}': {b:.4f},")
    print('}')

    caminho = os.path.join(RES, 'limiares_bluntez.csv')
    with open(caminho, 'w') as fid:
        fid.write('estacao,tc_ref,bluntez_nao_viesado,bluntez_conservador\n')
        for nome, est in ot.ESTACOES.items():
            fid.write(f"{nome},{est['tc_ref']:.6f},{limiares[nome]:.6f},"
                      f'{conservador:.6f}\n')
    print(f'\ngravado em {caminho}')


if __name__ == '__main__':
    main()
