'''
Analise do DOE de cl_max -- escolhe a restricao substituta.

Responde, com os dados de doe_clmax_xfoil.py, tres perguntas em ordem:

  1. Qual descritor geometrico melhor prediz cl_max? (correlacao de posto,
     que capta relacao monotona mesmo quando nao e linear)

  2. Existe um LIMIAR de um descritor so que garanta cl_max >= 1.8 com
     confianca? Isso e o ideal: vira um simples batente no otimizador, de
     custo zero. Medimos precisao (dos perfis aceitos, quantos de fato tem
     cl_max >= alvo) e cobertura (dos perfis bons, quantos sobrevivem ao
     limiar). Um limiar so serve se for preciso E nao jogar fora quase tudo.

  3. Se nenhum descritor isolado servir, uma superficie de resposta quadratica
     sobre poucos descritores da conta? Ela continua analitica e suave, logo
     continua utilizavel como g(x) num otimizador de gradiente, e continua
     custando microssegundos -- nenhuma chamada de XFoil dentro do laco.

Rodar com:  python analise_doe_clmax.py
'''

import itertools
import os
import sys

import numpy as np
from scipy import stats

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import descritores as dsc
from doe_clmax_xfoil import CLMAX_ALVO, CLMAX_VIABILIDADE

RES = os.path.join(AQUI, 'resultados')

# precisao minima exigida de um limiar para considera-lo utilizavel
PRECISAO_MIN = 0.95
# cobertura minima: abaixo disso o limiar e preciso mas inutil (corta demais)
COBERTURA_MIN = 0.30


def le_csv(caminho):
    '''Leitor simples: devolve dict de colunas (float quando possivel).'''
    with open(caminho) as fid:
        linhas = [l.rstrip('\n') for l in fid if l.strip()]
    cab = linhas[0].split(',')
    dados = {c: [] for c in cab}
    for l in linhas[1:]:
        partes = l.split(',')
        for c, v in zip(cab, partes):
            dados[c].append(v)
    saida = {}
    for c, vals in dados.items():
        try:
            saida[c] = np.array([float(v) if v != '' else np.nan for v in vals])
        except ValueError:
            saida[c] = np.array(vals, dtype=object)
    return saida


def pontos_confiaveis(d):
    '''So perfis avaliados com estol capturado de forma confiavel.'''
    ok = np.array([s == 'ok' for s in d['status']])
    conf = d['confiavel'] == 1
    fin = np.isfinite(d['clmax'])
    return ok & conf & fin


def qualidade_limiar(valores, bom, maior_melhor=True):
    '''
    Varre limiares e devolve o melhor com precisao >= PRECISAO_MIN.
    Devolve (limiar, precisao, cobertura) ou None se nenhum servir.
    '''
    candidatos = np.unique(valores)
    melhor = None
    for t in candidatos:
        aceito = valores >= t if maior_melhor else valores <= t
        if aceito.sum() < 10:
            continue
        precisao = bom[aceito].mean()
        cobertura = aceito[bom].mean() if bom.sum() else 0.0
        if precisao >= PRECISAO_MIN and (melhor is None or cobertura > melhor[2]):
            melhor = (float(t), float(precisao), float(cobertura))
    return melhor


def superficie_quadratica(X, y, n_folds=5, semente=23):
    '''
    Ajusta y ~ quadratica completa em X (padronizado) e devolve R2 e R2 de
    validacao cruzada. O R2 de treino sempre sobe com mais termos; o de
    validacao e o que diz se o modelo generaliza.
    '''
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd

    def base(Z):
        cols = [np.ones(len(Z))]
        cols += [Z[:, i] for i in range(Z.shape[1])]
        for i, j in itertools.combinations_with_replacement(range(Z.shape[1]), 2):
            cols.append(Z[:, i] * Z[:, j])
        return np.column_stack(cols)

    A = base(Z)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    r2 = 1 - np.sum((y - A @ coef) ** 2) / np.sum((y - y.mean()) ** 2)

    rng = np.random.default_rng(semente)
    ordem = rng.permutation(len(y))
    folds = np.array_split(ordem, n_folds)
    erros, totais = [], []
    for ii in range(n_folds):
        teste = folds[ii]
        treino = np.concatenate([folds[j] for j in range(n_folds) if j != ii])
        c, *_ = np.linalg.lstsq(A[treino], y[treino], rcond=None)
        erros.append(np.sum((y[teste] - A[teste] @ c) ** 2))
        totais.append(np.sum((y[teste] - y[treino].mean()) ** 2))
    r2_cv = 1 - np.sum(erros) / np.sum(totais)
    return r2, r2_cv, coef, mu, sd


def main():
    cam_lhs = os.path.join(RES, 'doe_clmax_lhs.csv')
    cam_corte = os.path.join(RES, 'doe_clmax_corte.csv')
    if not os.path.isfile(cam_lhs):
        print(f'ERRO: {cam_lhs} nao existe. Rode doe_clmax_xfoil.py antes.')
        return 1

    # ---------- estudo A: corte controlado ----------
    if os.path.isfile(cam_corte):
        c = le_csv(cam_corte)
        m = pontos_confiaveis(c)
        print('=== ESTUDO A: corte controlado (so o nariz varia) ===')
        print(f'{"r_LE_sup":>10s} {"delta_y":>9s} {"t_max":>8s} {"clmax":>8s} '
              f'{"a_stall":>8s}')
        for ii in np.where(m)[0]:
            print(f"{c['r_LE_sup'][ii]:10.5f} {c['delta_y'][ii]:9.3f} "
                  f"{c['t_max'][ii]:8.4f} {c['clmax'][ii]:8.4f} "
                  f"{c['alpha_clmax'][ii]:8.1f}")
        if m.sum() > 2:
            rho, p = stats.spearmanr(c['r_LE_sup'][m], c['clmax'][m])
            print(f'\n  correlacao de posto r_LE x cl_max no corte: '
                  f'rho = {rho:.3f} (p = {p:.2g})')

    # ---------- estudo B: LHS ----------
    d = le_csv(cam_lhs)
    m = pontos_confiaveis(d)
    n_tot = len(d['status'])
    print(f'\n=== ESTUDO B: LHS ===')
    print(f'{n_tot} perfis amostrados, {m.sum()} utilizaveis '
          f'(avaliados e com estol confiavel)')
    if m.sum() < 30:
        print('AMOSTRA PEQUENA DEMAIS para concluir. Aumente n do LHS.')
        return 1

    cl = d['clmax'][m]
    bom = cl >= CLMAX_ALVO
    print(f'cl_max: min={cl.min():.3f} mediana={np.median(cl):.3f} '
          f'max={cl.max():.3f}')
    print(f'perfis com cl_max >= {CLMAX_ALVO} (alvo):       {bom.sum()}/{len(cl)} '
          f'({bom.mean()*100:.1f}%)')
    viavel = cl >= CLMAX_VIABILIDADE
    print(f'perfis com cl_max >= {CLMAX_VIABILIDADE} (viabilidade): '
          f'{viavel.sum()}/{len(cl)} ({viavel.mean()*100:.1f}%)')

    # 1. correlacao de cada descritor com cl_max
    print(f'\n--- 1. correlacao de posto com cl_max ---')
    print(f'{"descritor":14s} {"rho":>8s} {"p":>10s}')
    ranking = []
    for nome in dsc.NOMES:
        if nome not in d:
            continue
        v = d[nome][m]
        if np.std(v) == 0:
            continue
        rho, p = stats.spearmanr(v, cl)
        ranking.append((abs(rho), rho, p, nome))
    ranking.sort(reverse=True)
    for _, rho, p, nome in ranking:
        print(f'{nome:14s} {rho:8.3f} {p:10.2g}')

    # 2. limiar de um descritor so, nos dois niveis
    algum = False
    for alvo, etiqueta in [(CLMAX_ALVO, 'ALVO: mantem a aeronave B como esta'),
                           (CLMAX_VIABILIDADE,
                            'VIABILIDADE: limite de empuxo na decolagem')]:
        bom_a = cl >= alvo
        print(f'\n--- 2. limiar de UM descritor para cl_max >= {alvo} '
              f'({etiqueta}) ---')
        print(f'    {int(bom_a.sum())}/{len(cl)} perfis atendem '
              f'({bom_a.mean():.0%}); precisao exigida >= {PRECISAO_MIN:.0%}')
        if bom_a.sum() < 5 or bom_a.all():
            print('    (fracao degenerada -- limiar sem sentido neste nivel)')
            continue
        print(f'{"descritor":14s} {"limiar":>11s} {"precisao":>9s} {"cobertura":>10s}')
        for _, rho, _, nome in ranking:
            v = d[nome][m]
            q = qualidade_limiar(v, bom_a, maior_melhor=(rho > 0))
            if q is None:
                print(f'{nome:14s} {"--":>11s}   nenhum limiar atinge a precisao')
            else:
                t, prec, cob = q
                sinal = '>=' if rho > 0 else '<='
                util = cob >= COBERTURA_MIN
                marca = '  <-- utilizavel' if util else '  (corta demais)'
                if util and alvo == CLMAX_ALVO:
                    algum = True
                print(f'{nome:14s} {sinal}{t:10.4f} {prec:9.1%} {cob:10.1%}{marca}')

    # 3. superficie de resposta com os melhores descritores
    print(f'\n--- 3. superficie de resposta quadratica ---')
    for k in (1, 2, 3, 4):
        nomes_k = [n for _, _, _, n in ranking[:k]]
        X = np.column_stack([d[n][m] for n in nomes_k])
        r2, r2cv, *_ = superficie_quadratica(X, cl)
        print(f'  {k} descritor(es) {str(nomes_k):58s} '
              f'R2={r2:.3f}  R2_cv={r2cv:.3f}')

    print(f'\n--- conclusao ---')
    if algum:
        print('Existe limiar de um descritor so com precisao e cobertura')
        print('aceitaveis: use-o como batente, e o mais barato possivel.')
    else:
        print('Nenhum descritor isolado separa bem. Use a superficie de')
        print('resposta com o menor k cujo R2_cv ja sature -- ela continua')
        print('analitica, suave e de custo desprezivel dentro do otimizador.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
