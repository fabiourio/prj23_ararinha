'''
Figuras do DOE de cl_max.

  1. doe_clmax_corte.png     -- corte controlado: so o nariz varia
  2. doe_clmax_dispersao.png -- LHS: cl_max contra os 4 melhores descritores
  3. doe_clmax_limiar.png    -- fronteira de decisao do descritor escolhido

Rodar depois de doe_clmax_xfoil.py:  python figuras_doe_clmax.py
'''

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import descritores as dsc
from analise_doe_clmax import (RES, le_csv, pontos_confiaveis,
                               qualidade_limiar)
from doe_clmax_xfoil import CLMAX_ALVO, CLMAX_VIABILIDADE
from estilo import INK, INK2, MUTED, PAL, ROTULO, salvar, style_axes, titulo


def ranking_descritores(d, m, cl):
    r = []
    for nome in dsc.NOMES:
        if nome not in d:
            continue
        v = d[nome][m]
        if np.std(v) == 0:
            continue
        rho, _ = stats.spearmanr(v, cl)
        if np.isfinite(rho):
            r.append((abs(rho), rho, nome))
    r.sort(reverse=True)
    return r


def fig_corte():
    cam = os.path.join(RES, 'doe_clmax_corte.csv')
    if not os.path.isfile(cam):
        return
    c = le_csv(cam)
    m = pontos_confiaveis(c)
    x, y = c['delta_y'][m], c['clmax'][m]
    o = np.argsort(x)
    x, y = x[o], y[o]

    fig, ax = plt.subplots(figsize=(8, 5))
    style_axes(ax)

    # os dois niveis: alvo (mantem a aeronave B) e viabilidade (empuxo)
    ax.axhspan(CLMAX_VIABILIDADE, CLMAX_ALVO, color=PAL[1], alpha=0.06, zorder=1)
    ax.axhline(CLMAX_ALVO, color=INK2, linewidth=1.1, linestyle='--', zorder=2)
    ax.annotate(f'$c_{{\\ell,max}} = {CLMAX_ALVO}$ — assumido no designTool',
                (x.max(), CLMAX_ALVO), xytext=(-4, 6),
                textcoords='offset points', ha='right', fontsize=9, color=INK2)
    ax.axhline(CLMAX_VIABILIDADE, color=PAL[1], linewidth=1.2,
               linestyle=':', zorder=2)
    ax.annotate(f'$c_{{\\ell,max}} = {CLMAX_VIABILIDADE}$ — abaixo disso os '
                'motores não decolam em 2.900 m',
                (x.max(), CLMAX_VIABILIDADE), xytext=(-4, 6),
                textcoords='offset points', ha='right', fontsize=9,
                color=PAL[1])

    ax.plot(x, y, '-', color=PAL[0], linewidth=1.3, alpha=0.55, zorder=3)
    ax.plot(x, y, 'o', color=PAL[0], markersize=6, markeredgecolor='white',
            markeredgewidth=1.0, zorder=4)

    # onde a curva cruza o alvo
    if y.min() < CLMAX_ALVO < y.max():
        xc = float(np.interp(CLMAX_ALVO, y, x))
        ax.plot([xc], [CLMAX_ALVO], 'o', color=PAL[1], markersize=9,
                markeredgecolor='white', markeredgewidth=1.3, zorder=6)
        ax.annotate(f'$\\Delta y = {xc:.2f}\\%$ da corda',
                    (xc, CLMAX_ALVO), xytext=(10, -20),
                    textcoords='offset points', fontsize=9.5, color=INK,
                    fontweight='bold',
                    arrowprops=dict(arrowstyle='-', color=MUTED, lw=1))

    ax.set_xlabel(ROTULO['delta_y'], color=INK2, fontsize=10)
    ax.set_ylabel(r'$c_{\ell,max}$  (XFoil viscoso)', color=INK2, fontsize=10)
    titulo(ax, 'Afiar o bordo de ataque destrói a sustentação máxima',
           'corte controlado a partir do NACA 1411: só o nariz varia; '
           r'Re = $4{,}2 \times 10^7$ e M = 0,266 (decolagem, na MAC)')
    salvar(fig, os.path.join(RES, 'doe_clmax_corte.png'))


def _dispersao(ax, x, y, bom, rotulo_x, mostrar_legenda):
    style_axes(ax)
    ax.axhline(CLMAX_ALVO, color=INK2, linewidth=1.0, linestyle='--', zorder=2)
    ax.plot(x[~bom], y[~bom], 'o', color=PAL[1], markersize=4.2, alpha=0.75,
            markeredgecolor='white', markeredgewidth=0.5, zorder=3,
            label=f'$c_{{\\ell,max}} < {CLMAX_ALVO}$')
    ax.plot(x[bom], y[bom], 'o', color=PAL[0], markersize=4.2, alpha=0.75,
            markeredgecolor='white', markeredgewidth=0.5, zorder=4,
            label=f'$c_{{\\ell,max}} \\geq {CLMAX_ALVO}$')
    ax.set_xlabel(rotulo_x, color=INK2, fontsize=9)
    if mostrar_legenda:
        leg = ax.legend(loc='lower right', fontsize=8.5, frameon=False)
        for t in leg.get_texts():
            t.set_color(INK2)


def fig_dispersao(d, m, cl, rank):
    bom = cl >= CLMAX_ALVO
    top = rank[:4]
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5))
    for k, (_, rho, nome) in enumerate(top):
        ax = axs.flat[k]
        _dispersao(ax, d[nome][m], cl, bom,
                   ROTULO.get(nome, nome) + f'   ($\\rho = {rho:+.2f}$)',
                   mostrar_legenda=(k == 0))
        if k % 2 == 0:
            ax.set_ylabel(r'$c_{\ell,max}$', color=INK2, fontsize=9)
    fig.suptitle('Qual descritor geométrico prediz a sustentação máxima?',
                 color=INK, fontsize=12, fontweight='bold', x=0.09, ha='left')
    fig.text(0.09, 0.945,
             f'{int(m.sum())} perfis do hipercubo latino nas 8 variáveis CST, '
             r'com estol capturado; $\rho$ = correlação de posto de Spearman',
             fontsize=9, color=INK2, ha='left')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    salvar(fig, os.path.join(RES, 'doe_clmax_dispersao.png'))


def fig_limiar(d, m, cl, rank):
    bom = cl >= CLMAX_ALVO
    # melhor descritor que admite um limiar utilizavel
    escolhido = None
    for _, rho, nome in rank:
        q = qualidade_limiar(d[nome][m], bom, maior_melhor=(rho > 0))
        if q is not None:
            escolhido = (nome, rho, q)
            break
    if escolhido is None:
        print('  nenhum limiar de um descritor so atinge a precisao exigida;')
        print('  figura de fronteira nao gerada (use a superficie de resposta)')
        return

    nome, rho, (t, prec, cob) = escolhido
    x = d[nome][m]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    style_axes(ax)

    maior_melhor = rho > 0
    lo, hi = (t, x.max() * 1.02) if maior_melhor else (x.min() * 0.98, t)
    ax.axvspan(lo, hi, color=PAL[0], alpha=0.07, zorder=1)

    ax.axhline(CLMAX_ALVO, color=INK2, linewidth=1.1, linestyle='--', zorder=2)
    ax.axvline(t, color=PAL[2], linewidth=1.8, zorder=3)

    ax.plot(x[~bom], cl[~bom], 'o', color=PAL[1], markersize=4.5, alpha=0.75,
            markeredgecolor='white', markeredgewidth=0.5, zorder=4,
            label=f'$c_{{\\ell,max}} < {CLMAX_ALVO}$')
    ax.plot(x[bom], cl[bom], 'o', color=PAL[0], markersize=4.5, alpha=0.75,
            markeredgecolor='white', markeredgewidth=0.5, zorder=5,
            label=f'$c_{{\\ell,max}} \\geq {CLMAX_ALVO}$')

    sinal = r'\geq' if maior_melhor else r'\leq'
    ax.annotate(f'limiar da restrição:  ${sinal}\\ {t:.4g}$',
                (t, cl.min()), xytext=(8 if maior_melhor else -8, 4),
                textcoords='offset points',
                ha='left' if maior_melhor else 'right',
                fontsize=9.5, color=PAL[2], fontweight='bold')
    ax.annotate(f'região aceita:\nprecisão {prec:.0%}, cobertura {cob:.0%}',
                (0.5 * (lo + hi), cl.max()), xytext=(0, -6),
                textcoords='offset points', ha='center', va='top',
                fontsize=9, color=INK2)

    ax.set_xlabel(ROTULO.get(nome, nome), color=INK2, fontsize=10)
    ax.set_ylabel(r'$c_{\ell,max}$  (XFoil viscoso)', color=INK2, fontsize=10)
    # canto inferior esquerdo: o unico livre, ja que a regiao aceita fica a
    # direita e o rotulo do limiar ocupa a base perto da linha verde
    leg = ax.legend(loc='lower left', fontsize=9, frameon=False)
    for txt in leg.get_texts():
        txt.set_color(INK2)
    titulo(ax, 'A restrição substituta que o otimizador vai carregar',
           'o limiar é analítico nos coeficientes CST, logo suave e de custo '
           'desprezível — ao contrário do próprio XFoil')
    salvar(fig, os.path.join(RES, 'doe_clmax_limiar.png'))


def main():
    cam = os.path.join(RES, 'doe_clmax_lhs.csv')
    if not os.path.isfile(cam):
        print(f'ERRO: {cam} nao existe. Rode doe_clmax_xfoil.py antes.')
        return 1
    d = le_csv(cam)
    m = pontos_confiaveis(d)
    cl = d['clmax'][m]
    rank = ranking_descritores(d, m, cl)

    fig_corte()
    fig_dispersao(d, m, cl, rank)
    fig_limiar(d, m, cl, rank)
    return 0


if __name__ == '__main__':
    sys.exit(main())
