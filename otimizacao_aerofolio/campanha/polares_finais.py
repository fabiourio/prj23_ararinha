'''
Polares do perfil otimizado -- itens 7 e 8 do roteiro.

  item 7: polar de arrasto TRANSONICA (Euler) do perfil otimizado, varrendo
          alpha de -2 a +2 graus em torno do otimo, sobreposta ao RAE2822
          -- um perfil supercritico de verdade, a referencia justa.

  item 8: curvas SUBSONICAS no XFoil viscoso (cl x alpha, cl x cd, cm x alpha)
          do perfil de partida e do otimizado, no Reynolds de decolagem.

O item 7 responde se a polar do otimizado tem alguma peculiaridade. Espera-se
que tenha: otimizacao transonica mono-ponto e conhecida por produzir um poco
de arrasto estreito -- excelente exatamente no ponto de projeto e pior
logo ao lado. E o que a rodada multiponto (etapa seguinte da campanha)
deveria corrigir.

Rodar com:  python polares_finais.py [estacao] [--transonica] [--subsonica]
'''

import os
import sys
import time
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import otimiza_secao as ot
import xfoil_runner as xr
from analisa_otimos import RES, carrega, roda_euler
from estilo import INK, INK2, MUTED, PAL, salvar, style_axes, titulo

# RAE2822 -- perfil supercritico classico, coeficientes do single_run.py
AU_RAE = [0.12569744, 0.14564925, 0.1520378, 0.21346977, 0.17867926,
          0.20872853]
AL_RAE = [-0.13334359, -0.11827498, -0.22179445, -0.12644608, -0.08233977,
          0.05193736]

DALPHA = np.arange(-2.0, 2.01, 0.5)      # graus, em torno do otimo

# Reynolds de decolagem por estacao (ver documento de projeto, secao 4)
RE_DECOLAGEM = {'raiz': 6.06e7, 'meio': 4.23e7, 'ponta': 1.45e7}


def _um_ponto(args):
    Al, Au, alpha, rotulo = args
    try:
        r = roda_euler(Al, Au, alpha, 1.0)
        return rotulo, alpha, r['CL'], r['CD'], r['CM']
    except Exception:                                          # noqa: BLE001
        return rotulo, alpha, np.nan, np.nan, np.nan


def polar_transonica(nome, com):
    Alo, Auo, a_ot = ot.desmonta(com['xx_otimo'])
    alphas = a_ot + np.radians(DALPHA)

    tarefas = [(Alo, Auo, a, 'otimizado') for a in alphas]
    tarefas += [(np.array(AL_RAE), np.array(AU_RAE), a, 'RAE2822')
                for a in alphas]

    print(f'  {len(tarefas)} avaliacoes do Euler...', flush=True)
    t0 = time.time()
    with Pool(4) as pool:
        saidas = pool.map(_um_ponto, tarefas)
    print(f'  {time.time()-t0:.0f} s')

    dados = {}
    for rot, a, cl, cd, cm in saidas:
        dados.setdefault(rot, []).append((a, cl, cd, cm))
    for rot in dados:
        dados[rot] = np.array(sorted(dados[rot]))

    fig, ax = plt.subplots(figsize=(8, 6))
    style_axes(ax)
    for rot, cor, marcador in [('otimizado', PAL[0], 'o'),
                               ('RAE2822', PAL[1], 's')]:
        d = dados[rot]
        m = np.isfinite(d[:, 2])
        ax.plot(d[m, 2], d[m, 1], '-', color=cor, linewidth=1.4, alpha=0.6)
        ax.plot(d[m, 2], d[m, 1], marcador, color=cor, markersize=5,
                markeredgecolor='white', markeredgewidth=0.8, label=rot)

    # destaca o ponto de projeto
    d = dados['otimizado']
    i0 = int(np.argmin(np.abs(d[:, 0] - a_ot)))
    ax.plot([d[i0, 2]], [d[i0, 1]], '*', color=PAL[2], markersize=18,
            markeredgecolor=INK, markeredgewidth=1.0, zorder=6)
    ax.annotate('ponto de projeto', (d[i0, 2], d[i0, 1]), xytext=(12, -4),
                textcoords='offset points', fontsize=9.5, color=INK,
                fontweight='bold')

    ax.set_xlabel(r'$c_d$', color=INK2, fontsize=10)
    ax.set_ylabel(r'$c_\ell$', color=INK2, fontsize=10)
    leg = ax.legend(fontsize=9, frameon=False, loc='lower right')
    for t in leg.get_texts():
        t.set_color(INK2)
    titulo(ax, f'Polar transônica — estação {nome}',
           f'$M_n$ = {ot.MACH_N}, varredura de ±2° em torno do ótimo; '
           'RAE2822 como referência supercrítica')
    salvar(fig, os.path.join(RES, f'polar_transonica_{nome}.png'))

    cam = os.path.join(RES, f'polar_transonica_{nome}.csv')
    with open(cam, 'w') as fid:
        fid.write('perfil,alpha_rad,alpha_deg,CL,CD,CM\n')
        for rot, d in dados.items():
            for a, cl, cd, cm in d:
                fid.write(f'{rot},{a:.6f},{np.degrees(a):.4f},'
                          f'{cl:.6f},{cd:.6f},{cm:.6f}\n')
    print(f'  dados: {cam}')


def curvas_subsonicas(nome, com):
    Ali, Aui, _ = ot.desmonta(com['xx_ini'])
    Alo, Auo, _ = ot.desmonta(com['xx_otimo'])
    Re = RE_DECOLAGEM[nome]

    print(f'  XFoil viscoso em Re = {Re:.2e} (decolagem)...', flush=True)
    curvas = {}
    for rot, (Au, Al) in [('partida', (Aui, Ali)), ('otimizado', (Auo, Alo))]:
        # ate 28 graus: a 20 a varredura parava ANTES do estol e o cl_max
        # saia como limite inferior ('nao_atingido'), escondendo justamente a
        # comparacao de sustentacao maxima que o item 8 pede
        r = xr.clmax_cst(Au, Al, Re=Re, Mach=0.0, alpha_seq=(-6.0, 28.0, 0.5),
                         timeout=900)
        curvas[rot] = r['polar']
        print(f'    {rot}: {r["n_pontos"]} pontos, '
              f'cl_max = {r["clmax"]:.4f} ({r["situacao"]})')

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.5))
    eixos = [('alpha', 'CL', r'$\alpha$ [graus]', r'$c_\ell$'),
             ('CD', 'CL', r'$c_d$', r'$c_\ell$'),
             ('alpha', 'CM', r'$\alpha$ [graus]', r'$c_m$')]
    for ax, (kx, ky, rx, ry) in zip(axs, eixos):
        style_axes(ax)
        for rot, cor in [('partida', MUTED), ('otimizado', PAL[0])]:
            p = curvas[rot]
            if len(p['CL']) == 0:
                continue
            ax.plot(p[kx], p[ky], '-', color=cor, linewidth=1.7, label=rot)
        ax.set_xlabel(rx, color=INK2, fontsize=10)
        ax.set_ylabel(ry, color=INK2, fontsize=10)
    leg = axs[0].legend(fontsize=9, frameon=False, loc='lower right')
    for t in leg.get_texts():
        t.set_color(INK2)

    fig.suptitle(f'Comportamento subsônico — estação {nome}',
                 color=INK, fontsize=12, fontweight='bold', x=0.06, ha='left',
                 y=1.06)
    fig.text(0.06, 0.99,
             f'XFoil viscoso, Re = {Re:.2e} (condição de decolagem), '
             'Mach desprezado como permite o roteiro',
             fontsize=9, color=INK2, ha='left', va='top')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    salvar(fig, os.path.join(RES, f'subsonico_{nome}.png'))


def main():
    quais = [a for a in sys.argv[1:] if not a.startswith('--')] or ['meio']
    fazer_t = '--transonica' in sys.argv or not any(
        a.startswith('--') for a in sys.argv[1:])
    fazer_s = '--subsonica' in sys.argv or not any(
        a.startswith('--') for a in sys.argv[1:])

    for nome in quais:
        com = carrega(f'otim_{nome}')
        if com is None:
            print(f'{nome}: sem resultado ainda')
            continue
        print(f'--- {nome} ---')
        if fazer_s:
            curvas_subsonicas(nome, com)
        if fazer_t:
            polar_transonica(nome, com)


if __name__ == '__main__':
    main()
