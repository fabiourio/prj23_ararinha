'''
Figuras das otimizacoes de secao -- itens 4, 5 e 6 do roteiro.

  1. otim_<estacao>_convergencia.png -- historico de CD, CL, t/c e bluntez
  2. otim_<estacao>_geometria.png    -- perfis sobrepostos, partida x otimo
  3. otim_<estacao>_cp_mach.png      -- Cp e Mach ao longo da corda

As duas ultimas exigem reavaliar o Euler nos dois pontos (~1,5 min cada).

Rodar depois de otimiza_secao.py:  python figuras_otimizacao.py [estacao]
'''

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import otimiza_secao as ot
import xfoil_runner as xr
from analisa_otimos import RES, carrega, roda_euler
from estilo import INK, INK2, MUTED, PAL, salvar, style_axes, titulo


def fig_convergencia(nome, com, sem=None):
    h = com['hist']
    n = np.arange(1, len(h['CD']) + 1)

    fig, axs = plt.subplots(2, 2, figsize=(11, 7))
    paineis = [
        ('CD', r'$c_d$', None, None),
        ('CL', r'$c_\ell$', com['cl_ref'], 'alvo'),
        ('maxt', r'$(t/c)_{max}$', com['tc_ref'], 'mínimo'),
        ('bluntez', r'$t_{01}/\sqrt{t/c}$', ot.limiar_bluntez(nome), 'limiar'),
    ]
    for ax, (chave, rot, ref, nome_ref) in zip(axs.flat, paineis):
        style_axes(ax)
        if chave == 'bluntez':
            y = np.array([ot.bluntez(t, m)
                          for t, m in zip(h['t01'], h['maxt'])])
        else:
            y = np.array(h[chave])
        ax.plot(n, y, '-', color=PAL[0], linewidth=1.2, alpha=0.6)
        ax.plot(n, y, 'o', color=PAL[0], markersize=3, zorder=3,
                label='com restrição de $c_{\\ell,max}$' if sem is not None else None)

        if sem is not None:
            hs = sem['hist']
            ns = np.arange(1, len(hs['CD']) + 1)
            ys = (np.array([ot.bluntez(t, m)
                            for t, m in zip(hs['t01'], hs['maxt'])])
                  if chave == 'bluntez' else np.array(hs[chave]))
            ax.plot(ns, ys, '-', color=PAL[1], linewidth=1.2, alpha=0.5)
            ax.plot(ns, ys, 's', color=PAL[1], markersize=3, zorder=2,
                    label='sem a restrição')

        if ref is not None:
            ax.axhline(ref, color=INK2, linestyle='--', linewidth=1.0)
            ax.annotate(nome_ref, (n[-1], ref), xytext=(-3, 4),
                        textcoords='offset points', ha='right',
                        fontsize=8.5, color=INK2)

        ax.plot([n[com['i_otimo']] + 1], [y[com['i_otimo']]], 'o',
                color=PAL[2], markersize=9, markeredgecolor='white',
                markeredgewidth=1.4, zorder=6)
        ax.set_ylabel(rot, color=INK2, fontsize=10)
        ax.set_xlabel('avaliação do Euler', color=INK2, fontsize=9)
        if chave == 'CD':
            ax.set_yscale('log')

    if sem is not None:
        leg = axs.flat[0].legend(fontsize=8.5, frameon=False, loc='upper right')
        for t in leg.get_texts():
            t.set_color(INK2)

    cd0, cdf = com['hist']['CD'][0], com['hist']['CD'][com['i_otimo']]
    fig.suptitle(f'Convergência da otimização — estação {nome}',
                 color=INK, fontsize=12, fontweight='bold', x=0.07, ha='left')
    fig.text(0.07, 0.945,
             f'{com["n_aval"]} avaliações em {com["tempo_min"]:.0f} min; '
             f'$c_d$ de {cd0:.5f} para {cdf:.5f} ({(cdf/cd0-1)*100:+.1f}%); '
             f'ponto verde = ótimo',
             fontsize=9, color=INK2, ha='left')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    salvar(fig, os.path.join(RES, f'otim_{nome}_convergencia.png'))


def fig_geometria(nome, com):
    Ali, Aui, _ = ot.desmonta(com['xx_ini'])
    Alo, Auo, _ = ot.desmonta(com['xx_otimo'])
    fi = xr.cst_coords(Aui, Ali)
    fo = xr.cst_coords(Auo, Alo)

    fig, ax = plt.subplots(figsize=(10, 3.6))
    style_axes(ax)
    ax.plot(np.real(fi['x_coord']), np.real(fi['y_coord']), '-',
            color=MUTED, linewidth=1.8, label='partida (NACA 1411 reescalado)')
    ax.plot(np.real(fo['x_coord']), np.real(fo['y_coord']), '-',
            color=PAL[0], linewidth=2.0, label='otimizado')
    ax.axis('equal')
    ax.set_xlabel('$x/c$', color=INK2, fontsize=10)
    ax.set_ylabel('$y/c$', color=INK2, fontsize=10)
    leg = ax.legend(fontsize=9, frameon=False, loc='upper right')
    for t in leg.get_texts():
        t.set_color(INK2)
    titulo(ax, f'Geometria — estação {nome}',
           f'mesma espessura máxima ({com["tc_ref"]:.4f}); '
           'o otimizador redistribui, não engorda')
    salvar(fig, os.path.join(RES, f'otim_{nome}_geometria.png'))


def fig_cp_mach(nome, com):
    print('  reavaliando o Euler nos dois pontos...', flush=True)
    Ali, Aui, ai = ot.desmonta(com['xx_ini'])
    Alo, Auo, ao = ot.desmonta(com['xx_otimo'])
    ri = roda_euler(Ali, Aui, ai, 1.0)
    ro = roda_euler(Alo, Auo, ao, 1.0)

    fig, axs = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    for ax, chave, rot in [(axs[0], 'Cp', '$C_p$'),
                           (axs[1], 'Mach', 'Mach local')]:
        style_axes(ax)
        ax.plot(ri['x'], ri[chave], '-', color=MUTED, linewidth=1.6,
                label=f'partida  ($c_d$ = {ri["CD"]:.5f})')
        ax.plot(ro['x'], ro[chave], '-', color=PAL[0], linewidth=1.8,
                label=f'otimizado  ($c_d$ = {ro["CD"]:.5f})')
        ax.set_ylabel(rot, color=INK2, fontsize=10)
        if chave == 'Cp':
            ax.invert_yaxis()
        else:
            ax.axhline(1.0, color=PAL[1], linestyle='--', linewidth=1.1)
            ax.annotate('sônico (M = 1)', (0.98, 1.0), xytext=(-4, 5),
                        textcoords='offset points', ha='right',
                        fontsize=9, color=PAL[1])
        leg = ax.legend(fontsize=9, frameon=False)
        for t in leg.get_texts():
            t.set_color(INK2)
    axs[1].set_xlabel('$x/c$', color=INK2, fontsize=10)

    titulo(axs[0], f'Como o otimizador reduziu o arrasto — estação {nome}',
           f'$M_n$ = {ot.MACH_N}, $c_\\ell$ = {com["cl_ref"]:.4f}; '
           'a região acima de M = 1 é onde mora o choque')
    salvar(fig, os.path.join(RES, f'otim_{nome}_cp_mach.png'))


def main():
    quais = [a for a in sys.argv[1:] if not a.startswith('--')] or \
        list(ot.ESTACOES)
    for nome in quais:
        com = carrega(f'otim_{nome}')
        if com is None:
            print(f'{nome}: sem resultado ainda')
            continue
        sem = carrega(f'otim_{nome}_sem_bluntez')
        print(f'--- {nome} ---')
        fig_convergencia(nome, com, sem)
        fig_geometria(nome, com)
        fig_cp_mach(nome, com)


if __name__ == '__main__':
    main()
