'''
Figuras compartilhadas do Lab 02 (equipe).

style_axes : grade/spines usados nos historicos e Pareto.
planform / sideview : desenho esquematico a partir do dict do designTool
                      (apos analyze). Nao substituem plot_geometry 3D.
'''

import numpy as np

PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100',
       '#e87ba4', '#008300', '#4a3aa7', '#e34948']
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#c3c2b7'
GRID = '#e1e0d9'


def style_axes(ax):
    ax.grid(color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9)


def planform(ax, ap, color, label):
    '''Vista de planta: asa, HT, fuselagem e trens.'''

    i = ap['inputs']
    g = ap['geometry']

    for s in (1, -1):
        ax.plot([i['xr_w'], g['xt_w'], g['xt_w'] + g['ct_w'], i['xr_w'] + g['cr_w'], i['xr_w']],
                s * np.array([0, g['yt_w'], g['yt_w'], 0, 0]), color=color, lw=1.5,
                label=label if s == 1 else None)
        ax.plot([g['xr_h'], g['xt_h'], g['xt_h'] + g['ct_h'], g['xr_h'] + g['cr_h'], g['xr_h']],
                s * np.array([0, g['yt_h'], g['yt_h'], 0, 0]), color=color, lw=1.5)
        ax.plot([0, i['L_f']], [s * i['D_f'] / 2, s * i['D_f'] / 2], color=color, lw=1.0)

    ax.plot([i['x_mlg']] * 2, [-i['y_mlg'], i['y_mlg']], 'o', color=color, ms=5)
    ax.plot([i['x_nlg']], [0], 'o', color=color, ms=5)


def sideview(ax, ap, color, label):
    '''Vista lateral: fuselagem, asa, EV e linha do solo (z_lg).'''

    i = ap['inputs']
    g = ap['geometry']

    th = np.linspace(0, 2 * np.pi, 120)
    ax.plot(i['L_f'] / 2 + i['L_f'] / 2 * np.cos(th), i['D_f'] / 2 * np.sin(th),
            color=color, lw=1.0, label=label)

    ax.plot([g['xr_v'], g['xt_v'], g['xt_v'] + g['ct_v'], g['xr_v'] + g['cr_v'], g['xr_v']],
            [i['zr_v'], g['zt_v'], g['zt_v'], i['zr_v'], i['zr_v']], color=color, lw=1.5)
    ax.plot([i['xr_w'], i['xr_w'] + g['cr_w']], [i['zr_w'], i['zr_w']], color=color, lw=2.5)
    ax.plot([g['xr_h'], g['xr_h'] + g['cr_h']], [i['zr_h'], i['zr_h']], color=color, lw=2.0)

    ax.plot([i['x_mlg'], i['x_nlg']], [i['z_lg'], i['z_lg']], 'o', color=color, ms=5)
    ax.plot([i['x_tailstrike']], [i['z_tailstrike']], 's', color=color, ms=5)
    ax.axhline(i['z_lg'], color=color, lw=0.6, ls=':')
