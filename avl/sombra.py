'''
Sombreamento, com geometria REAL do designTool (sem o deslocamento em Z do
AVL, cuja esteira e reta e nao enxerga o fenomeno).

Sombra da asa sobre a EH: faixa entre as retas que partem do bordo de
ataque e do bordo de fuga de cada secao da asa na direcao do escoamento
local, inclinada de gama = alfa - eps = (1 - dε/dα)·alfa em relacao ao eixo
do corpo (raiz sem incidencia, eps0 desprezado). A EH esta fora da sombra
quando o extradorso do seu bordo de ataque fica abaixo da reta do bordo de
fuga. Como gama cresce com alfa, a reta so sobe: depois de sair, a EH nao
volta.

Leme encoberto (recuperacao de parafuso): regiao entre a reta a 60 graus
pelo bordo de ataque da raiz da EH e a reta a 30 graus pelo bordo de fuga,
as mesmas linhas tracejadas do plots.py do designTool. So informativo.
'''

import numpy as np


def alfa_saida(av, n=41):
    '''alfa [graus] a partir do qual toda a EH fica abaixo da sombra da asa.'''
    w, h = av['asa'], av['EH']
    y = np.linspace(0.0, h['yt'], n)
    s_h = y/h['yt']
    x_le_h = h['xr'] + s_h*(h['xt'] - h['xr'])
    c_h = h['cr'] + s_h*(h['ct'] - h['cr'])
    z_topo = h['zr'] + s_h*(h['zt'] - h['zr']) + h['tc']*c_h/2

    s_w = y/w['yt']
    x_te_w = w['xr'] + s_w*(w['xt'] - w['xr']) + w['cr'] + s_w*(w['ct'] - w['cr'])
    z_w = w['zr'] + s_w*(w['zt'] - w['zr'])

    gama = np.arctan2(z_topo - z_w, x_le_h - x_te_w).max()
    return float(np.degrees(gama)/(1.0 - av['deda']))


def fracao_leme_encoberto(av, n=400, charneira=0.70):
    h, v = av['EH'], av['EV']
    z = np.linspace(v['zr'], v['zt'], n)
    s = (z - v['zr'])/(v['zt'] - v['zr'])
    x_le = v['xr'] + s*(v['xt'] - v['xr'])
    c = v['cr'] + s*(v['ct'] - v['cr'])
    x_ch, x_te = x_le + charneira*c, x_le + c

    dz = z - h['zr']
    lim1 = h['xr'] + dz/np.tan(np.radians(60.0))
    lim2 = h['xr'] + h['cr'] + dz/np.tan(np.radians(30.0))
    sobre = np.clip(np.minimum(x_te, lim2) - np.maximum(x_ch, lim1), 0.0, None)
    sobre[dz < 0] = 0.0
    return float(sobre.sum()/(x_te - x_ch).sum())
