'''
Figuras do relatorio -- desenha a partir dos CSVs, sem tocar em solver.

As curvas vem de curvas_relatorio.py, que roda o Euler e o XFoil uma vez e
grava tudo em relatorio/tables/. Aqui so se le CSV, entao replotar custa
segundos e da para iterar no estilo a vontade.

Gera em relatorio/images/:

  04_otimizacao/historico_convergencia.png   item 4
  05_transonico/geometria.png                item 5 -- os dois perfis da Tab. 2
  05_transonico/cp.png                       item 5
  05_transonico/mach.png                     item 5
  05_transonico/geometria_item9.png          item 9 -- otimizado x cusp
  06_polar/polar_transonica.png              item 7
  07_subsonico/cl_alpha.png                  item 8
  07_subsonico/cl_cd.png                     item 8
  07_subsonico/cm_alpha.png                  item 8

ESCALA: as figuras sao desenhadas na largura final em polegadas, com corpo 9,
de modo que \\includegraphics[width=0.8\\textwidth] as coloque na pagina em
1:1 -- sem reescala, e com o texto da figura do mesmo tamanho do texto do
relatorio. Nao ha titulo dentro da figura: quem titula e a legenda do LaTeX.

Rodar com:  python figuras_relatorio.py
'''

import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                # noqa: E402
import numpy as np                                             # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RELATORIO = os.path.normpath(os.path.join(AQUI, '..', '..', 'relatorio'))
TABLES = os.path.join(RELATORIO, 'tables')
IMG = os.path.join(RELATORIO, 'images')

# largura util do relatorio (a4, margem de 1in) x a fracao usada no
# \includegraphics -- ver ESCALA no cabecalho
TEXTWIDTH_IN = 6.5
L80 = 0.80 * TEXTWIDTH_IN
L90 = 0.90 * TEXTWIDTH_IN

CL_REF, TC_REF, MACH_N = 0.7319, 0.1772, 0.7196
LIMIAR_BLUNTEZ = 0.0938
RE_DECOLAGEM = 4.23e7

# azul e laranja da paleta da equipe (verificada para daltonismo); o cinza
# carrega o perfil de partida, que e referencia e nao protagonista
C_PART, C_OPT, C_CUSP, C_REF = '#8c8c8c', '#0072B2', '#D55E00', '#b0b0b0'

ROT = {'partida': 'NACA 1411 reescalado', 'otimizado': 'otimizado',
       'cusp_liberado': 'cusp liberado'}

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'font.size': 9,
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'legend.fontsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'axes.linewidth': 0.6,
    'lines.linewidth': 1.3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'legend.frameon': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
})


def le(arq):
    '''CSV longo -> {chave: {coluna: array}}, agrupado pela 1a coluna.'''
    cam = os.path.join(TABLES, arq)
    if not os.path.isfile(cam):
        print(f'  AUSENTE: {arq} -- rode curvas_relatorio.py')
        return None
    bruto = defaultdict(lambda: defaultdict(list))
    with open(cam, encoding='utf-8') as fid:
        r = csv.reader(fid)
        cab = next(r)
        for lin in r:
            for col, v in zip(cab[1:], lin[1:]):
                bruto[lin[0]][col].append(float(v))
    return {k: {c: np.array(v) for c, v in d.items()} for k, d in bruto.items()}


def eixo(ax):
    for lado in ('top', 'right'):
        ax.spines[lado].set_visible(False)


def salva(fig, sub, arq):
    destino = os.path.join(IMG, sub)
    os.makedirs(destino, exist_ok=True)
    fig.savefig(os.path.join(destino, arq), facecolor='white')
    plt.close(fig)
    print(f'  {sub}/{arq}')


# ---------------------------------------------------------------------------
def fig_geometria(d):
    for arq, a, b, ca, cb in (
            ('geometria.png', 'partida', 'otimizado', C_PART, C_OPT),
            ('geometria_item9.png', 'otimizado', 'cusp_liberado', C_OPT, C_CUSP)):
        fig, ax = plt.subplots(figsize=(L80, L80 * 0.30))
        ax.plot(d[a]['x'], d[a]['y'], '--', color=ca, label=ROT[a])
        ax.plot(d[b]['x'], d[b]['y'], '-', color=cb, label=ROT[b])
        ax.set_xlabel('$x/c$')
        ax.set_ylabel('$y/c$')
        ax.set_aspect('equal')
        ax.set_xlim(-0.02, 1.02)
        eixo(ax)
        ax.legend(loc='lower right', handlelength=2.4)
        salva(fig, '05_transonico', arq)


def fig_cp_mach(d):
    lab = {k: ROT[k] for k in d}
    fig, ax = plt.subplots(figsize=(L80, L80 * 0.62))
    ax.plot(d['partida']['x'], d['partida']['cp'], '--', color=C_PART,
            label=lab['partida'])
    ax.plot(d['otimizado']['x'], d['otimizado']['cp'], '-', color=C_OPT,
            label=lab['otimizado'])
    ax.invert_yaxis()
    ax.set_xlabel('$x/c$')
    ax.set_ylabel('$C_p$')
    ax.set_xlim(0, 1)
    eixo(ax)
    ax.legend(loc='lower right', handlelength=2.4)
    salva(fig, '05_transonico', 'cp.png')

    fig, ax = plt.subplots(figsize=(L80, L80 * 0.62))
    ax.axhline(1.0, color=C_REF, lw=0.7, ls=':', zorder=1)
    ax.text(0.995, 1.012, 'sônico', ha='right', va='bottom', fontsize=7.5,
            color='#707070')
    ax.plot(d['partida']['x'], d['partida']['mach'], '--', color=C_PART,
            label=lab['partida'], zorder=2)
    ax.plot(d['otimizado']['x'], d['otimizado']['mach'], '-', color=C_OPT,
            label=lab['otimizado'], zorder=3)
    ax.set_xlabel('$x/c$')
    ax.set_ylabel('Mach local')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.65)
    eixo(ax)
    ax.legend(loc='lower right', handlelength=2.4)
    salva(fig, '05_transonico', 'mach.png')


def fig_polar(d):
    fig, ax = plt.subplots(figsize=(L80, L80 * 0.72))
    for chave, cor, marca, rot in (('RAE2822', C_PART, 's', 'RAE 2822'),
                                   ('otimizado', C_OPT, 'o', 'otimizado')):
        if chave not in d:
            continue
        p = d[chave]
        o = np.argsort(p['alpha_deg'])
        ax.plot(p['cd'][o], p['cl'][o], '-', color=cor, marker=marca,
                markersize=3, markerfacecolor='white', markeredgewidth=0.8,
                label=rot)
    # ponto de projeto: o c_l da otimizacao
    p = d['otimizado']
    i = int(np.argmin(np.abs(p['cl'] - CL_REF)))
    ax.plot(p['cd'][i], p['cl'][i], '*', color=C_CUSP, markersize=11,
            markeredgecolor='white', markeredgewidth=0.6, zorder=5,
            label='ponto de projeto')
    ax.set_xlabel('$c_d$')
    ax.set_ylabel('$c_\\ell$')
    ax.set_xlim(0, None)
    eixo(ax)
    ax.legend(loc='lower right', handlelength=2.4)
    salva(fig, '06_polar', 'polar_transonica.png')


def fig_subsonico(d):
    eixos = (('alpha_deg', 'cl', r'$\alpha$ [graus]', '$c_\\ell$',
              'cl_alpha.png', None, None),
             ('cd', 'cl', '$c_d$', '$c_\\ell$', 'cl_cd.png', (0, 0.035), None),
             ('alpha_deg', 'cm', r'$\alpha$ [graus]', '$c_m$',
              'cm_alpha.png', None, None))
    for kx, ky, rx, ry, arq, xlim, ylim in eixos:
        fig, ax = plt.subplots(figsize=(L80, L80 * 0.62))
        for chave, cor, ls in (('partida', C_PART, '--'),
                               ('otimizado', C_OPT, '-')):
            p = d[chave]
            ax.plot(p[kx], p[ky], ls, color=cor, label=ROT[chave])
        ax.set_xlabel(rx)
        ax.set_ylabel(ry)
        if xlim:
            ax.set_xlim(*xlim)
        if ylim:
            ax.set_ylim(*ylim)
        if kx == 'alpha_deg':
            ax.set_xlim(-6, 26)
        eixo(ax)
        ax.legend(loc='lower right' if ky == 'cm' else 'upper left',
                  handlelength=2.4)
        salva(fig, '07_subsonico', arq)


def fig_convergencia(d):
    p = d['todos'] if 'todos' in d else None
    if p is None:                       # CSV de convergencia nao e agrupado
        return
    n = np.arange(1, len(p['cd']) + 1)
    # o otimo e o ponto VIAVEL de menor c_d. Sem o filtro, o argmin cai numa
    # avaliacao de busca em linha que tem c_d baixo so porque nao sustenta --
    # a avaliacao 4, por exemplo, da c_d = 0,0114 com c_l = 0,451.
    viavel = ((np.abs(p['cl'] - CL_REF) < 5e-3)
              & (p['tc_max'] >= TC_REF - 1e-4))
    i_ot = int(np.argmin(np.where(viavel, p['cd'], np.inf)))
    paineis = (('cd', '$c_d$', None, True),
               ('cl', '$c_\\ell$', CL_REF, False),
               ('tc_max', '$(t/c)_{\\max}$', TC_REF, False),
               ('bluntez', '$t_{01}/\\sqrt{t/c}$', LIMIAR_BLUNTEZ, False))

    fig, axs = plt.subplots(2, 2, figsize=(L90, L90 * 0.62))
    for ax, (col, rot, alvo, log) in zip(axs.ravel(), paineis):
        ax.plot(n, p[col], '-', color=C_OPT, marker='o', markersize=2.4,
                markerfacecolor='white', markeredgewidth=0.6)
        if alvo is not None:
            ax.axhline(alvo, color=C_PART, lw=0.7, ls=':')
        ax.plot(n[i_ot], p[col][i_ot], 'o', color=C_CUSP, markersize=4.5,
                markeredgecolor='white', markeredgewidth=0.6, zorder=5)
        if log:
            ax.set_yscale('log')
        ax.set_ylabel(rot)
        eixo(ax)
    for ax in axs[1]:
        ax.set_xlabel('avaliação do Euler')
    fig.align_ylabels(axs)
    fig.tight_layout(pad=0.4, w_pad=1.4, h_pad=0.8)
    salva(fig, '04_otimizacao', 'historico_convergencia.png')


def le_convergencia():
    cam = os.path.join(TABLES, 'curva_convergencia.csv')
    if not os.path.isfile(cam):
        print('  AUSENTE: curva_convergencia.csv -- rode curvas_relatorio.py')
        return None
    col = defaultdict(list)
    with open(cam, encoding='utf-8') as fid:
        for lin in csv.DictReader(fid):
            for k, v in lin.items():
                col[k].append(float(v))
    return {'todos': {k: np.array(v) for k, v in col.items()}}


def main():
    print(f'desenhando em {IMG}')
    g = le('curva_geometria.csv')
    if g:
        fig_geometria(g)
    cm = le('curva_cp_mach.csv')
    if cm:
        fig_cp_mach(cm)
    pol = le('tab4_polar_transonica.csv')
    if pol:
        fig_polar(pol)
    sub = le('curva_subsonico.csv')
    if sub:
        fig_subsonico(sub)
    cv = le_convergencia()
    if cv:
        fig_convergencia(cv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
