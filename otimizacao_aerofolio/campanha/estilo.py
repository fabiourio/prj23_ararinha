'''
Estilo das figuras do Lab 03, herdado do Lab 02 para o relatorio ficar
coerente. A paleta foi verificada para daltonismo: pior par adjacente com
dE 9.2 (deutan) e 27.6 (visao normal), acima dos pisos exigidos.

O verde PAL[2] tem contraste 2.74:1 contra o fundo, abaixo de 3:1 -- por isso
so aparece acompanhado de rotulo direto, nunca sozinho como unica pista.
'''

import matplotlib.pyplot as plt

PAL = ['#2a78d6', '#eb6834', '#1baf7a']
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


def titulo(ax, texto, subtitulo=None):
    ax.set_title(texto, color=INK, fontsize=11.5, fontweight='bold',
                 loc='left', pad=26 if subtitulo else 8)
    if subtitulo:
        ax.text(0.0, 1.018, subtitulo, transform=ax.transAxes,
                fontsize=9, color=INK2, ha='left', va='bottom')


def salvar(fig, caminho):
    fig.savefig(caminho, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('  figura:', caminho)


# rotulos legiveis para os descritores
ROTULO = {
    'r_LE_sup': r'raio de bordo de ataque  $r_{LE}/c$  (extradorso)',
    'r_LE_inf': r'raio de bordo de ataque  $r_{LE}/c$  (intradorso)',
    'delta_y': 'afiamento de bordo de ataque  $\\Delta y$  [% da corda]',
    't_01': 'espessura em $x/c = 1$%',
    't_05': 'espessura em $x/c = 5$%',
    't_max': r'espessura máxima  $(t/c)_{max}$',
    'x_tmax': r'posição da espessura máxima  $x_{t,max}/c$',
    't_min': 'espessura mínima',
    'c_max': r'arqueamento máximo  $(h/c)_{max}$',
    'x_cmax': r'posição do arqueamento máximo  $x_{h,max}/c$',
    'camber_085': r'arqueamento em $x/c = 0{,}85$  (carregamento traseiro)',
    'dcamber_te': 'inclinação da linha média no bordo de fuga',
}
