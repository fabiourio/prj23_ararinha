'''Estilo comum das figuras: marcas finas, grade discreta, paleta fixa.'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

AZUL, LARANJA, VERDE = '#2a78d6', '#eb6834', '#1baf7a'
TINTA, TINTA2, GRADE, CINZA = '#0b0b0b', '#52514e', '#e1e0d9', '#b9b8b3'

plt.rcParams.update({
    'figure.dpi': 150, 'savefig.dpi': 200, 'savefig.bbox': 'tight',
    'axes.grid': True, 'grid.color': GRADE, 'grid.linewidth': 0.7,
    'axes.edgecolor': CINZA, 'axes.labelcolor': TINTA2,
    'xtick.color': TINTA2, 'ytick.color': TINTA2, 'text.color': TINTA,
    'axes.spines.top': False, 'axes.spines.right': False,
    'lines.linewidth': 2.0, 'lines.markersize': 5, 'font.size': 9,
    'legend.frameon': False,
})
