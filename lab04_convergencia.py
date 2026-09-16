'''
Lab 04 -- estudo de convergencia de malha do modelo AVL.

Refina simultaneamente o numero de paineis de todas as superficies do
fwd.avl e roda o AVL na condicao do ponto de projeto (M = 0,85 com meta
de CL = 0,5053) em cada nivel. As metricas acompanhadas sao o alpha de
equilibrio, o arrasto induzido, o fator de Oswald e o ponto neutro.

Gera:
  - relatorio_lab04/tables/convergencia_malha.csv
  - relatorio_lab04/images/03_convergencia/convergencia_malha.png

Rodar da raiz do repo:  python lab04_convergencia.py
(requer avl/fwd.avl atualizado; rode lab04_gera_avl.py antes se preciso)
'''

# IMPORTS
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.join('otimizacao_aerofolio', 'campanha'))
import matplotlib.pyplot as plt          # noqa: E402
from estilo import PAL, INK2, style_axes, salvar  # noqa: E402

#=========================================

# Niveis de refino: (asa Nc, asa Ns, EH Nc, EH Ns, EV Nc, EV Ns).
# O nivel 3 e a malha adotada nos arquivos da entrega.
NIVEIS = [
    (6, 20, 4, 8, 4, 6),
    (9, 30, 6, 12, 6, 9),
    (12, 40, 8, 16, 8, 12),
    (18, 60, 12, 24, 12, 18),
    (24, 80, 16, 32, 16, 24),
]
NIVEL_ADOTADO = 2                        # indice em NIVEIS (malha da entrega)

MALHA_BASE = {'wing': '12 1.0 40 1.0', 'ht': '8 1.0 16 1.0',
              'vt': '8 1.0 12 1.0'}

COMANDOS = ('load saidas/conv_{i}.avl\n'
            'oper\nm\nmn 0.85\n\n'
            'a c 0.5053\nx\nst\n\n\nquit\n')


def roda_nivel(i, nivel, base):
    '''Escreve a variante de malha, roda o AVL e extrai as metricas.'''
    wc, ws, hc, hs, vc, vs = nivel
    texto = base.replace(MALHA_BASE['wing'], f'{wc} 1.0 {ws} 1.0')
    texto = texto.replace(MALHA_BASE['ht'], f'{hc} 1.0 {hs} 1.0')
    texto = texto.replace(MALHA_BASE['vt'], f'{vc} 1.0 {vs} 1.0')
    with open(f'avl/saidas/conv_{i}.avl', 'w', encoding='ascii') as f:
        f.write(texto)

    r = subprocess.run([os.path.abspath(os.path.join('avl', 'avl.exe'))],
                       cwd='avl', input=COMANDOS.format(i=i),
                       capture_output=True, text=True, timeout=580)
    saida = r.stdout

    def pega(padrao):
        m = re.findall(padrao, saida)
        return float(m[-1]) if m else np.nan

    return {'vortices': pega(r'(\d+)\s+Vortices'),
            'alpha': pega(r'Alpha =\s+([-\d.]+)'),
            'CLtot': pega(r'CLtot =\s+([-\d.]+)'),
            'CDind': pega(r'CDind =\s*([-\d.Ee+]+)'),
            'e': pega(r'\se =\s+([-\d.]+)'),
            'xnp': pega(r'Xnp =\s+([-\d.]+)')}


# EXECUTION
with open('avl/fwd.avl', encoding='ascii') as f:
    base = f.read()
for chave, linha in MALHA_BASE.items():
    if linha not in base:
        raise RuntimeError(f'linha de malha "{linha}" nao achada no fwd.avl '
                           '-- rode lab04_gera_avl.py e confira MALHA_BASE')

os.makedirs('avl/saidas', exist_ok=True)
os.makedirs('relatorio_lab04/tables', exist_ok=True)
os.makedirs('relatorio_lab04/images/03_convergencia', exist_ok=True)

resultados = []
for i, nivel in enumerate(NIVEIS):
    res = roda_nivel(i, nivel, base)
    res['malha'] = f'{nivel[0]}x{nivel[1]}'
    resultados.append(res)
    print(f"nivel {i} (asa {res['malha']}): N = {res['vortices']:.0f}, "
          f"alpha = {res['alpha']:.4f}, CDind = {res['CDind']:.5f}, "
          f"e = {res['e']:.4f}, xnp = {res['xnp']:.4f}")

# Tabela CSV
with open('relatorio_lab04/tables/convergencia_malha.csv', 'w',
          encoding='ascii') as f:
    f.write('malha_asa,vortices,alpha_deg,CL,CDind,e_oswald,xnp_m\n')
    for res in resultados:
        f.write(f"{res['malha']},{res['vortices']:.0f},{res['alpha']:.4f},"
                f"{res['CLtot']:.4f},{res['CDind']:.6f},{res['e']:.4f},"
                f"{res['xnp']:.4f}\n")

# Figura: tres metricas contra o numero de vortices
NN = [res['vortices'] for res in resultados]
fig, eixos = plt.subplots(1, 3, figsize=(12, 3.6))
paineis = [('alpha', r'$\alpha$ para $C_L = 0{,}5053$  [graus]'),
           ('CDind', r'$C_{D,\mathrm{ind}}$'),
           ('xnp', r'$x_{np}$  [m]')]
for ax, (chave, rotulo) in zip(eixos, paineis):
    valores = [res[chave] for res in resultados]
    style_axes(ax)
    ax.plot(NN, valores, '-o', color=PAL[0], linewidth=2, markersize=6)
    ax.axvline(NN[NIVEL_ADOTADO], color=INK2, linewidth=0.8,
               linestyle=(0, (4, 3)))
    ax.set_xlabel('numero de vortices', fontsize=9, color=INK2)
    ax.set_title(rotulo, fontsize=10, loc='left')
eixos[0].annotate('malha adotada', xy=(NN[NIVEL_ADOTADO], eixos[0].get_ylim()[1]),
                  fontsize=8, color=INK2, ha='left', va='top',
                  xytext=(NN[NIVEL_ADOTADO]*1.08, eixos[0].get_ylim()[1]))
fig.tight_layout()
salvar(fig, 'relatorio_lab04/images/03_convergencia/convergencia_malha.png')

# Desvios do nivel mais fino, para o texto do relatorio
fino = resultados[-1]
adotado = resultados[NIVEL_ADOTADO]
print('\nDesvios da malha adotada para a mais fina:')
for chave in ('alpha', 'CDind', 'xnp'):
    desvio = 100*abs(adotado[chave] - fino[chave])/abs(fino[chave])
    print(f'  {chave}: {desvio:.2f}%')
