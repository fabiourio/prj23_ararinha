'''
Lab 04 -- estudo de convergencia de malha do modelo AVL.

Refina simultaneamente o numero de paineis de todas as superficies do
fwd.avl e roda o AVL na condicao do ponto de projeto (M = 0,85 com meta
de CL = 0,5053) em cada nivel. As metricas acompanhadas sao o alpha de
equilibrio, o arrasto induzido, o fator de Oswald e o ponto neutro.

Gera relatorio_lab04/tables/convergencia_malha.csv. A figura do relatorio
sai desse CSV pelo script Julia lab04_convergencia_fig.jl.

Rodar da raiz do repo:  python lab04_convergencia.py
(requer avl/fwd.avl atualizado; rode lab04_gera_avl.py antes se preciso)
'''

# IMPORTS
import os
import re
import subprocess

import numpy as np

#=========================================

# Niveis de refino: (asa Nc, asa Ns, EH Nc, EH Ns, EV Nc, EV Ns), com
# totais de vortices proximos de 100, 200, 500, 1000, 1300 (malha da
# entrega) e 2000. Alem disso a varredura ja esta convergida e os desvios
# so flutuam no ruido. O ultimo nivel, quatro vezes mais fino que a malha
# adotada, nao entra na figura: serve de referencia para o calculo do erro.
NIVEIS = [
    (3, 12, 3, 5, 2, 4),
    (5, 16, 3, 7, 3, 5),
    (8, 25, 5, 10, 4, 8),
    (10, 35, 7, 14, 7, 10),
    (12, 40, 8, 16, 8, 12),
    (15, 50, 10, 20, 10, 15),
    (24, 80, 16, 32, 16, 24),
]
NIVEL_ADOTADO = 4                        # indice em NIVEIS (malha da entrega)

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

    # missing_ok: nao falha se o temporario ja tiver sido removido por fora
    try:
        os.remove(f'avl/saidas/conv_{i}.avl')
    except FileNotFoundError:
        pass
    return {'vortices': pega(r'(\d+)\s+Vortices'),
            'alpha': pega(r'Alpha =\s+([-\d.]+)'),
            'CLtot': pega(r'CLtot =\s+([-\d.]+)'),
            'CLa': pega(r'CLa =\s+([-\d.]+)'),
            'CDind': pega(r'CDind =\s*([-\d.Ee+]+)'),
            'CDff': pega(r'CDff\s*=\s*([-\d.Ee+]+)'),
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
malha_adotada = f'{NIVEIS[NIVEL_ADOTADO][0]}x{NIVEIS[NIVEL_ADOTADO][1]}'
for i, nivel in enumerate(NIVEIS):
    res = roda_nivel(i, nivel, base)
    res['malha'] = f'{nivel[0]}x{nivel[1]}'
    if np.isnan(res['CLtot']):
        print(f"nivel {i} (asa {res['malha']}): DESCARTADO "
              '(provavel limite de vortices do executavel)')
        continue
    resultados.append(res)
    print(f"nivel {i} (asa {res['malha']}): N = {res['vortices']:.0f}, "
          f"alpha = {res['alpha']:.4f}, CLa = {res['CLa']:.4f}, "
          f"CDind = {res['CDind']:.5f}, CDff = {res['CDff']:.5f}, "
          f"e = {res['e']:.4f}, xnp = {res['xnp']:.4f}")

# Tabela CSV
with open('relatorio_lab04/tables/convergencia_malha.csv', 'w',
          encoding='ascii') as f:
    f.write('malha_asa,vortices,alpha_deg,cla,CDind,CDff,e_oswald,xnp_m\n')
    for res in resultados:
        f.write(f"{res['malha']},{res['vortices']:.0f},{res['alpha']:.4f},"
                f"{res['CLa']:.4f},{res['CDind']:.6f},{res['CDff']:.6f},"
                f"{res['e']:.4f},{res['xnp']:.4f}\n")

# Desvios do nivel mais fino, para o texto do relatorio
fino = resultados[-1]
i_adot = next(i for i, res in enumerate(resultados)
              if res['malha'] == malha_adotada)
adotado = resultados[i_adot]
print(f"\nMalha mais fina que rodou: {fino['malha']} "
      f"({fino['vortices']:.0f} vortices)")
print('Desvios da malha adotada para a mais fina:')
for chave in ('alpha', 'CLa', 'CDind', 'CDff', 'xnp'):
    desvio = 100*abs(adotado[chave] - fino[chave])/abs(fino[chave])
    print(f'  {chave}: {desvio:.2f}%')
