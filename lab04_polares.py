'''
Lab 04, itens 5 a 7 -- polares da aeronave no ponto de projeto.

Para cada caso (CG dianteiro/traseiro, sem deflexao de profundor ou com
arfagem compensada), roda o AVL em M = 0,85 com a incidencia de empenagem
do item 3, varrendo metas de CL de -0,5 ate o CLmax do metodo da secao
critica daquele caso. O ponto de projeto (CL = 0,5053) entra exatamente na
grade.

Gera:
  - relatorio_lab04/tables/polar_<caso>.csv (alpha, CL, CD, delta_e)
  - avl/saidas/polares.json (CD no CL de projeto por caso)

Rodar da raiz do repo:  python lab04_polares.py
(requer trim.json e secao_critica.json em avl/saidas/)
'''

# IMPORTS
import json
import os

import numpy as np

from avl_batch import roda_avl, pega_todos

#=========================================

MACH = 0.85
CL_PROJ = 0.5053
CL_MIN = -0.5

CASOS = {'fwd_livre': ('fwd', 'd2 d2 0'),
         'fwd_trim': ('fwd', 'd2 pm 0'),
         'aft_livre': ('aft', 'd2 d2 0'),
         'aft_trim': ('aft', 'd2 pm 0')}


def varre_cl(cg, restricao, it, metas):
    '''Roda a sequencia de metas de CL num unico processo do AVL.'''
    cmds = (f'load {cg}.avl\noper\n'
            f'm\nmn {MACH}\n\n'
            f'de\n1 {it}\n\n'
            f'{restricao}\n')
    for cl in metas:
        cmds += f'a c {cl:.4f}\nx\n'
    cmds += '\nquit\n'
    saida = roda_avl(cmds, timeout=580)
    dados = []
    for bloco in saida.split('Vortex Lattice Output')[1:]:
        dados.append({
            'alpha': pega_todos(bloco, r'Alpha =\s+([-\d.]+)')[0],
            'CL': pega_todos(bloco, r'CLtot =\s+([-\d.]+)')[0],
            'CD': pega_todos(bloco, r'CDtot =\s+([-\d.]+)')[0],
            'CDind': pega_todos(bloco, r'CDind =\s*([-\d.Ee+]+)')[0],
            'de': pega_todos(bloco, r'elevator\s+=\s+([-\d.]+)')[0]})
    return dados


with open('avl/saidas/trim.json', encoding='ascii') as f:
    trim = json.load(f)
with open('avl/saidas/secao_critica.json', encoding='ascii') as f:
    critica = json.load(f)

os.makedirs('relatorio_lab04/tables', exist_ok=True)
resumo = {}
for caso, (cg, restricao) in CASOS.items():
    it = trim[cg]['it_deg']
    cl_max = critica[caso]['CLmax']
    metas = np.arange(CL_MIN, cl_max, 0.05)
    metas = np.unique(np.round(np.append(metas, [CL_PROJ, cl_max]), 4))

    dados = varre_cl(cg, restricao, it, metas)

    with open(f'relatorio_lab04/tables/polar_{caso}.csv', 'w',
              encoding='ascii') as f:
        f.write('alpha_deg,CL,CD,CDind,delta_e_deg\n')
        for d in dados:
            f.write(f"{d['alpha']:.4f},{d['CL']:.4f},{d['CD']:.5f},"
                    f"{d['CDind']:.5f},{d['de']:.4f}\n")

    ponto = min(dados, key=lambda d: abs(d['CL'] - CL_PROJ))
    resumo[caso] = {'CD_projeto': round(ponto['CD'], 5),
                    'alpha_projeto_deg': round(ponto['alpha'], 4),
                    'delta_e_projeto_deg': round(ponto['de'], 4),
                    'n_pontos': len(dados)}
    print(f"{caso}: {len(dados)} pontos, CD(CL=0,5053) = {ponto['CD']:.5f}, "
          f"alpha = {ponto['alpha']:.3f}, delta_e = {ponto['de']:.3f}")

with open('avl/saidas/polares.json', 'w', encoding='ascii') as f:
    json.dump(resumo, f, indent=2)
print('gravado: avl/saidas/polares.json')
