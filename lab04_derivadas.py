'''
Lab 04, secao 3 -- derivadas de estabilidade para o CG traseiro.

Roda o aft.avl no ponto de projeto (M = 0,85, CL = 0,5053, it do item 3,
profundor fixo em zero, sem restricao de trimagem) e extrai as derivadas
dos comandos ft, st e sb, conforme o roteiro:
  - longitudinais e CY do st (eixos de estabilidade);
  - momentos latero-direcionais (Clp, Clr, Cnp, Cnr e controles) do sb;
  - CL0 e CM0 de uma rodada em alpha = 0 com it = 0;
  - CD0, CDa e CDa2 do ajuste parabolico da polar nao trimada do CG
    traseiro (item 5.b);
  - CDq, CDit e CDde por diferencas finitas centradas em alpha fixo.

As conversoes para o MVO seguem o enunciado: derivadas de controle (it,
de, da, dr) multiplicadas por 180/pi, e sinais invertidos em CYb, CYp,
CYr, Clda, Cldr, Cnda e Cndr.

Gera:
  - relatorio_lab04/tables/derivadas.csv (valor do AVL e valor do MVO)
  - avl/saidas/derivadas.json

Rodar da raiz do repo:  python lab04_derivadas.py
(requer trim.json e polar_aft_livre.csv gerados pelos scripts anteriores)
'''

# IMPORTS
import json
import os
import re

import numpy as np

from avl_batch import roda_avl, pega, pega_todos

#=========================================

MACH = 0.85
CL_PROJ = 0.5053
DQ = 0.05                                # perturbacao em qc/2V
DIT = 0.5                                # perturbacao em it [graus]
DDE = 1.0                                # perturbacao em delta_e [graus]


def tokens(bloco):
    '''Dicionario token -> valor de um bloco st ou sb.'''
    return {m.group(1): float(m.group(2))
            for m in re.finditer(r'(\w+)\s*=\s*([-+\dEe.]+)', bloco)}


with open('avl/saidas/trim.json', encoding='ascii') as f:
    it_aft = json.load(f)['aft']['it_deg']

# Rodada A: alpha = 0, it = 0, profundor 0 -> CL0 e CM0
saida_a = roda_avl('load aft.avl\noper\n'
                   f'm\nmn {MACH}\n\n'
                   'd2 d2 0\nde\n1 0.0\n\n'
                   'a a 0\nx\n\nquit\n')
CL0 = pega(saida_a, r'CLtot =\s+([-\d.]+)')
CM0 = pega(saida_a, r'Cmtot =\s+([-\d.]+)')

# Rodada B: ponto de projeto -> st e sb
saida_b = roda_avl('load aft.avl\noper\n'
                   f'm\nmn {MACH}\n\n'
                   f'd2 d2 0\nde\n1 {it_aft}\n\n'
                   f'a c {CL_PROJ}\nx\nst\n\nsb\n\n\nquit\n')
alpha_b = pega(saida_b, r'Alpha =\s+([-\d.]+)')
xnp = pega(saida_b, r'Xnp =\s+([-\d.]+)')
# O bloco st termina antes da linha do ponto neutro, para a razao de
# estabilidade espiral (que termina em "Cnb = ...") nao contaminar o token
st = tokens(saida_b.split('Stability-axis derivatives')[1]
            .split('Neutral point')[0])
sb = tokens(saida_b.split('Geometry-axis derivatives')[1])

# Rodada C: diferencas finitas de CD em alpha fixo. A ordem dos blocos de
# execucao segue a ordem dos comandos abaixo.
cmds = ('load aft.avl\noper\n'
        f'm\nmn {MACH}\n\n'
        f'd2 d2 0\nde\n1 {it_aft}\n\n'
        f'a a {alpha_b}\n'
        f'p p {DQ}\nx\n'
        f'p p {-DQ}\nx\n'
        'p p 0\n'
        f'de\n1 {it_aft + DIT}\n\nx\n'
        f'de\n1 {it_aft - DIT}\n\nx\n'
        f'de\n1 {it_aft}\n\n'
        f'd2 d2 {DDE}\nx\n'
        f'd2 d2 {-DDE}\nx\n'
        '\nquit\n')
saida_c = roda_avl(cmds)
cds = []
cls_fd = []
for bloco in saida_c.split('Vortex Lattice Output')[1:]:
    cds.append(pega_todos(bloco, r'CDtot =\s+([-\d.]+)')[0])
    cls_fd.append(pega_todos(bloco, r'CLtot =\s+([-\d.]+)')[0])
if len(cds) != 6:
    raise RuntimeError(f'esperava 6 execucoes na rodada C, veio {len(cds)}')
CDq = (cds[0] - cds[1])/(2*DQ)
CDit_deg = (cds[2] - cds[3])/(2*DIT)
CDde_deg = (cds[4] - cds[5])/(2*DDE)
CLq_fd = (cls_fd[0] - cls_fd[1])/(2*DQ)  # conferencia contra o st

# Ajuste parabolico da polar nao trimada do CG traseiro (item 5.b)
polar = np.genfromtxt('relatorio_lab04/tables/polar_aft_livre.csv',
                      delimiter=',', names=True)
alfa_rad = np.radians(polar['alpha_deg'])
c2, c1, c0 = np.polyfit(alfa_rad, polar['CD'], 2)

RAD = 180.0/np.pi

# (nome, valor no AVL, fator para o MVO, fonte)
LINHAS = [
    ('CL0', CL0, 1.0, 'ft, alpha = 0'),
    ('CLa', st['CLa'], 1.0, 'st'),
    ('CLq', st['CLq'], 1.0, 'st'),
    ('CLit', st['CLg1'], RAD, 'st (g1)'),
    ('CLde', st['CLd2'], RAD, 'st (d2)'),
    ('CD0', c0, 1.0, 'ajuste item 5.b'),
    ('CDa', c1, 1.0, 'ajuste item 5.b'),
    ('CDa2', c2, 1.0, 'ajuste item 5.b'),
    ('CDq', CDq, 1.0, 'dif. finita'),
    ('CDit', CDit_deg, RAD, 'dif. finita'),
    ('CDde', CDde_deg, RAD, 'dif. finita'),
    ('CM0', CM0, 1.0, 'ft, alpha = 0'),
    ('CMa', st['Cma'], 1.0, 'st'),
    ('CMq', st['Cmq'], 1.0, 'st'),
    ('CMit', st['Cmg1'], RAD, 'st (g1)'),
    ('CMde', st['Cmd2'], RAD, 'st (d2)'),
    ('CYb', st['CYb'], -1.0, 'st'),
    ('CYp', st['CYp'], -1.0, 'st'),
    ('CYr', st['CYr'], -1.0, 'st'),
    ('CYdr', st['CYd3'], RAD, 'st (d3)'),
    ('Clb', st['Clb'], 1.0, 'st'),
    ('Clp', sb['Clp'], 1.0, 'sb'),
    ('Clr', sb['Clr'], 1.0, 'sb'),
    ('Clda', sb['Cld1'], -RAD, 'sb (d1)'),
    ('Cldr', sb['Cld3'], -RAD, 'sb (d3)'),
    ('Cnb', st['Cnb'], 1.0, 'st'),
    ('Cnp', sb['Cnp'], 1.0, 'sb'),
    ('Cnr', sb['Cnr'], 1.0, 'sb'),
    ('Cnda', sb['Cnd1'], -RAD, 'sb (d1)'),
    ('Cndr', sb['Cnd3'], -RAD, 'sb (d3)'),
]

os.makedirs('relatorio_lab04/tables', exist_ok=True)
saidas = {}
with open('relatorio_lab04/tables/derivadas.csv', 'w',
          encoding='ascii') as f:
    f.write('parametro,valor_avl,valor_mvo,fonte\n')
    for nome, valor, fator, fonte in LINHAS:
        mvo = valor*fator
        saidas[nome] = {'avl': valor, 'mvo': mvo, 'fonte': fonte}
        f.write(f'{nome},{valor:.6g},{mvo:.6g},{fonte}\n')
        print(f'{nome:6s} AVL = {valor:12.6f}   MVO = {mvo:12.6f}   {fonte}')

extras = {'alpha_projeto_deg': alpha_b, 'it_aft_deg': it_aft,
          'xnp_st_m': xnp, 'CLq_fd_conferencia': CLq_fd}
saidas['_condicao'] = extras
with open('avl/saidas/derivadas.json', 'w', encoding='ascii') as f:
    json.dump(saidas, f, indent=2)

print(f'\nconferencia: CLq (st) = {st["CLq"]:.4f}  '
      f'CLq (dif. finita) = {CLq_fd:.4f}')
print(f'xnp (st) = {xnp:.4f} m; alpha de projeto = {alpha_b:.4f} graus')
print('gravado: relatorio_lab04/tables/derivadas.csv e '
      'avl/saidas/derivadas.json')
