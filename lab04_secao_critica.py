'''
Lab 04, item 4 -- CLmax de asa limpa pelo metodo da secao critica.

Para cada caso (CG dianteiro/traseiro, com e sem trimagem de arfagem), em
M = 0,2 e com a incidencia de empenagem do item 3, varre o angulo de
ataque ate a distribuicao de cl_norm da asa alcancar o limite de clmax dos
perfis do Lab 03 (plano normal, valores corrigidos pela calibracao:
revisao_plano_clmax.csv). Como o VLM e linear, o cruzamento e resolvido
por interpolacao entre os pontos da varredura.

Gera:
  - avl/saidas/secao_critica.json (alpha_max, CLmax e delta_e por caso)
  - relatorio_lab04/tables/clxy_<caso>.csv (distribuicao cl_norm x eta na
    condicao de estol, com a curva limite, para as figuras)

Rodar da raiz do repo:  python lab04_secao_critica.py
(requer avl/saidas/trim.json gerado por lab04_trim.py)
'''

# IMPORTS
import json
import os
import re

import numpy as np

from avl_batch import roda_avl, pega_todos

#=========================================

MACH = 0.2
ALPHAS = np.arange(4.0, 24.1, 2.0)
SEMI_ENV = 30.0759                       # b/2 [m]

# Limite de clmax por estacao (plano normal, corrigido -- Lab 03)
ETA_LIM = np.array([0.1011, 0.398, 0.90])
CLMAX_LIM = np.array([1.774, 1.7985, 1.7338])

CASOS = {'fwd_livre': ('fwd', 'd2 d2 0'),
         'fwd_trim': ('fwd', 'd2 pm 0'),
         'aft_livre': ('aft', 'd2 d2 0'),
         'aft_trim': ('aft', 'd2 pm 0')}

RE_LINHA = re.compile(r'^\s*\d+\s+(' + r'[-\dEe.+]+\s+'*11 +
                      r'[-\dEe.+]+)\s*$', re.M)


def limite(eta):
    return np.interp(eta, ETA_LIM, CLMAX_LIM)


def blocos_asa(saida):
    '''Distribuicoes (eta, cl_norm) da asa direita, uma por comando fs.'''
    blocos = []
    for parte in saida.split('Surface # 1     Wing')[1:]:
        parte = parte.split('Surface # 2')[0]
        linhas = [[float(v) for v in m.group(1).split()]
                  for m in RE_LINHA.finditer(parte)]
        dados = np.array(linhas)
        eta = dados[:, 0]/SEMI_ENV       # Yle -> eta
        blocos.append((eta, dados[:, 5]))  # coluna cl_norm
    return blocos


def varre(cg, restricao, it, alphas):
    '''Roda a sequencia de alphas num unico processo do AVL e devolve as
    listas de alpha, CL, delta_e e as distribuicoes da asa. Os valores sao
    extraidos bloco a bloco de execucao, para nao confundir com as linhas
    de menu do AVL.'''
    cmds = (f'load {cg}.avl\noper\n'
            f'm\nmn {MACH}\n\n'
            f'de\n1 {it}\n\n'
            f'{restricao}\n')
    for alpha in alphas:
        cmds += f'a a {alpha}\nx\nfs\n\n'
    cmds += '\nquit\n'
    saida = roda_avl(cmds, timeout=580)
    alphas_out, cls, des, dists = [], [], [], []
    for bloco in saida.split('Vortex Lattice Output')[1:]:
        alphas_out.append(pega_todos(bloco, r'Alpha =\s+([-\d.]+)')[0])
        cls.append(pega_todos(bloco, r'CLtot =\s+([-\d.]+)')[0])
        des.append(pega_todos(bloco, r'elevator\s+=\s+([-\d.]+)')[0])
        dists.extend(blocos_asa(bloco))
    return alphas_out, cls, des, dists


with open('avl/saidas/trim.json', encoding='ascii') as f:
    trim = json.load(f)

os.makedirs('relatorio_lab04/tables', exist_ok=True)
resultado = {}
for caso, (cg, restricao) in CASOS.items():
    it = trim[cg]['it_deg']
    alphas, cls, des, dists = varre(cg, restricao, it, ALPHAS)
    margens = np.array([np.max(cln - limite(eta)) for eta, cln in dists])

    if margens[0] > 0 or margens[-1] < 0:
        raise RuntimeError(f'{caso}: estol fora da faixa varrida')
    k = int(np.argmax(margens > 0))      # primeiro alpha com estouro
    frac = -margens[k-1]/(margens[k] - margens[k-1])
    alpha_max = ALPHAS[k-1] + frac*(ALPHAS[k] - ALPHAS[k-1])

    # Rodada de verificacao exatamente no alpha de estol
    a2, cl2, de2, dist2 = varre(cg, restricao, it, [round(alpha_max, 3)])
    eta, cln = dist2[0]
    margem_res = float(np.max(cln - limite(eta)))
    eta_crit = float(eta[np.argmax(cln - limite(eta))])

    with open(f'relatorio_lab04/tables/clxy_{caso}.csv', 'w',
              encoding='ascii') as f:
        f.write('eta,cl_norm,clmax_lim\n')
        for e, c in zip(eta, cln):
            f.write(f'{e:.4f},{c:.4f},{limite(e):.4f}\n')

    resultado[caso] = {'alpha_max_deg': round(float(alpha_max), 3),
                       'CLmax': round(cl2[0], 4),
                       'delta_e_deg': round(de2[0], 3),
                       'eta_critica': round(eta_crit, 3),
                       'margem_residual': round(margem_res, 5)}
    print(f"{caso}: alpha_max = {alpha_max:.2f} graus, "
          f"CLmax = {cl2[0]:.4f}, delta_e = {de2[0]:.3f} graus, "
          f"estacao critica eta = {eta_crit:.3f}")

with open('avl/saidas/secao_critica.json', 'w', encoding='ascii') as f:
    json.dump(resultado, f, indent=2)
print('gravado: avl/saidas/secao_critica.json')
