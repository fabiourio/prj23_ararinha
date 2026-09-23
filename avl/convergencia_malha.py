'''
Etapa 0a -- convergencia de malha.

Geometria de base (Lc_h = 4,6, sem torcao, CG traseiro), ponto de projeto
compensado com it que zera o profundor. Refina uma direcao por vez, na
ordem de SEQ, ja usando o que foi adotado nas anteriores. Adota a menor malha
cuja diferenca para a mais fina da sequencia fica abaixo de TOL.

A regra de adocao percorre os valores da sequencia da mais fina para a mais
grossa e adota o menor N tal que ele e TODOS os mais finos que ele ficam
dentro da tolerancia da malha mais fina. Isso evita adotar uma malha grossa
por coincidencia, quando alguma grandeza monitorada oscila (nao-monotonica)
e um ponto intermediario furou a tolerancia.

Tres correcoes de modelagem, feitas antes de rodar este estudo (ver
docstrings de gera_avl.py e estol.py para os detalhes):
  1. estol.py exclui da secao critica a(s) faixa(s) da asa na juncao com o
     winglet (singularidade de vórtice de canto): sem isso alfa_estol caia
     de ~10 para ~5 graus so por refinar asa_ns.
  2. gera_avl.py reamostra fuselage_nondim.dat (fuselage_reamostrada.dat):
     o arquivo original faz a spline de comprimento de arco do AVL laçar
     perto da cauda, e CDff/it so convergem com a malha do corpo a partir
     de nbody ~= 80 com o arquivo reamostrado (o original nao converge).
  3. gera_avl.py poe o anel da nacele em COMPONENT proprio (nacele_componente
     default 2): no mesmo componente da asa, CDff oscila +-0,2 count com a
     malha da envergadura da asa.
Com as tres, CDff, it e alfa_estol convergem nas tres sequencias.

Rodar de dentro de avl/:
    python convergencia_malha.py            # estudo completo
    python convergencia_malha.py --lc 4.8   # confere a malha adotada em outro Lc_h
'''

import argparse
import json
import os

from estilo import plt, AZUL, LARANJA, TINTA2
from aeronave import aeronave, LC_H_BASE
from analises import it_para_de_zero
from estol import estol
from gera_avl import escreve_avl, MALHA_PADRAO, ARQ_MALHA, malha_adotada, AQUI

# asa_ns comeca em 20: com as 11 secoes definidas em ETAS (gera_avl.py), o
# AVL 3.37 recusa malhas mais grossas ("Cannot adjust spanwise spacing... /
# Insufficient number of spanwise vortices to work with"), pois nao consegue
# alocar vortices suficientes para respeitar o espacamento pedido em cada
# secao. Confirmado tentando 10, 12 e 14 antes de fixar o piso em 20.
#
# nbody vai ate 100: em 120 o AVL 3.37 estoura um limite de array fixo no
# binario ("MAKEBODY: Array overflow. Increase KLMAX to 120" / "...to 102"),
# que so se resolve recompilando o AVL com KLMAX maior (fora do escopo
# deste lab). O teto foi confirmado tentando 120, 150 e 200 (todos falham
# com a mesma mensagem) antes de fixar o teto em 100. 10 e 20 foram tirados
# da sequencia: com a fuselagem reamostrada CDff/it so convergem a partir
# de nbody ~= 80, entao a malha inicial deste estudo ja usa nbody=80 (ver
# `estudo`) para as sequencias de asa/EH nao serem poluidas por um corpo
# nao convergido.
SEQ = [('asa_ns', [20, 30, 40, 60, 80]), ('asa_nc', [4, 8, 12, 16]),
       ('winglet_ns', [4, 8, 12]), ('eh_ns', [5, 10, 15, 20, 30]),
       ('eh_nc', [8, 16, 24]), ('nbody', [40, 60, 80, 100])]
TOL = {'CDff': 0.1, 'it': 0.02, 'alfa_estol': 0.15}   # [count], [graus], [graus]
ARQ = 'resultados/_tmp/malha.avl'


def metricas(av, malha):
    escreve_avl(av, ARQ, cg='aft', malha=malha)
    it, r = it_para_de_zero(ARQ, av['M'], av['CL'])
    e = estol(ARQ, it=it, trim=True, av=av)
    return {'CDff': 1e4*r['CDff'], 'it': it, 'alfa_estol': e['alfa']}


def escolhe_adotado(valores, res):
    '''
    Menor N tal que ele e todos os mais finos que ele ficam dentro da
    tolerancia da malha mais fina (percorre de tras para frente, ao inves de
    aceitar o primeiro da frente que passar, para nao adotar uma malha
    grossa por coincidencia quando a grandeza monitorada oscila).
    '''
    fino = res[-1]
    adotado = valores[-1]
    for n, r in zip(valores, res):
        dentro = all(abs(r[k] - fino[k]) < TOL[k] for k in TOL)
        if dentro:
            adotado = n
            break
    return adotado


def estudo():
    av = aeronave(LC_H_BASE)
    # nbody=80 desde o inicio (ver comentario de SEQ): as sequencias de asa
    # e EH nao ficam poluidas por um corpo nao convergido.
    malha = dict(MALHA_PADRAO, nbody=80)
    historico = {}
    for chave, valores in SEQ:
        res = []
        for n in valores:
            m = dict(malha, **{chave: n})
            res.append(metricas(av, m))
            print(f'  {chave:10s} = {n:3d}   CDff = {res[-1]["CDff"]:8.3f} count'
                  f'   it = {res[-1]["it"]:7.3f}   alfa_estol = {res[-1]["alfa_estol"]:6.2f}')
        adotado = escolhe_adotado(valores, res)
        malha[chave] = adotado
        historico[chave] = {'valores': valores, 'res': res, 'adotado': adotado}
        print(f'  -> {chave} adotado = {adotado}\n')

    os.makedirs(os.path.dirname(ARQ_MALHA), exist_ok=True)
    with open(ARQ_MALHA, 'w') as f:
        json.dump(malha, f, indent=2)
    with open(os.path.join(AQUI, 'resultados', 'convergencia_malha.json'), 'w') as f:
        json.dump(historico, f, indent=2)
    figura(historico)
    print('malha adotada:', malha)


def figura(hist):
    fig, eixos = plt.subplots(3, len(hist), figsize=(2.3*len(hist), 6.0), sharex='col')
    for j, (chave, h) in enumerate(hist.items()):
        v = h['valores']
        for i, (k, cor, rot, tol) in enumerate((
                ('CDff', AZUL, 'CDff [count]', TOL['CDff']),
                ('it', LARANJA, 'it [°]', TOL['it']),
                ('alfa_estol', TINTA2, 'alfa_estol [°]', TOL['alfa_estol']))):
            y = [r[k] for r in h['res']]
            ax = eixos[i, j]
            ax.axhspan(y[-1] - tol, y[-1] + tol, color='#e1e0d9', lw=0)
            ax.plot(v, y, '-o', color=cor)
            ax.axvline(h['adotado'], color=TINTA2, lw=0.8, ls='--')
            if j == 0:
                ax.set_ylabel(rot)
        eixos[2, j].set_xlabel(chave)
    fig.suptitle('Convergência de malha (faixa cinza: tolerância em torno da malha mais fina)')
    fig.savefig(os.path.join(AQUI, 'resultados', 'convergencia_malha.png'))


def confere(lc):
    av = aeronave(lc)
    m = malha_adotada()
    base = metricas(av, m)
    fina = metricas(av, dict(m, asa_ns=2*m['asa_ns'], eh_ns=2*m['eh_ns']))
    d = {k: fina[k] - base[k] for k in base}
    ok = all(abs(d[k]) < TOL[k] for k in TOL)
    print(f'Lc_h = {lc}: dCDff = {d["CDff"]:+.3f} count, dit = {d["it"]:+.4f} graus,'
          f' dalfa_estol = {d["alfa_estol"]:+.4f} graus -> {"ok" if ok else "REVER MALHA"}')
    return ok


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--lc', type=float, default=None)
    a = p.parse_args()
    confere(a.lc) if a.lc is not None else estudo()
