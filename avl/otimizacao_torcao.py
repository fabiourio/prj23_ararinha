'''
Otimizacao de torcao da asa.

Variaveis: Ainc das estacoes da asa, com a raiz fixa em 0 (a raiz define o
eixo do corpo). Objetivo: CDff (Trefftz) no CL de projeto, M = 0,85 -- o
CDp e constante, entao minimizar CDff e minimizar o CD.

modelo='completo': aft.avl com compensacao pelo profundor (d2 pm 0) e it
fixo; restricao de estol: toda faixa externa a ETA_ESTOL tem de estolar
pelo menos DESEMPATE graus depois da primeira faixa interna.
modelo='asa_limpa': so a asa, plana, sem winglet, sem compensacao
(verificacao eliptica).

O mesmo codigo serve as duas coisas, de proposito: a verificacao eliptica
testa exatamente o otimizador usado no projeto.

Restricao de estol (res['g']): modelo_estol exclui as faixas da juncao com
o winglet e da raiz (mascara 'excluidas', ver docstring de estol.py) --
essas faixas nao entram na restricao. As faixas remanescentes sao fixas
para uma dada malha/etas (a mascara depende so da geometria, nao de alfa),
entao o comprimento do vetor de restricao e estavel entre avaliacoes, como
o SLSQP exige.
'''

import numpy as np
from scipy.optimize import minimize

from avl_run import caso
from estol import modelo_estol
from gera_avl import escreve_avl, ETAS

LIMITES = (-8.0, 3.0)
ETA_ESTOL = 0.56        # raiz do aileron: o estol tem de comecar para dentro
DESEMPATE = 0.2         # [graus]
EPS_FD = 0.2            # [graus] passo das diferencas finitas


def elipse_normalizada(eta):
    return np.sqrt(np.clip(1.0 - np.asarray(eta)**2, 0.0, None))


def carga_normalizada(av, r):
    '''c·cl dividido pela carga eliptica na raiz, 4 Sref CL / (pi b).'''
    f = r['faixas']
    return f['y']/av['asa']['yt'], f['ccl']/(4*av['Sref']*r['CL']/(np.pi*av['Bref']))


class Avaliador:
    def __init__(self, av, modelo='completo', etas=ETAS, it=0.0,
                 com_estol=True, nome='ot'):
        self.av, self.modelo, self.etas = av, modelo, list(etas)
        self.it, self.com_estol = it, com_estol
        self.arquivo = f'resultados/_tmp/{nome}.avl'
        self.cache, self.n_rodadas, self.ultimo = {}, 0, None
        self.historico = []

    def torcao(self, x):
        return np.concatenate([[0.0], np.asarray(x, float)])

    def avalia(self, x):
        chave = tuple(np.round(np.asarray(x, float), 8))
        if chave in self.cache:
            self.ultimo = self.cache[chave]
            return self.ultimo
        limpa = self.modelo == 'asa_limpa'
        escreve_avl(self.av, self.arquivo, cg='aft', torcao=self.torcao(x),
                    etas=self.etas, so_asa=limpa, diedro=not limpa)
        r = caso(self.arquivo, self.av['M'], cl=self.av['CL'], it=self.it,
                 trim=not limpa, empenagem=not limpa, faixas=True)
        self.n_rodadas += 1
        res = {'CDff': r['CDff'], 'e': r['e'], 'r': r}
        if self.com_estol:
            m = modelo_estol(self.arquivo, self.it, trim=not limpa, av=self.av,
                             empenagem=not limpa)
            self.n_rodadas += 2
            validas = ~m['excluidas']
            eta_v, alfa_v = m['eta'][validas], m['alfa_i'][validas]
            dentro = eta_v <= ETA_ESTOL
            res['g'] = alfa_v[~dentro] - alfa_v[dentro].min() - DESEMPATE
            res['estol'] = m
        self.cache[chave] = res
        self.ultimo = res
        self.historico.append((np.array(x, float), 1e4*res['CDff']))
        return res

    def objetivo(self, x):
        return 1e4*self.avalia(x)['CDff']          # [count]

    def restricoes(self, x):
        return self.avalia(x)['g']


def otimiza(avaliador, x0, maxiter=100):
    cons = ([{'type': 'ineq', 'fun': avaliador.restricoes}]
            if avaliador.com_estol else [])
    r = minimize(avaliador.objetivo, np.asarray(x0, float), method='SLSQP',
                 bounds=[LIMITES]*len(x0), constraints=cons,
                 options={'ftol': 1e-3, 'eps': EPS_FD, 'maxiter': maxiter})
    return r
