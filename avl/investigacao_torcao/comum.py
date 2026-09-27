'''
Investigacao da otimizacao de torcao (auditoria). Nada aqui altera os
arquivos entregues: os .avl temporarios vao para resultados/_tmp/inv_<tag>_*.avl
e as saidas para investigacao_torcao/.

Avaliador generalizado: mesmo objetivo e mesma restricao de
otimizacao_torcao.Avaliador (reproduz exatamente o modo 'completo' com
cg='aft', parametrizacao 'full', margem 0,2), com opcoes de parametrizacao,
margem, CG, modelo (sem empenagem) e registro de cada avaliacao/iteracao.
'''

import json
import os
import sys
import time

AQUI_INV = os.path.dirname(os.path.abspath(__file__))
AVL_DIR = os.path.dirname(AQUI_INV)
if AVL_DIR not in sys.path:
    sys.path.insert(0, AVL_DIR)

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize

from aeronave import aeronave
from analises import it_para_de_zero
from avl_run import caso
from estol import modelo_estol, estol, clmax_local
import gera_avl
from gera_avl import escreve_avl, ETAS
from otimizacao_torcao import LIMITES, ETA_ESTOL, DESEMPATE, EPS_FD
from sombra import alfa_saida

LC = 4.8237
ETAS_A = np.array(ETAS)
ETA_SOB = ETAS[1]          # lado da fuselagem (0,1011)
RUNS = os.path.join(AQUI_INV, 'runs')
os.makedirs(RUNS, exist_ok=True)


# ------------------------------------------------------------------ parametrizacoes
def param_para_torcao(nome, p):
    '''Devolve a torcao nas 11 estacoes (raiz = 0).'''
    p = np.asarray(p, float)
    if nome == 'full':
        return np.concatenate([[0.0], p])
    if nome == 'pl2':          # incidencia no lado da fuselagem + washout linear ate a ponta
        nos = [ETA_SOB, 1.0]
        return np.concatenate([[0.0], np.interp(ETAS_A[1:], nos, p)])
    if nome.startswith('pl3'):   # 'pl3' (quebra em 0,56) ou 'pl3_048'
        k = 0.48 if nome == 'pl3_048' else 0.56
        nos = [ETA_SOB, k, 1.0]
        return np.concatenate([[0.0], np.interp(ETAS_A[1:], nos, p)])
    if nome == 'pchip4':       # PCHIP pelos nos 0,1011 / 0,40 / 0,70 / 1,0
        nos = [ETA_SOB, 0.40, 0.70, 1.0]
        return np.concatenate([[0.0], PchipInterpolator(nos, p)(ETAS_A[1:])])
    raise ValueError(nome)


N_PARAM = {'full': 10, 'pl2': 2, 'pl3': 3, 'pl3_048': 3, 'pchip4': 4}


def ajusta_param(nome, torcao):
    '''Minimos quadrados: parametros cuja torcao mais se aproxima de `torcao`.'''
    from scipy.optimize import least_squares
    n = N_PARAM[nome]
    f = lambda p: param_para_torcao(nome, p)[1:] - np.asarray(torcao)[1:]
    return least_squares(f, np.zeros(n)).x


# ------------------------------------------------------------------ modelos
def escreve_modelo(av, caminho, torcao, cg='aft', modelo='completo', malha=None):
    if modelo == 'completo':
        return escreve_avl(av, caminho, cg=cg, torcao=torcao, malha=malha)
    # asa + winglet + fuselagem + nacele, sem EH/EV (sem compensacao)
    malha = gera_avl.malha_adotada() if malha is None else malha
    secoes = gera_avl.secoes_asa(av, torcao, ETAS, True)
    xref = av['xcg_fwd'] if cg == 'fwd' else av['xcg_aft']
    txt = (f'Ararinha Lc_h={av["Lc_h"]:.4f} CG {cg} {modelo}\n'
           f'{av["M"]:.4f}          # Mach\n0 0 0.0       # iYsym iZsym Zsym\n'
           f'{av["Sref"]:.4f} {av["Cref"]:.4f} {av["Bref"]:.4f}   # Sref Cref Bref\n'
           f'{xref:.4f} 0.0 0.0   # Xref Yref Zref\n{av["CD0"]:.5f}       # CDp\n')
    txt += gera_avl._asa(secoes, malha, controles=True)
    txt += gera_avl._winglet(secoes, malha)
    if modelo in ('asa_fus', 'asa_fus_sem_nac', 'asa_sem_nac'):
        pass
    if modelo in ('asa_fus', 'asa_fus_sem_nac'):
        txt += gera_avl._corpo(av, malha)
    if modelo in ('asa_fus',):
        txt += gera_avl._nacele(av, 2)
    destino = os.path.join(AVL_DIR, caminho)
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    with open(destino, 'w', encoding='utf-8') as f:
        f.write(txt)
    return caminho


# ------------------------------------------------------------------ avaliador
class Aval:
    def __init__(self, av, tag, param='full', margem=DESEMPATE, cg='aft',
                 modelo='completo', it=None, excl_frac=0.25, com_estol=True,
                 malha=None, alfa_min=None):
        self.av, self.tag, self.param, self.margem = av, tag, param, margem
        self.cg, self.modelo, self.malha = cg, modelo, malha
        self.emp = modelo == 'completo'
        self.excl_frac, self.com_estol, self.alfa_min = excl_frac, com_estol, alfa_min
        self.arquivo = f'resultados/_tmp/inv_{tag}_ot.avl'
        if it is None and self.emp:
            base = escreve_modelo(av, f'resultados/_tmp/inv_{tag}_base.avl',
                                  np.zeros(11), cg=cg, modelo=modelo, malha=malha)
            it, _ = it_para_de_zero(base, av['M'], av['CL'])
        self.it = 0.0 if it is None else it
        self.cache, self.evals, self.iters = {}, [], []
        self.t0 = time.time()

    def mascara(self, y):
        semi, d = self.av['asa']['yt'], self.excl_frac*self.av['asa']['ct']
        return ((semi - y) <= d) | (y <= d)

    def avalia(self, p):
        chave = tuple(np.round(np.asarray(p, float), 8))
        if chave in self.cache:
            return self.cache[chave]
        tw = param_para_torcao(self.param, p)
        escreve_modelo(self.av, self.arquivo, tw, cg=self.cg, modelo=self.modelo,
                       malha=self.malha)
        r = caso(self.arquivo, self.av['M'], cl=self.av['CL'], it=self.it,
                 trim=self.emp, empenagem=self.emp, faixas=True)
        res = {'CDff': 1e4*r['CDff'], 'e': r['e'], 'r': r}
        if self.com_estol:
            m = modelo_estol(self.arquivo, self.it, trim=self.emp, av=self.av,
                             empenagem=self.emp)
            y = m['eta']*self.av['asa']['yt']
            validas = ~self.mascara(y)
            eta_v, alfa_v = m['eta'][validas], m['alfa_i'][validas]
            dentro = eta_v <= ETA_ESTOL
            a_on = alfa_v[dentro].min()
            g = alfa_v[~dentro] - a_on - self.margem
            if self.alfa_min is not None:          # piso de alfa de estol (opcional)
                g = np.concatenate([g, [alfa_v.min() - self.alfa_min]])
            res.update(g=g, alfa_on=float(a_on),
                       eta_on=float(eta_v[dentro][np.argmin(alfa_v[dentro])]),
                       alfa_min_all=float(alfa_v.min()),
                       eta_min_all=float(eta_v[np.argmin(alfa_v)]),
                       margem_ext=float((alfa_v[~dentro] - a_on).min()),
                       modelo_estol=m)
        self.cache[chave] = res
        self.evals.append({'k': len(self.evals), 't': time.time() - self.t0,
                           'p': [float(v) for v in p], 'twist': tw.tolist(),
                           'CDff': res['CDff'], 'e': res['e'],
                           'gmin': float(np.min(res['g'])) if self.com_estol else None,
                           'alfa_on': res.get('alfa_on'), 'eta_on': res.get('eta_on'),
                           'margem_ext': res.get('margem_ext')})
        return res

    def objetivo(self, p):
        return self.avalia(p)['CDff']

    def restricoes(self, p):
        return self.avalia(p)['g']

    def salva(self, extra=None):
        d = {'tag': self.tag, 'param': self.param, 'margem': self.margem,
             'cg': self.cg, 'modelo': self.modelo, 'it': self.it,
             'evals': self.evals, 'iters': self.iters}
        if extra:
            d.update(extra)
        with open(os.path.join(RUNS, f'{self.tag}.json'), 'w') as f:
            json.dump(d, f, indent=1, default=float)


def otimiza(a, p0, eps=EPS_FD, maxiter=100, ftol=1e-3, limites=None):
    '''Mesmas opcoes do otimiza() de producao, com callback de iteracao.'''
    lim = limites or [LIMITES]*len(p0)
    cons = [{'type': 'ineq', 'fun': a.restricoes}] if a.com_estol else []

    def cb(xk):
        r = a.avalia(xk)
        a.iters.append({'it': len(a.iters) + 1, 'nev': len(a.evals),
                        'p': [float(v) for v in xk],
                        'twist': param_para_torcao(a.param, xk).tolist(),
                        'CDff': r['CDff'],
                        'gmin': float(np.min(r['g'])) if a.com_estol else None,
                        'alfa_on': r.get('alfa_on'), 'eta_on': r.get('eta_on')})
        a.salva()
        print(f'[{a.tag}] it {len(a.iters):3d} nev {len(a.evals):4d} '
              f'CDff {r["CDff"]:.3f} gmin {a.iters[-1]["gmin"]} '
              f't {time.time()-a.t0:.0f}s', flush=True)

    r = minimize(a.objetivo, np.asarray(p0, float), method='SLSQP', bounds=lim,
                 constraints=cons, callback=cb,
                 options={'ftol': ftol, 'eps': eps, 'maxiter': maxiter})
    return r


# ------------------------------------------------------------------ verificacao (como main() de producao)
def verifica(twist, tag, lc=LC, excl_frac=0.25):
    '''
    Mesma verificacao do main() de producao para uma torcao: aft/fwd com it
    retrimado, CDff/e no cruzeiro (aft), estol nos dois CG, folga de sombra,
    margens de estol sobre a regiao externa e sobre o aileron.
    '''
    av = aeronave(lc)
    tw = np.asarray(twist, float)
    out = {'twist': tw.tolist()}
    for cg in ('aft', 'fwd'):
        arq = escreve_avl(av, f'resultados/_tmp/inv_{tag}_v{cg}.avl', cg=cg, torcao=tw)
        it, r = it_para_de_zero(arq, av['M'], av['CL'])
        e = estol(arq, it=it, trim=True, av=av)
        m = modelo_estol(arq, it, trim=True, av=av)
        y = m['eta']*av['asa']['yt']
        semi, d = av['asa']['yt'], excl_frac*av['asa']['ct']
        val = ~(((semi - y) <= d) | (y <= d))
        eta, al = m['eta'], m['alfa_i']
        ail = val & (eta >= 0.56) & (eta <= 0.90)
        ext = val & (eta > ETA_ESTOL)
        out[cg] = {'it': it, 'CDff': 1e4*r['CDff'], 'e': r['e'],
                   'alfa': e['alfa'], 'CL': e['CL'], 'de': e['de'],
                   'eta_crit': e['eta_crit'],
                   'margem_aileron': float(al[ail].min() - e['alfa']),
                   'margem_externa': float(al[ext].min() - e['alfa']),
                   'eta': e['eta'].tolist(), 'cl_norm': e['cl_norm'].tolist(),
                   'clmax': e['clmax'].tolist(), 'alfa_i': m['alfa_i'].tolist(),
                   'excluidas_mask': (~val).tolist()}
    a_s = alfa_saida(av)
    a_min = min(out['aft']['alfa'], out['fwd']['alfa'])
    out['sombra'] = {'alfa_estol_min': a_min, 'alfa_saida': a_s, 'folga': a_min - a_s}
    return out
