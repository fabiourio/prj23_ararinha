'''
Otimizacao de uma secao da asa -- Lab 03, PRJ-23, equipe Ararinha.

    min   c_d(Al, Au, alpha)
    s.a.  c_l          = cl_ref        condicao de projeto da secao
          (t/c)_max   >= (t/c)_ref     espessura estrutural / tanque
          mint        <= 0.01          bordo de fuga fino
          mint        >= 0             superficies nao se cruzam
          t_01        >= 0.0368        substituta de cl_max (ver abaixo)
          batentes em Al, Au, alpha

Objetivo e restricoes de escoamento vem do eulerblock com metodo adjunto
(verificado contra diferencas finitas em verificacao_adjunto.py: erro < 0,4%
em dCD/dalpha e dCD/dAu0; a pior componente e dCL/dAu0, com 4,3%).

SOBRE A RESTRICAO DE SUSTENTACAO MAXIMA
---------------------------------------
A formulacao de Euler e nao-viscosa e nao enxerga cl_max. Sem uma restricao,
nada impede o otimizador de afiar o bordo de ataque para reduzir arrasto de
onda e destruir a sustentacao maxima -- o DOE mostrou variacao de 0,888 a
2,132 mexendo SO no nariz.

O cl_max do XFoil nao serve como restricao de gradiente: ele vem de varrer
alpha ate o solver viscoso divergir, e medimos saltos de 0,6 causados apenas
por falha de convergencia. Diferenca finita sobre isso da lixo.

Usamos entao um batente geometrico medido: t_01 (espessura em x/c = 1%) e o
descritor de maior correlacao com cl_max no DOE (rho = 0,756), e o limiar
t_01 >= 0,0368 separa cl_max >= 1,80 com 95,2% de precisao e 52,7% de
cobertura. E linear nos coeficientes CST, logo a derivada e exata e de graca.
O cl_max verdadeiro e verificado no XFoil ao final.

Rodar com:  python otimiza_secao.py <estacao>     (meio | raiz | ponta)
'''

import os
import pickle
import shutil
import sys
import time

import numpy as np
from scipy.optimize import Bounds, minimize
from scipy.special import comb

AQUI = os.path.dirname(os.path.abspath(__file__))
PACOTE = os.path.normpath(os.path.join(AQUI, '..', 'eulerblock', 'package'))
sys.path.insert(0, PACOTE)
sys.path.insert(0, AQUI)

from eulerblock import euler_mod as eb                          # noqa: E402

RES = os.path.join(AQUI, 'resultados')

# --- condicao de projeto (documento de projeto, secoes 3 e 4) ---------------
MACH_N = 0.7196          # M * cos(Lambda_c/2) -- normal ao enflechamento

ESTACOES = {
    # nome    eta     cl_ref   (t/c)_ref
    'raiz':  dict(eta=0.101, cl_ref=0.5997, tc_ref=0.2179),
    'meio':  dict(eta=0.398, cl_ref=0.7319, tc_ref=0.1772),
    'ponta': dict(eta=0.900, cl_ref=0.7677, tc_ref=0.1082),
}

# --- restricoes geometricas -------------------------------------------------
X_T01 = 0.01             # estacao do descritor substituto
MINT_MAX = 0.01          # bordo de fuga fino (formulacao do professor)
MINT_MIN = 0.0           # superficies nao se cruzam (acrescentado por nos)

# Substituta de cl_max: t_01 / sqrt(t_max) >= BLUNTEZ_MIN.
#
# O limiar em t_01 puro (>= 0,0368) tem 95,2% de precisao mas so 52,7% de
# cobertura, e a cobertura despenca para 25% nos perfis finos -- ele exige do
# perfil de ponta a mesma bluntez absoluta que do perfil de raiz, o que e
# conservador demais. Normalizar pela espessura INTEIRA (t_01/t_max) e pior
# ainda (rho cai de 0,756 para 0,619): o que governa o pico de succao e a
# bluntez absoluta, nao a relativa. A raiz quadrada e o meio-termo que o DOE
# escolheu: rho = 0,742, precisao 95,8% e cobertura 60,7%.
#
# Continua diferenciavel: t_01 e linear nos coeficientes CST e t_max vem da
# funcao KS do airfoil_mod, suave por construcao, com gradiente ja disponivel.
BLUNTEZ_MIN = 0.0979

# --- parametros do solver ---------------------------------------------------
NCHORD, NJ, S0 = 31, 49, 0.5e-2      # malha nivel 1,0
ITER, DT, CFL = 20000, 0.001, 0.2
RES_NK, RES_TOL = 1e-5, 1e-8
FTOL, MAXITER = 1e-6, 100

# NACA 1411 -- ponto de partida indicado no roteiro
AU_1411 = np.array([0.16146332, 0.18349204, 0.14126241, 0.18194397])
AL_1411 = np.array([-0.1489439, -0.10330027, -0.10305128, -0.10514982])

NVAR = 4
AL_LOWER, AL_UPPER = [-1.00] * NVAR, [-0.05] * NVAR
AU_LOWER, AU_UPPER = [0.05] * NVAR, [1.00] * NVAR
ALPHA_MIN, ALPHA_MAX = -3 * np.pi / 180, 9 * np.pi / 180


# ---------------------------------------------------------------------------
# t_01 e sua derivada -- analiticos, porque a superficie CST e LINEAR nos
# coeficientes para um x fixo:
#     y(x) = C(x) * sum_i A_i * comb(N,i) * x^i * (1-x)^(N-i)
# entao d t_01 / dA_i e apenas o termo de Bernstein, constante.
def _bernstein_cst(x, n_coef):
    N = n_coef - 1
    C = x ** 0.5 * (1.0 - x) ** 1.0
    return np.array([C * comb(N, i) * x ** i * (1.0 - x) ** (N - i)
                     for i in range(n_coef)])


B_T01 = _bernstein_cst(X_T01, NVAR)


def t01_de(Al, Au):
    return float(np.dot(Au - Al, B_T01))


def grad_t01():
    '''Gradiente de t_01 em relacao a xx = [Al, Au, alpha].'''
    return np.hstack([-B_T01, B_T01, [0.0]])


def bluntez(t01, t_max):
    '''Substituta de cl_max: t_01 / sqrt(t_max).'''
    return t01 / np.sqrt(t_max)


def grad_bluntez(t01, t_max, dt_max):
    '''
    d/dx [ t01 * t_max^(-1/2) ]
      = (dt01/dx) * t_max^(-1/2)  -  (1/2) * t01 * t_max^(-3/2) * (dt_max/dx)
    '''
    return (grad_t01() / np.sqrt(t_max)
            - 0.5 * t01 * t_max ** (-1.5) * np.asarray(dt_max))


# ---------------------------------------------------------------------------
def desmonta(xx):
    return np.asarray(xx[:NVAR]), np.asarray(xx[NVAR:2 * NVAR]), float(xx[-1])


def ponto_de_partida(tc_alvo):
    '''
    NACA 1411 com a espessura escalada para perto do alvo, mantendo o
    arqueamento. Em CST, espessura ~ (Au - Al) e arqueamento ~ (Au + Al)/2,
    entao escalar a diferenca preserva a linha media.

    Sem isso o ponto de partida do roteiro (t/c = 0,115) violaria a restricao
    de espessura em todas as tres estacoes.
    '''
    from descritores import descritores
    d0 = descritores(AU_1411, AL_1411)
    k = tc_alvo / d0['t_max']
    meio = (AU_1411 + AL_1411) / 2
    semi = (AU_1411 - AL_1411) / 2
    return meio - k * semi, meio + k * semi          # Al, Au


class Avaliador:
    '''
    Roda o eulerblock uma unica vez por ponto de projeto e serve o resultado
    ao objetivo, as restricoes e aos respectivos gradientes -- o SLSQP pede
    cada um separadamente, e sem esse cache o custo quadruplicaria.
    '''

    def __init__(self, cl_ref, tc_ref, pasta):
        self.cl_ref, self.tc_ref = cl_ref, tc_ref
        self.pasta = pasta
        self.hist = {'xx': [], 'CD': [], 'CL': [], 'maxt': [], 'mint': [],
                     't01': [], 'tempo': []}
        self.dados = []          # dicionario completo de cada ponto avaliado
        self.reinit = 0
        self.t0 = time.time()

    def __call__(self, xx, forcar=False):
        xx = np.asarray(xx, dtype=float)
        if not forcar:
            # devolve o dicionario DAQUELE ponto, nao o do ultimo avaliado --
            # o SLSQP reconsulta pontos antigos durante a busca em linha
            for i, xh in enumerate(self.hist['xx']):
                if np.linalg.norm(xh - xx) < 1e-11:
                    return self.dados[i]

        Al, Au, alpha = desmonta(xx)
        cwd = os.getcwd()
        try:
            os.chdir(self.pasta)
            r = eb.run_cst(Al, Au, NCHORD, alpha, MACH_N,
                           gamma=1.4, order=2,
                           iter=ITER, dt=DT, CFL=CFL, use_local_dt=1,
                           res_NK=RES_NK, res_tol=RES_TOL,
                           reinitialize=self.reinit, plot=False,
                           adj_funcs=['cl_jlow', 'cd_jlow'], NJ=NJ, s0=S0)
        finally:
            os.chdir(cwd)

        CL, CD = float(r['CL']), float(r['CD'])
        # se a solucao divergiu, nao reaproveite o campo na proxima chamada
        self.reinit = 0 if (np.isnan(CL) or np.isnan(CD)) else 1

        g = r['grads']
        dados = {
            'CL': CL, 'CD': CD,
            'maxt': float(np.real(r['maxt'])), 'mint': float(np.real(r['mint'])),
            't01': t01_de(Al, Au),
            'dCL': np.hstack([g['cl_jlow']['dAl'], g['cl_jlow']['dAu'],
                              [g['cl_jlow']['alpha']]]),
            'dCD': np.hstack([g['cd_jlow']['dAl'], g['cd_jlow']['dAu'],
                              [g['cd_jlow']['alpha']]]),
            'dmaxt': np.hstack([g['maxt']['dAl'], g['maxt']['dAu'], [0.0]]),
            'dmint': np.hstack([g['mint']['dAl'], g['mint']['dAu'], [0.0]]),
        }
        self._grava(xx, dados)
        return dados

    def _grava(self, xx, d):
        h = self.hist
        h['xx'].append(xx.copy())
        for k in ('CD', 'CL', 'maxt', 'mint', 't01'):
            h[k].append(d[k])
        h['tempo'].append(time.time() - self.t0)
        self.dados.append(d)
        n = len(h['CD'])
        print(f"[{n:3d}] {h['tempo'][-1]/60:6.1f} min  CD={d['CD']:.6f}  "
              f"CL={d['CL']:.5f} (alvo {self.cl_ref:.4f})  "
              f"t/c={d['maxt']:.4f} (min {self.tc_ref:.4f})  "
              f"t01={d['t01']:.5f}  alpha={xx[-1]*180/np.pi:.3f} deg",
              flush=True)
        with open(os.path.join(self.pasta, 'historico.pickle'), 'wb') as fid:
            pickle.dump({'hist': h, 'cl_ref': self.cl_ref,
                         'tc_ref': self.tc_ref, 'mach_n': MACH_N}, fid)


def otimiza(nome_estacao):
    est = ESTACOES[nome_estacao]
    cl_ref, tc_ref = est['cl_ref'], est['tc_ref']

    pasta = os.path.join(RES, f'otim_{nome_estacao}')
    shutil.rmtree(pasta, ignore_errors=True)
    os.makedirs(pasta)

    Al0, Au0 = ponto_de_partida(tc_ref)
    av = Avaliador(cl_ref, tc_ref, pasta)

    print(f'=== ESTACAO {nome_estacao.upper()} (eta = {est["eta"]}) ===')
    print(f'  M_n = {MACH_N}   cl_ref = {cl_ref}   (t/c)_ref = {tc_ref}')
    print(f'  t_01/sqrt(t/c) >= {BLUNTEZ_MIN} (substituta de cl_max >= 1,80)')
    print(f'  partida: NACA 1411 reescalado, t/c = {tc_ref:.4f}, '
          f't01 = {t01_de(Al0, Au0):.5f}\n', flush=True)

    # --- acha o alpha de partida com um passo de Newton usando o adjunto ---
    alpha = 3.0 * np.pi / 180
    for _ in range(2):
        d = av(np.hstack([Al0, Au0, alpha]))
        erro = d['CL'] - cl_ref
        if abs(erro) < 5e-3 or not np.isfinite(d['dCL'][-1]):
            break
        alpha -= erro / d['dCL'][-1]
        alpha = float(np.clip(alpha, ALPHA_MIN, ALPHA_MAX))
    print(f'  alpha de partida: {alpha*180/np.pi:.3f} deg\n', flush=True)

    xx0 = np.hstack([Al0, Au0, alpha])

    def objfun(xx):
        return av(xx)['CD']

    def objgrad(xx):
        return av(xx)['dCD']

    def eqfun(xx):
        return av(xx)['CL'] - cl_ref

    def eqgrad(xx):
        return av(xx)['dCL']

    def ineqfun(xx):
        d = av(xx)
        return np.array([
            d['maxt'] - tc_ref,                            # espessura exigida
            MINT_MAX - d['mint'],                          # bordo de fuga fino
            d['mint'] - MINT_MIN,                          # sem cruzar
            bluntez(d['t01'], d['maxt']) - BLUNTEZ_MIN,    # substituta clmax
        ])

    def ineqgrad(xx):
        d = av(xx)
        return np.vstack([d['dmaxt'], -d['dmint'], d['dmint'],
                          grad_bluntez(d['t01'], d['maxt'], d['dmaxt'])])

    cons = [{'type': 'ineq', 'fun': ineqfun, 'jac': ineqgrad},
            {'type': 'eq', 'fun': eqfun, 'jac': eqgrad}]
    bounds = Bounds(AL_LOWER + AU_LOWER + [ALPHA_MIN],
                    AL_UPPER + AU_UPPER + [ALPHA_MAX], keep_feasible=True)

    t0 = time.time()
    res = minimize(objfun, xx0, jac=objgrad, constraints=cons, bounds=bounds,
                   method='SLSQP',
                   options={'maxiter': MAXITER, 'ftol': FTOL, 'disp': True})
    dt = time.time() - t0

    d = av(res.x, forcar=True)
    Al, Au, alpha = desmonta(res.x)

    print(f'\n=== OTIMO -- {nome_estacao} ===')
    print(f'  {res.message}  ({res.nit} iteracoes, {len(av.hist["CD"])} '
          f'avaliacoes, {dt/60:.1f} min)')
    print(f'  CD  = {d["CD"]:.6f}   (partida {av.hist["CD"][0]:.6f}, '
          f'{(d["CD"]/av.hist["CD"][0]-1)*100:+.2f}%)')
    print(f'  CL  = {d["CL"]:.5f}   (alvo {cl_ref})')
    print(f'  t/c = {d["maxt"]:.5f}  (min {tc_ref})')
    bl = bluntez(d['t01'], d['maxt'])
    print(f'  t01 = {d["t01"]:.5f}   bluntez = {bl:.5f} '
          f'(min {BLUNTEZ_MIN})'
          f'{"  ATIVA" if abs(bl - BLUNTEZ_MIN) < 1e-4 else "  folgada"}')
    print(f'  Al  = {np.array2string(Al, precision=6)}')
    print(f'  Au  = {np.array2string(Au, precision=6)}')
    print(f'  alpha = {alpha*180/np.pi:.4f} deg')
    print(f'\n  historico em {os.path.join(pasta, "historico.pickle")}')
    return res


if __name__ == '__main__':
    estacao = sys.argv[1] if len(sys.argv) > 1 else 'meio'
    if estacao not in ESTACOES:
        print(f'estacao invalida: {estacao}. Use: {list(ESTACOES)}')
        sys.exit(1)
    os.makedirs(RES, exist_ok=True)
    otimiza(estacao)
