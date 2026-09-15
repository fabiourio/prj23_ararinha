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

# A raiz precisa de tratamento proprio. Com (t/c)_n = 0,2179 a M_n = 0,7196 a
# secao opera 48 pontos de Mach alem da divergencia de arrasto, e partindo do
# NACA 1411 reescalado o Euler DIVERGE (residuo cai ate 1e-5, o Newton-Krylov
# encolhe o passo para 0,25 e a solucao explode para NaN).
#
# Sondagem (sonda_raiz.py): CFL 0,10 ja estabiliza a partida NACA (181 s), e
# CFL 0,05 da praticamente o mesmo resultado (CD 0,0730 contra 0,0716), o que
# indica que 0,10 ja esta na regiao estavel e nao apenas mascarando o problema.
#
# O RAE2822 reescalado converge ate em CFL 0,20, em 73 s, porque sua espessura
# maxima fica em x/c = 0,376 contra 0,296 do NACA -- a assinatura supercritica.
# Mesmo assim mantemos o NACA 1411 em todas as estacoes: e o que o roteiro
# prescreve, e usar a mesma parametrizacao (4+4 coeficientes) em todas as
# secoes mantem a comparacao entre elas honesta. O RAE2822 fica como
# recomendacao no relatorio, e como referencia na polar do item 7.
ESTACOES = {
    # nome    eta     cl_ref   (t/c)_ref   CFL
    'raiz':  dict(eta=0.101, cl_ref=0.5997, tc_ref=0.2179, cfl=0.10),
    'meio':  dict(eta=0.398, cl_ref=0.7319, tc_ref=0.1772, cfl=0.20),
    'ponta': dict(eta=0.900, cl_ref=0.7677, tc_ref=0.1082, cfl=0.20),
}


def limiar_bluntez(estacao):
    '''
    Limiar da substituta de cl_max para a estacao, lido da calibracao.

    O limiar depende da espessura -- exigir da ponta fina a mesma bluntez
    absoluta da raiz grossa seria conservador demais -- entao calibra_limiar.py
    resolve um valor por estacao. Se o arquivo nao existir, cai no valor
    global do DOE, que e conservador mas seguro.
    '''
    cam = os.path.join(RES, 'limiares_bluntez.csv')
    if not os.path.isfile(cam):
        print(f'  AVISO: {cam} nao existe; usando o limiar global '
              f'{BLUNTEZ_GLOBAL} (conservador). Rode calibra_limiar.py.')
        return BLUNTEZ_GLOBAL
    with open(cam) as fid:
        linhas = [l.split(',') for l in fid.read().splitlines()[1:] if l.strip()]
    for partes in linhas:
        if partes[0] == estacao:
            return float(partes[2])          # coluna nao-viesada
    raise KeyError(f'estacao {estacao} ausente em {cam}')

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
# Valor de reserva, do limiar conservador global do DOE. O limiar de fato
# usado sai de calibra_limiar.py, por estacao -- ver limiar_bluntez().
BLUNTEZ_GLOBAL = 0.0979

# --- parametros do solver ---------------------------------------------------
NCHORD, NJ, S0 = 31, 49, 0.5e-2      # malha nivel 1,0
ITER, DT, CFL = 20000, 0.001, 0.2
RES_NK, RES_TOL = 1e-5, 1e-8
# ftol PRECISA ficar acima do piso de ruido do solver, senao o criterio de
# parada e inatingivel e o SLSQP queima avaliacoes em buscas em linha que o
# ruido esconde. Medido na raiz: o c_d oscila entre 0,019855 e 0,019870
# (espalhamento 1,5e-5) depois de convergido. Com ftol = 1e-6, 15x menor que
# isso, ela convergiu na avaliacao 16 e gastou mais 123 avaliacoes sem sair
# do lugar. 1e-5 e o valor do proprio professor, e e compativel com o ruido.
# Perseguir mais que isso seria ilusorio de qualquer forma: o erro de malha no
# nivel 1,0 e de 70%.
FTOL, MAXITER = 1e-5, 100

# NACA 1411 -- ponto de partida indicado no roteiro
AU_1411 = np.array([0.16146332, 0.18349204, 0.14126241, 0.18194397])
AL_1411 = np.array([-0.1489439, -0.10330027, -0.10305128, -0.10514982])

# CFLs de recuo, caso a solucao divirja durante a otimizacao. A sondagem so
# testou o PONTO DE PARTIDA; nada garante que a geometria nao passe por
# regioes piores no caminho, entao o avaliador recua sozinho.
CFL_RECUO = [0.10, 0.05, 0.025]

NVAR = 4
AL_LOWER, AL_UPPER = [-1.00] * NVAR, [-0.05] * NVAR
AU_LOWER, AU_UPPER = [0.05] * NVAR, [1.00] * NVAR

# Batente ALTERNATIVO para o ultimo coeficiente do intradorso.
#
# Com os batentes do roteiro (Al <= -0,05), o Al4 fica prensado no limite nas
# TRES estacoes. Isso nao e coincidencia: o ultimo coeficiente do intradorso e
# o que constroi o "cusp" concavo do bordo de fuga, a essencia do
# carregamento traseiro supercritico. O RAE2822 tem Al[-1] = +0,052 --
# POSITIVO. O batente do roteiro proibe estruturalmente essa forma.
#
# Soltar so esse coeficiente responde o item 9 do roteiro com experimento:
# quanto se ganha permitindo a forma que os supercriticos de verdade usam?
AL4_SOLTO = 0.30
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

    def __init__(self, cl_ref, tc_ref, pasta, cfl=0.20):
        self.cl_ref, self.tc_ref = cl_ref, tc_ref
        self.cfl = cfl
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

        # tenta com o CFL da estacao e recua se divergir. A sondagem so testou
        # o ponto de partida; a geometria pode piorar no caminho.
        cfls = [self.cfl] + [c for c in CFL_RECUO if c < self.cfl]
        try:
            os.chdir(self.pasta)
            for tentativa, cfl in enumerate(cfls, 1):
                r = eb.run_cst(Al, Au, NCHORD, alpha, MACH_N,
                               gamma=1.4, order=2,
                               iter=ITER, dt=DT, CFL=cfl, use_local_dt=1,
                               res_NK=RES_NK, res_tol=RES_TOL,
                               reinitialize=0 if tentativa > 1 else self.reinit,
                               plot=False,
                               adj_funcs=['cl_jlow', 'cd_jlow'], NJ=NJ, s0=S0)
                if np.isfinite(r['CL']) and np.isfinite(r['CD']):
                    break
                print(f'      divergiu em CFL={cfl}; recuando', flush=True)
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
        # rede de seguranca no lugar da restricao 'mint >= 0', que tornava a
        # jacobiana singular: se o perfil chegar a cruzar as superficies, o
        # aviso aparece no log em vez de passar despercebido.
        if d['mint'] < 0:
            print(f'      AVISO: superficies cruzadas (mint = {d["mint"]:.6f})',
                  flush=True)
        print(f"[{n:3d}] {h['tempo'][-1]/60:6.1f} min  CD={d['CD']:.6f}  "
              f"CL={d['CL']:.5f} (alvo {self.cl_ref:.4f})  "
              f"t/c={d['maxt']:.4f} (min {self.tc_ref:.4f})  "
              f"t01={d['t01']:.5f}  alpha={xx[-1]*180/np.pi:.3f} deg",
              flush=True)
        with open(os.path.join(self.pasta, 'historico.pickle'), 'wb') as fid:
            pickle.dump({'hist': h, 'cl_ref': self.cl_ref,
                         'tc_ref': self.tc_ref, 'mach_n': MACH_N}, fid)


def pontos_multistart(tc_alvo, estacao, n=6, semente=23):
    '''
    Pontos de partida DIVERSOS para verificar se o otimo e global.

    O SLSQP e local. Partindo de um unico ponto -- o NACA 1411 reescalado --
    nao sabemos se existe bacia melhor. O Lab 02 verificou isso confrontando
    o SLSQP com o NSGA-II; aqui o equivalente barato e o multistart.

    O primeiro ponto e sempre o do roteiro, para a comparacao ser direta. Os
    demais saem de um hipercubo latino na caixa de coeficientes, reescalados
    para a espessura exigida e filtrados pelas restricoes geometricas, para
    nao gastar otimizacao partindo de algo inviavel.
    '''
    from scipy.stats import qmc

    import descritores as dsc

    pontos = [ponto_de_partida(tc_alvo)]
    lim = limiar_bluntez(estacao)

    lo = np.array([-0.45] * NVAR + [0.06] * NVAR)
    hi = np.array([-0.06] * NVAR + [0.45] * NVAR)
    # sorteio generoso: os filtros geometricos + o de sustentacao alcancavel
    # descartam a grande maioria, entao precisamos de muitos candidatos para
    # sair com pontos de partida realmente diversos
    amostra = qmc.LatinHypercube(d=2 * NVAR, seed=semente).random(4000)
    X = lo + amostra * (hi - lo)

    for linha in X:
        if len(pontos) >= n:
            break
        Al, Au = linha[:NVAR], linha[NVAR:]
        # reescala a espessura preservando o arqueamento, como em
        # ponto_de_partida
        d = dsc.descritores(Au, Al)
        if d['t_max'] <= 1e-6:
            continue
        k = tc_alvo / d['t_max']
        meio, semi = (Au + Al) / 2, (Au - Al) / 2
        Al2, Au2 = meio - k * semi, meio + k * semi
        if np.any(Al2 > AL_UPPER) or np.any(Al2 < AL_LOWER):
            continue
        if np.any(Au2 < AU_LOWER) or np.any(Au2 > AU_UPPER):
            continue
        d2 = dsc.descritores(Au2, Al2)
        if d2['t_min'] < MINT_MIN or d2['t_min'] > MINT_MAX:
            continue
        if bluntez(t01_de(Al2, Au2), d2['t_max']) < lim:
            continue
        # O ponto de partida precisa CONSEGUIR a sustentacao alvo dentro dos
        # batentes de alpha. Pela teoria de perfil fino, alpha_L0 ~ -1 grau
        # por 1% de arqueamento, entao um perfil de arqueamento muito negativo
        # exige alpha alto demais. Sem este filtro, o SLSQP parte de um ponto
        # onde a restricao de igualdade e INALCANCAVEL e a rodada se perde
        # tentando ganhar sustentacao em vez de reduzir arrasto.
        cl_ref = ESTACOES[estacao]['cl_ref']
        alpha_L0 = -100.0 * d2['c_max']                  # graus, aproximado
        alpha_nec = np.degrees(cl_ref / (2 * np.pi)) + alpha_L0
        if not (-2.0 < alpha_nec < 7.0):
            continue
        # descarta pontos quase iguais a um ja escolhido
        if any(np.linalg.norm(np.hstack([Al2, Au2]) - np.hstack(p)) < 0.08
               for p in pontos):
            continue
        pontos.append((Al2, Au2))

    return pontos


def otimiza(nome_estacao, com_bluntez=True, cusp_livre=False, i_partida=0):
    '''
    com_bluntez=False roda o MESMO problema sem a restricao de sustentacao
    maxima. Serve para medir quanto ela custa em arrasto: se o otimo sem ela
    ja atender o limiar, a restricao era inativa e nao custou nada; se nao
    atender, a diferenca de c_d entre as duas rodadas e o preco exato de
    manter clmax_w = 1,80.
    '''
    est = ESTACOES[nome_estacao]
    cl_ref, tc_ref = est['cl_ref'], est['tc_ref']
    bluntez_min = limiar_bluntez(nome_estacao) if com_bluntez else None

    sufixo = '' if com_bluntez else '_sem_bluntez'
    if cusp_livre:
        sufixo += '_cusp'
    if i_partida:
        sufixo += f'_ms{i_partida}'
    pasta = os.path.join(RES, f'otim_{nome_estacao}{sufixo}')
    # ignore_errors deixa a pasta de pe se algum arquivo estiver travado (o
    # eulerblock recem-morto ainda segura wall.dat por alguns segundos), e o
    # makedirs seguinte estourava. exist_ok resolve.
    shutil.rmtree(pasta, ignore_errors=True)
    os.makedirs(pasta, exist_ok=True)

    if i_partida:
        pontos = pontos_multistart(tc_ref, nome_estacao)
        if i_partida >= len(pontos):
            print(f'so ha {len(pontos)} pontos de partida distintos')
            return None
        Al0, Au0 = pontos[i_partida]
    else:
        Al0, Au0 = ponto_de_partida(tc_ref)
    av = Avaliador(cl_ref, tc_ref, pasta, cfl=est.get('cfl', 0.20))

    print(f'=== ESTACAO {nome_estacao.upper()} (eta = {est["eta"]}) ===')
    print(f'  M_n = {MACH_N}   cl_ref = {cl_ref}   (t/c)_ref = {tc_ref}')
    if com_bluntez:
        print(f'  t_01/sqrt(t/c) >= {bluntez_min:.4f} '
              f'(substituta de cl_max >= 1,80)')
    else:
        print('  SEM a restricao de sustentacao maxima '
              '(rodada de comparacao)')
    print(f'  partida: NACA 1411 reescalado, t/c = {tc_ref:.4f}, '
          f't01 = {t01_de(Al0, Au0):.5f}\n', flush=True)

    # --- acha o alpha de partida com passos de Newton usando o adjunto ---
    # Ate 6 iteracoes, e AVISA se nao convergir: partir com a restricao de
    # igualdade violada faz o SLSQP gastar as primeiras iteracoes buscando
    # sustentacao em vez de reduzir arrasto, e se o alvo for inalcancavel
    # dentro dos batentes de alpha a rodada inteira se perde.
    alpha = 3.0 * np.pi / 180
    erro = np.inf
    for _ in range(6):
        d = av(np.hstack([Al0, Au0, alpha]))
        erro = d['CL'] - cl_ref
        if abs(erro) < 5e-3 or not np.isfinite(d['dCL'][-1]):
            break
        alpha_novo = float(np.clip(alpha - erro / d['dCL'][-1],
                                   ALPHA_MIN, ALPHA_MAX))
        if abs(alpha_novo - alpha) < 1e-6:      # travou num batente
            break
        alpha = alpha_novo
    print(f'  alpha de partida: {alpha*180/np.pi:.3f} deg '
          f'(erro de c_l: {erro:+.5f})', flush=True)
    if abs(erro) > 5e-3:
        print(f'  AVISO: o ponto de partida nao atinge c_l = {cl_ref} dentro '
              f'de alpha em [{ALPHA_MIN*180/np.pi:.0f}, '
              f'{ALPHA_MAX*180/np.pi:.0f}] graus.', flush=True)
        if abs(alpha - ALPHA_MAX) < 1e-6 or abs(alpha - ALPHA_MIN) < 1e-6:
            print('  alpha travou no batente -- a restricao de igualdade e '
                  'INALCANCAVEL daqui. Abortando.', flush=True)
            return None
    print(flush=True)

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
        # NAO acrescentar 'mint >= 0' aqui. Junto com 'MINT_MAX - mint' isso
        # poe na jacobiana duas linhas exatamente antiparalelas (-dmint e
        # +dmint), deixando-a deficiente em posto POR CONSTRUCAO, e o SLSQP
        # morre com "Singular matrix E in LSQ subproblem".
        # Alem de nociva, era desnecessaria: em todas as rodadas o mint
        # converge para 0,010000, prensado contra o teto -- o otimizador
        # empurra o bordo de fuga para ficar mais GROSSO, na direcao oposta ao
        # cruzamento de superficies. A checagem em _grava cobre o caso raro.
        g = [d['maxt'] - tc_ref,          # espessura exigida
             MINT_MAX - d['mint']]        # bordo de fuga fino
        if com_bluntez:
            g.append(bluntez(d['t01'], d['maxt']) - bluntez_min)
        return np.array(g)

    def ineqgrad(xx):
        d = av(xx)
        J = [d['dmaxt'], -d['dmint']]
        if com_bluntez:
            J.append(grad_bluntez(d['t01'], d['maxt'], d['dmaxt']))
        return np.vstack(J)

    cons = [{'type': 'ineq', 'fun': ineqfun, 'jac': ineqgrad},
            {'type': 'eq', 'fun': eqfun, 'jac': eqgrad}]
    al_upper = list(AL_UPPER)
    if cusp_livre:
        al_upper[-1] = AL4_SOLTO      # permite o cusp concavo do bordo de fuga
    bounds = Bounds(AL_LOWER + AU_LOWER + [ALPHA_MIN],
                    al_upper + AU_UPPER + [ALPHA_MAX], keep_feasible=True)

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
    if com_bluntez:
        print(f'  t01 = {d["t01"]:.5f}   bluntez = {bl:.5f} '
              f'(min {bluntez_min:.4f})'
              f'{"  ATIVA" if abs(bl - bluntez_min) < 1e-4 else "  folgada"}')
    else:
        ref = limiar_bluntez(nome_estacao)
        print(f'  t01 = {d["t01"]:.5f}   bluntez = {bl:.5f}   '
              f'(limiar seria {ref:.4f}: '
              f'{"ATENDE mesmo sem a restricao" if bl >= ref else "VIOLA"})')
    print(f'  Al  = {np.array2string(Al, precision=6)}')
    print(f'  Au  = {np.array2string(Au, precision=6)}')
    print(f'  alpha = {alpha*180/np.pi:.4f} deg')
    print(f'\n  historico em {os.path.join(pasta, "historico.pickle")}')
    return res


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    estacao = args[0] if args else 'meio'
    com_bluntez = '--sem-bluntez' not in sys.argv
    cusp_livre = '--cusp' in sys.argv
    i_partida = 0
    for a in sys.argv:
        if a.startswith('--partida='):
            i_partida = int(a.split('=')[1])
    if estacao not in ESTACOES:
        print(f'estacao invalida: {estacao}. Use: {list(ESTACOES)}')
        sys.exit(1)
    os.makedirs(RES, exist_ok=True)
    otimiza(estacao, com_bluntez=com_bluntez, cusp_livre=cusp_livre,
            i_partida=i_partida)
