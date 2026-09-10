'''
DOE de cl_max no XFoil -- Lab 03, PRJ-23, equipe Ararinha.

PERGUNTA QUE ESTE DOE RESPONDE
------------------------------
A otimizacao transonica (eulerblock, nao-viscoso) nao enxerga cl_max. Nada
impede o otimizador de afiar o bordo de ataque para reduzir o arrasto de onda
e, no caminho, destruir a sustentacao maxima do perfil -- que o designTool
assumiu como clmax_w = 1.8 e da qual dependem CLmaxTO, as distancias de
decolagem e pouso e, por consequencia, a area da asa da aeronave.

Como cl_max do XFoil e nao-diferenciavel (vem de varrer alpha ate o solver
viscoso divergir), ele nao serve como g(x) num otimizador de gradiente.
A saida e uma RESTRICAO SUBSTITUTA: um descritor geometrico analitico e suave
que preveja cl_max bem o bastante para ser imposto como restricao barata.

Este DOE nao pressupoe QUAL descritor serve. Ele varre o espaco de projeto,
mede cl_max no XFoil e registra doze candidatos (ver descritores.py) para que
a analise decida -- inclusive a possibilidade de que nenhum sirva sozinho e
seja preciso combinar dois ou tres.

DOIS ESTUDOS
------------
A) Corte controlado: varia so o nariz (Au[0], Al[0]) em torno do NACA 1411.
   Isola o efeito do raio de bordo de ataque. E a figura didatica.

B) LHS nas 8 variaveis CST: mede se os descritores continuam predizendo
   quando espessura, arqueamento e carregamento traseiro tambem variam.
   E o experimento que decide a restricao.

Rodar com:  python doe_clmax_xfoil.py
Saida:      resultados/doe_clmax_corte.csv
            resultados/doe_clmax_lhs.csv
'''

import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.stats import qmc

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import descritores as dsc
import xfoil_runner as xr

RES = os.path.join(AQUI, 'resultados')

# --- condicao de analise -----------------------------------------------------
# Reynolds e Mach da DECOLAGEM na corda media aerodinamica da aeronave B.
# Ver ponto_projeto: V2 = 90.59 m/s, c_MAC = 6.899 m, nivel do mar.
RE_DECOLAGEM = 4.2e7
MACH_DECOLAGEM = 0.266
ALPHA_SEQ = (0.0, 30.0, 0.5)

# Dois niveis, ver secao 1.1 do documento de projeto:
#   ALVO       -- valor assumido pelo designTool (standard_airplane.py:215).
#                 Manter isso deixa a aeronave B exatamente como o Lab 02 a
#                 deixou.
#   VIABILIDADE -- abaixo disso o empuxo requerido para a decolagem em 2900 m
#                 excede o Tmax dos motores (T0req/T0 = 1.000 em clmax_w=1.20).
#                 Verificado varrendo clmax_w na aeronave B; a restricao de
#                 pouso nunca fica ativa -- quem limita e o empuxo.
CLMAX_ALVO = 1.8
CLMAX_VIABILIDADE = 1.2

# NACA 1411 -- ponto de partida da otimizacao do roteiro
AU_REF = np.array([0.16146332, 0.18349204, 0.14126241, 0.18194397])
AL_REF = np.array([-0.1489439, -0.10330027, -0.10305128, -0.10514982])

# --- espaco de projeto do LHS ------------------------------------------------
# Faixas escolhidas para cobrir perfis transonicos plausiveis. Al[3] pode ser
# POSITIVO: e assim que se obtem o carregamento traseiro (cusp) dos perfis
# supercriticos -- o RAE2822 tem Al[-1] = +0.052.
LIM_AU = [(0.06, 0.26), (0.06, 0.30), (0.06, 0.30), (0.06, 0.30)]
LIM_AL = [(-0.26, -0.06), (-0.30, -0.02), (-0.30, -0.02), (-0.25, 0.15)]

# perfis inviaveis ou fora de proposito sao descartados antes de gastar XFoil
T_MIN_ACEITAVEL = 0.002
T_MAX_FAIXA = (0.07, 0.20)


def avalia(caso):
    '''Avalia um perfil: descritores geometricos + cl_max no XFoil.'''
    idx, Au, Al = caso[:3]
    faixa_t = caso[3] if len(caso) > 3 else T_MAX_FAIXA
    Au, Al = np.asarray(Au), np.asarray(Al)

    try:
        d = dsc.descritores(Au, Al)
    except Exception as exc:                                   # noqa: BLE001
        return {'idx': idx, 'status': f'geom_erro:{type(exc).__name__}'}

    # filtros de viabilidade -- evitam rodar XFoil em geometria degenerada
    if d['t_min'] < T_MIN_ACEITAVEL:
        return {'idx': idx, 'status': 'descartado_t_min', **d}
    if not (faixa_t[0] <= d['t_max'] <= faixa_t[1]):
        return {'idx': idx, 'status': 'descartado_t_max', **d}

    try:
        r = xr.clmax_cst(Au, Al, Re=RE_DECOLAGEM, Mach=MACH_DECOLAGEM,
                         alpha_seq=ALPHA_SEQ, timeout=600)
    except Exception as exc:                                   # noqa: BLE001
        return {'idx': idx, 'status': f'xfoil_erro:{type(exc).__name__}', **d}

    linha = {'idx': idx, 'status': 'ok', **d,
             'clmax': r['clmax'], 'alpha_clmax': r['alpha_clmax'],
             'situacao': r['situacao'], 'confiavel': int(r['confiavel']),
             'n_pontos': r['n_pontos'], 'n_buracos': r['n_buracos'],
             'alpha_max_conv': r['alpha_max_convergido'],
             'incl_linear': r['incl_linear'], 'incl_fim': r['incl_fim'],
             'razao_incl': r['razao_incl'],
             'tentativas': r.get('tentativas', 1),
             'passo_final': r.get('passo_final', ALPHA_SEQ[2])}
    for ii in range(len(Au)):
        linha[f'Au{ii}'] = float(Au[ii])
    for ii in range(len(Al)):
        linha[f'Al{ii}'] = float(Al[ii])
    return linha


def casos_corte(n=17):
    '''Estudo A: escala o nariz mantendo o resto do NACA 1411.'''
    casos = []
    for ii, k in enumerate(np.linspace(0.40, 1.55, n)):
        Au = AU_REF.copy()
        Al = AL_REF.copy()
        Au[0] *= k
        Al[0] *= k
        casos.append((ii, Au, Al))
    return casos


def casos_lhs(n=256, semente=23):
    '''Estudo B: amostragem por hipercubo latino nas 8 variaveis.'''
    limites = LIM_AU + LIM_AL
    amostra = qmc.LatinHypercube(d=len(limites), seed=semente).random(n)
    baixo = np.array([a for a, _ in limites])
    alto = np.array([b for _, b in limites])
    X = baixo + amostra * (alto - baixo)
    return [(ii, X[ii, :4], X[ii, 4:], T_MAX_FAIXA) for ii in range(n)]


# Caixa deslocada para perfis GROSSOS. O estudo B foi so ate t/c = 0,194, mas
# a estacao da raiz precisa de (t/c)_n = 0,218 -- ficaria fora da amostra, e
# extrapolar o limiar para la seria chute.
LIM_AU_GROSSO = [(0.12, 0.42), (0.12, 0.46), (0.12, 0.46), (0.12, 0.46)]
LIM_AL_GROSSO = [(-0.42, -0.12), (-0.46, -0.10), (-0.46, -0.10), (-0.40, 0.05)]
T_MAX_FAIXA_GROSSO = (0.16, 0.28)


def casos_grossos(n=96, semente=2309):
    '''Estudo C: cobre a faixa de espessura da estacao da raiz.'''
    limites = LIM_AU_GROSSO + LIM_AL_GROSSO
    amostra = qmc.LatinHypercube(d=len(limites), seed=semente).random(n)
    baixo = np.array([a for a, _ in limites])
    alto = np.array([b for _, b in limites])
    X = baixo + amostra * (alto - baixo)
    return [(ii, X[ii, :4], X[ii, 4:], T_MAX_FAIXA_GROSSO) for ii in range(n)]


def casos_corte_com_faixa():
    return [(i, au, al, T_MAX_FAIXA) for i, au, al in casos_corte()]


def escreve_csv(linhas, caminho):
    colunas = []
    for l in linhas:
        for k in l:
            if k not in colunas:
                colunas.append(k)
    with open(caminho, 'w') as fid:
        fid.write(','.join(colunas) + '\n')
        for l in linhas:
            fid.write(','.join(
                '' if k not in l else
                (f'{l[k]:.6g}' if isinstance(l[k], (int, float, np.floating))
                 else str(l[k]))
                for k in colunas) + '\n')


def roda(casos, nome, n_proc):
    print(f'\n=== {nome}: {len(casos)} perfis, {n_proc} processos ===', flush=True)
    t0 = time.time()
    with Pool(n_proc) as pool:
        linhas = []
        for ii, linha in enumerate(pool.imap_unordered(avalia, casos), 1):
            linhas.append(linha)
            if ii % 20 == 0 or ii == len(casos):
                dt = time.time() - t0
                print(f'  {ii}/{len(casos)}  ({dt:.0f} s, '
                      f'{dt/ii:.1f} s/perfil)', flush=True)

    linhas.sort(key=lambda l: l['idx'])
    caminho = os.path.join(RES, f'doe_clmax_{nome}.csv')
    escreve_csv(linhas, caminho)

    ok = [l for l in linhas if l.get('status') == 'ok']
    conf = [l for l in ok if l.get('confiavel')]
    print(f'  gravado em {caminho}')
    print(f'  {len(ok)}/{len(linhas)} avaliados, {len(conf)} com estol confiavel')
    if conf:
        cl = np.array([l['clmax'] for l in conf])
        print(f'  cl_max: min={cl.min():.3f}  mediana={np.median(cl):.3f}  '
              f'max={cl.max():.3f}')
        print(f'  fracao com cl_max >= {CLMAX_ALVO}: '
              f'{np.mean(cl >= CLMAX_ALVO)*100:.1f}%')
    return linhas


if __name__ == '__main__':
    os.makedirs(RES, exist_ok=True)
    n_proc = max(1, min(6, (os.cpu_count() or 4) - 2))

    quais = sys.argv[1:] or ['corte', 'lhs', 'grossos']
    if 'corte' in quais:
        roda(casos_corte_com_faixa(), 'corte', n_proc)
    if 'lhs' in quais:
        roda(casos_lhs(), 'lhs', n_proc)
    if 'grossos' in quais:
        roda(casos_grossos(), 'grossos', n_proc)
