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

import json
import os
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

from aeronave import aeronave
from analises import it_para_de_zero
from avl_run import caso
from estilo import plt, AZUL, LARANJA, VERDE, TINTA2, CINZA
from estol import modelo_estol, estol
from gera_avl import escreve_avl, ETAS, AQUI
from sombra import alfa_saida

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


# --------------------------------------------------------------------------
# Etapa 2 -- otimizacao de torcao na posicao escolhida da EH, com conferencia
# de independencia em dois Lc_h vizinhos (rodados em paralelo, cada um com o
# seu proprio Avaliador/arquivo .avl, ver docstring do modulo).
#
# it fixo: calculado uma unica vez, no modelo SEM torcao (asa plana), e
# mantido constante durante toda a otimizacao -- refazer o trim a cada
# avaliacao mudaria o tamanho/sentido do gradiente numerico do SLSQP com a
# torcao (o d2 pm 0 arrasta o it junto) e tornaria a superficie de CDff(x)
# menos suave. A verificacao final (main()) retrima o modelo JA torcido
# (it proprio) para conferir que a diferenca e desprezivel.
# --------------------------------------------------------------------------

FOLGA_SOMBRA_MIN = 2.0   # [graus] exigido pela etapa 2 (ver varredura_eh.py)


def _projeto(lc):
    '''
    Otimiza a torcao da asa para um Lc_h (modelo completo, CG traseiro).
    Roda dentro de um Pool: devolve so dados picklaveis (o Avaliador guarda
    strings da saida do AVL, que nao interessam fora do processo).
    '''
    av = aeronave(lc)
    base = escreve_avl(av, f'resultados/_tmp/base_{lc:.4f}.avl', cg='aft')
    it, _ = it_para_de_zero(base, av['M'], av['CL'])
    a = Avaliador(av, modelo='completo', it=it, com_estol=True, nome=f'ot_{lc:.4f}')
    x0 = np.zeros(len(ETAS) - 1)
    base_res = a.avalia(x0)
    r = otimiza(a, x0)
    fim = a.avalia(r.x)
    viavel = bool(np.all(fim['g'] >= -1e-3))
    tol_limite = 1e-3
    nos_limites = [i for i, v in enumerate(r.x)
                   if min(abs(v - LIMITES[0]), abs(v - LIMITES[1])) < tol_limite]
    return {
        'Lc_h': lc, 'it': it,
        'twist': [float(v) for v in r.x],
        'CDff': 1e4*fim['CDff'], 'e': fim['e'],
        'CDff_base': 1e4*base_res['CDff'], 'e_base': base_res['e'],
        'g_min': float(np.min(fim['g'])),
        'sucesso': bool(r.success), 'mensagem': str(r.message),
        'n_rodadas': a.n_rodadas, 'viavel': viavel, 'nos_limites': nos_limites,
    }


def _figura_torcao(resultados, e_aft, nome):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0))

    cores = [AZUL, LARANJA, VERDE]
    for res, cor in zip(resultados, cores):
        tw = np.concatenate([[0.0], res['twist']])
        rot = f'Lc_h {res["Lc_h"]:.4f}' + (' (escolhido)' if res is resultados[0] else '')
        ax1.plot(ETAS, tw, '-o', color=cor, label=rot, markersize=3.0)
    ax1.set_xlabel('η = 2y/b')
    ax1.set_ylabel('torção, Ainc [°]')
    ax1.set_title('Torção ótima: independência da posição da EH')
    ax1.legend(fontsize=7)

    eta, cl_norm, clmax = e_aft['eta'], e_aft['cl_norm'], e_aft['clmax']
    excl = np.isin(eta, e_aft['excluidas'])
    ax2.plot(eta, clmax, color=TINTA2, lw=1.2, ls='--', label='cl_max do perfil')
    ax2.plot(eta[~excl], cl_norm[~excl], color=AZUL, label='cl_norm (válidas)')
    if np.any(excl):
        ax2.scatter(eta[excl], cl_norm[excl], color=CINZA, marker='x', zorder=3,
                   label='excluída (raiz/junção)')
    ax2.axvline(ETA_ESTOL, color=TINTA2, lw=0.8, ls=':')
    ax2.set_xlabel('η = 2y/b')
    ax2.set_ylabel('cl_norm')
    ax2.set_title(f'Estol da asa torcida (CG traseiro): η_crit = {e_aft["eta_crit"]:.3f}')
    ax2.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, 'resultados', nome))


def _figura_carga(av, eta_b, carga_b, eta_a, carga_a, nome):
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    eta_e = np.linspace(0.0, 1.0, 201)
    ax.plot(eta_e, elipse_normalizada(eta_e), color=TINTA2, lw=1.2, ls='--', label='elíptica')
    ax.plot(eta_b, carga_b, color=LARANJA, label='sem torção')
    ax.plot(eta_a, carga_a, color=AZUL, label='torção ótima')
    ax.set_xlabel('η = 2y/b')
    ax.set_ylabel('c·cl / (4 Sref CL / π b)')
    ax.set_title('Carga no cruzeiro (informativo: com EH/fuselagem o\n'
                  'ótimo de CDff não é a elíptica)')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, 'resultados', nome))


def main():
    with open(os.path.join(AQUI, 'resultados', 'lc_h_escolhido.json')) as f:
        Lc_h = json.load(f)['Lc_h']
    lcs = [Lc_h, Lc_h - 0.2, Lc_h - 0.4]
    print(f'Lc_h: {lcs[0]:.4f} (escolhido), vizinhos {lcs[1]:.4f} e {lcs[2]:.4f} '
          '(conferencia de independencia)')

    print('medindo o tempo de uma avaliacao (3 rodadas do AVL, malha adotada)...')
    t0 = time.time()
    av_t = aeronave(lcs[0])
    base_t = escreve_avl(av_t, 'resultados/_tmp/timing.avl', cg='aft')
    it_t, _ = it_para_de_zero(base_t, av_t['M'], av_t['CL'])
    Avaliador(av_t, modelo='completo', it=it_t, com_estol=True,
              nome='timing').objetivo(np.zeros(len(ETAS) - 1))
    t_avaliacao = time.time() - t0
    # ordem de grandeza observada em verifica_eliptica.py (asa_limpa, mesmo
    # numero de variaveis, so sem a restricao de estol): ~220 avaliacoes ate
    # o SLSQP convergir.
    evals_estimados = 220
    eta_min_caso = t_avaliacao*evals_estimados/60
    print(f'  1 avaliacao = {t_avaliacao:.1f} s  ->  ETA por caso ~{eta_min_caso:.0f} min '
          '(parede: os 3 casos rodam em paralelo, Pool(3))')
    if eta_min_caso > 120:
        print('AVISO: ETA > 2 h -- reportar antes de continuar (NEEDS_CONTEXT)')
        return

    print('otimizando os 3 casos em paralelo...')
    t0 = time.time()
    with Pool(3) as p:
        resultados = p.map(_projeto, lcs)
    dt_total = time.time() - t0
    print(f'  levou {dt_total/60:.1f} min de parede')

    principal, viz1, viz2 = resultados
    for res in resultados:
        marca = []
        if not res['viavel']:
            marca.append('INVIAVEL (g_min=%.3f)' % res['g_min'])
        if res['nos_limites']:
            marca.append(f'limites em {res["nos_limites"]}')
        print(f'  Lc_h {res["Lc_h"]:.4f}: it {res["it"]:.3f}  CDff '
              f'{res["CDff_base"]:.2f} -> {res["CDff"]:.2f} count  '
              f'e {res["e_base"]:.4f} -> {res["e"]:.4f}  {res["n_rodadas"]} rodadas  '
              f'"{res["mensagem"]}"' + (('  [' + '; '.join(marca) + ']') if marca else ''))

    # 3. fwd.avl e aft.avl finais, na posicao escolhida, com a torcao otima
    av = aeronave(principal['Lc_h'])
    torcao = np.concatenate([[0.0], principal['twist']])
    escreve_avl(av, 'aft.avl', cg='aft', torcao=torcao)
    escreve_avl(av, 'fwd.avl', cg='fwd', torcao=torcao)
    print('escritos aft.avl e fwd.avl (torção ótima, malha adotada)')

    # estol nos dois CG, cada um trimado com o seu proprio it de cruzeiro
    it_aft, _ = it_para_de_zero('aft.avl', av['M'], av['CL'])
    it_fwd, _ = it_para_de_zero('fwd.avl', av['M'], av['CL'])
    e_aft = estol('aft.avl', it=it_aft, trim=True, av=av)
    e_fwd = estol('fwd.avl', it=it_fwd, trim=True, av=av)
    stall = {}
    for tag, e in (('aft', e_aft), ('fwd', e_fwd)):
        ok = e['eta_crit'] <= ETA_ESTOL
        print(f'  estol {tag} (it {(it_aft if tag == "aft" else it_fwd):.3f}): '
              f'alfa {e["alfa"]:.2f}°  CLmax {e["CL"]:.3f}  de {e["de"]:.2f}°  '
              f'eta_crit {e["eta_crit"]:.3f}  {"ok" if ok else "FALHOU"}')
        stall[tag] = {'alfa': e['alfa'], 'CL': e['CL'], 'de': e['de'],
                      'eta_crit': e['eta_crit'], 'ok': ok}

    # sombreamento: folga entre o alfa de estol (o menor dos dois CG) e o
    # alfa em que a EH sai da sombra da asa.
    alfa_estol_min = min(e_aft['alfa'], e_fwd['alfa'])
    a_saida = alfa_saida(av)
    folga_sombra = alfa_estol_min - a_saida
    ok_sombra = folga_sombra >= FOLGA_SOMBRA_MIN
    print(f'  sombra: alfa_estol_min {alfa_estol_min:.2f}°  alfa_saida {a_saida:.2f}°  '
          f'folga {folga_sombra:.2f}°  {"ok" if ok_sombra else "FAILED"}')
    if not ok_sombra:
        print('  AVISO: folga da sombra abaixo do exigido -- registrado, decisao do '
              'controlador do projeto')

    # CDff/e antes/depois, cada um com o seu proprio it que zera o profundor
    # (o "antes" reusa o mesmo it fixo da otimizacao; o "depois" retrima o
    # modelo JA torcido, para conferir que a diferenca de it e desprezivel).
    base_arq = f'resultados/_tmp/base_{principal["Lc_h"]:.4f}.avl'
    it_antes, r_antes = it_para_de_zero(base_arq, av['M'], av['CL'])
    it_depois, r_depois = it_para_de_zero('aft.avl', av['M'], av['CL'])
    print(f'  CDff (CG traseiro): {1e4*r_antes["CDff"]:.2f} -> {1e4*r_depois["CDff"]:.2f} '
          f'count  (it {it_antes:.3f} -> {it_depois:.3f})')
    print(f'  e    (CG traseiro): {r_antes["e"]:.4f} -> {r_depois["e"]:.4f}')

    # independencia: diferenca de torcao com os vizinhos
    vizinhos = []
    for viz in (viz1, viz2):
        d = float(np.max(np.abs(np.array(viz['twist']) - np.array(principal['twist']))))
        print(f'  independencia Lc_h {viz["Lc_h"]:.4f}: max|Δtorção| = {d:.3f}°  '
              f'(CDff = {viz["CDff"]:.2f} count)')
        vizinhos.append({'Lc_h': viz['Lc_h'], 'twist': viz['twist'],
                         'CDff': viz['CDff'], 'diff_max_deg': d})

    # cargas no cruzeiro (informativo): antes/depois, com faixas
    r_carga_antes = caso(base_arq, av['M'], cl=av['CL'], it=it_antes, trim=True, faixas=True)
    r_carga_depois = caso('aft.avl', av['M'], cl=av['CL'], it=it_depois, trim=True, faixas=True)
    eta_b, carga_b = carga_normalizada(av, r_carga_antes)
    eta_a, carga_a = carga_normalizada(av, r_carga_depois)

    _figura_torcao(resultados, e_aft, 'otimizacao_torcao.png')
    _figura_carga(av, eta_b, carga_b, eta_a, carga_a, 'carga_torcao.png')

    saida = {
        'Lc_h': principal['Lc_h'], 'it': principal['it'],
        'etas': list(ETAS), 'twist': principal['twist'],
        'CDff_antes': 1e4*r_antes['CDff'], 'CDff_depois': 1e4*r_depois['CDff'],
        'e_antes': r_antes['e'], 'e_depois': r_depois['e'],
        'it_antes': it_antes, 'it_depois': it_depois,
        'estol': stall,
        'sombra': {'alfa_estol_min': alfa_estol_min, 'alfa_saida': a_saida,
                   'folga': folga_sombra, 'ok': ok_sombra},
        'vizinhos': vizinhos,
        'sucesso_slsqp': principal['sucesso'], 'mensagem_slsqp': principal['mensagem'],
        'nos_limites': principal['nos_limites'], 'g_min': principal['g_min'],
        'n_rodadas': principal['n_rodadas'], 'tempo_total_min': dt_total/60,
    }
    with open(os.path.join(AQUI, 'resultados', 'torcao_otimizada.json'), 'w') as f:
        json.dump(saida, f, indent=2, ensure_ascii=False)
    print('salvo resultados/torcao_otimizada.json')


if __name__ == '__main__':
    main()
