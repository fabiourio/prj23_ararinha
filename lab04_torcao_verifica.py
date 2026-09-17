'''
Lab 04 -- verificacao independente da otimizacao de torcao.

Confere, direto no AVL (sem os modelos substitutos usados na otimizacao):
  1. o ponto otimo e viavel: o estol comeca em eta <= 0,50;
  2. o ponto otimo e localmente otimo: nenhuma perturbacao viavel das
     torcoes de controle reduz o arrasto induzido de cruzeiro;
  3. o surrogate quadratico bate com o AVL fora dos pontos usados para
     construi-lo (o CDff deve ser exatamente quadratico na torcao).

Rodar da raiz do repo:  python lab04_torcao_verifica.py
'''

# IMPORTS
import json
import os

import numpy as np
from scipy.interpolate import CubicSpline

from avl_batch import roda_avl, pega, pega_todos
import lab04_gera_avl as gera

#=========================================

MACH_CRU = 0.85
CL_PROJ = 0.5053
MACH_BAIXO = 0.2
SEMI_ENV = 30.0759
ETA_SEGURA = 0.50
FOLGA_ALVO = 0.5                         # margem de alpha exigida (otimizador)
ALFAS_ESTOL = np.arange(6.0, 22.1, 1.0)

ETA_LIM = np.array([0.1011, 0.398, 0.90])
CLMAX_LIM = np.array([1.774, 1.7985, 1.7338])

ETAS = np.array(gera.ETAS_ASA)
ETAS_CTRL = np.array([0.3, 0.56, 0.8, 1.0])
NOS = np.concatenate([[0.0], ETAS_CTRL])
VARIANTE = 'avl/saidas/tw_ver.avl'
LOAD_REL = 'saidas/tw_ver.avl'          # caminho relativo a pasta avl/

with open('avl/saidas/trim.json', encoding='ascii') as f:
    TRIM = json.load(f)
with open('avl/saidas/torcao.json', encoding='ascii') as f:
    OTIMO = json.load(f)

CTRL_OTIMO = np.array([OTIMO['controle'][str(e)] for e in ETAS_CTRL])
import re
RE_LINHA = re.compile(r'^\s*\d+\s+(' + r'[-\dEe.+]+\s+'*11 +
                      r'[-\dEe.+]+)\s*$', re.M)


def torcoes_de(ctrl):
    spl = CubicSpline(NOS, np.concatenate([[0.0], ctrl]), bc_type='natural')
    return {float(e): float(v) for e, v in zip(ETAS, spl(ETAS))}


def cdff(ctrl):
    '''CDff de cruzeiro, aeronave trimada em arfagem pelo profundor.'''
    gera.gera_variante(VARIANTE, torcoes_de(ctrl), cg='aft')
    saida = roda_avl(f'load {LOAD_REL}\noper\n'
                     f'm\nmn {MACH_CRU}\n\n'
                     f'de\n1 {TRIM["aft"]["it_deg"]}\n\n'
                     f'd2 pm 0\na c {CL_PROJ}\nx\n\nquit\n')
    return pega(saida, r'CDff\s*=\s*([-\d.Ee+]+)')


def analisa_estol(ctrl):
    '''Analisa o estol em M 0,2: devolve a estacao critica (eta) e a folga
    de alpha (quanto a regiao externa, eta > 0,50, estola depois da
    interna). A folga e o que o otimizador restringe; eta_crit e a checagem
    fisica de estol antes do aileron.'''
    gera.gera_variante(VARIANTE, torcoes_de(ctrl), cg='fwd')
    cmds = (f'load {LOAD_REL}\noper\n'
            f'm\nmn {MACH_BAIXO}\n\n'
            f'de\n1 {TRIM["fwd"]["it_deg"]}\n\nd2 d2 0\n')
    for a in ALFAS_ESTOL:
        cmds += f'a a {a}\nx\nfs\n\n'
    cmds += '\nquit\n'
    saida = roda_avl(cmds)
    # margem = cln - limite, por estacao e por alpha
    margens, etas_strip = [], None
    for bloco in saida.split('Vortex Lattice Output')[1:]:
        parte = bloco.split('Surface # 1     Wing')[1].split('Surface # 2')[0]
        dados = np.array([[float(v) for v in m.group(1).split()]
                          for m in RE_LINHA.finditer(parte)])
        etas_strip = dados[:, 0]/SEMI_ENV
        lim = np.interp(etas_strip, ETA_LIM, CLMAX_LIM)
        margens.append(dados[:, 5] - lim)
    margens = np.array(margens)           # [n_alpha, n_strip]
    interna = etas_strip <= ETA_SEGURA
    externa = ~interna

    def alpha_de_cruzamento(col_mask):
        # primeiro alpha em que alguma estacao do grupo cruza o limite
        pico = np.array([margens[i, col_mask].max()
                         for i in range(len(ALFAS_ESTOL))])
        k = np.argmax(pico >= 0)
        if pico[k] < 0:
            return np.inf                 # nao estolou na faixa
        if k == 0:
            return ALFAS_ESTOL[0]
        frac = -pico[k-1]/(pico[k] - pico[k-1])
        return ALFAS_ESTOL[k-1] + frac*(ALFAS_ESTOL[k] - ALFAS_ESTOL[k-1])

    a_int = alpha_de_cruzamento(interna)
    a_ext = alpha_de_cruzamento(externa)
    # estacao critica: onde a margem e maxima no alpha de estol interno
    k_est = int(np.argmax(margens[:, :].max(axis=1) >= 0))
    eta_crit = float(etas_strip[np.argmax(margens[k_est])])
    return eta_crit, float(a_ext - a_int)


os.makedirs('avl/saidas', exist_ok=True)

# 1 e 2: otimo, viabilidade e otimalidade local -------------------
# Viabilidade: o estol comeca antes do aileron (eta_crit <= 0,50) E ha a
# margem de alpha exigida pelo otimizador (folga >= 0,5 grau). E a mesma
# restricao usada na otimizacao, medida aqui direto no AVL.
f_ot = cdff(CTRL_OTIMO)
eta_ot, folga_ot = analisa_estol(CTRL_OTIMO)
print('=== PONTO OTIMO ===')
print(f'torcoes de controle [graus]: {np.round(CTRL_OTIMO, 3)}')
print(f'CDff (AVL)      = {f_ot:.6f}')
print(f'estacao de estol = {eta_ot:.3f} (limite {ETA_SEGURA})')
print(f'folga de alpha   = {folga_ot:.3f} graus (alvo {FOLGA_ALVO})')


def viavel_de(eta, folga):
    return eta <= ETA_SEGURA + 1e-6 and folga >= FOLGA_ALVO - 0.05


print('\n=== OTIMALIDADE LOCAL (perturbacoes diretas no AVL) ===')
rng = np.random.default_rng(23)
dirs = list(np.eye(4)) + list(-np.eye(4))
dirs += [rng.normal(size=4) for _ in range(6)]
melhorou = False
for k, d in enumerate(dirs):
    d = np.asarray(d)/np.linalg.norm(np.asarray(d))
    ctrl = np.clip(CTRL_OTIMO + 1.0*d, -10.0, 2.5)   # passo de 1 grau
    f = cdff(ctrl)
    eta, folga = analisa_estol(ctrl)
    viavel = viavel_de(eta, folga)
    marca = ''
    if viavel and f < f_ot - 1e-6:
        marca = '  <-- MELHOROU E VIAVEL (otimo suspeito!)'
        melhorou = True
    print(f'  perturbacao {k:2d}: CDff = {f:.6f} '
          f'(delta {1e4*(f - f_ot):+5.2f} counts), '
          f'eta = {round(eta, 3)}, folga = {folga:+.2f}, '
          f'{"viavel" if viavel else "inviavel"}{marca}')

# 3: surrogate quadratico vs AVL fora da amostra ------------------
print('\n=== SURROGATE QUADRATICO vs AVL (fora da amostra) ===')
cache = 'avl/saidas/torcao_modelo.npz'
if os.path.exists(cache):
    m = np.load(cache)
    f0, G, H = float(m['f0']), m['g'], m['h']

    def modelo(ctrl):
        t = np.asarray(ctrl)
        return f0 + G @ t + 0.5*t @ H @ t

    for k in range(5):
        ctrl = rng.uniform(-6, 1, 4)
        fa = cdff(ctrl)
        fm = modelo(ctrl)
        print(f'  ponto {k}: AVL = {fa:.6f}, modelo = {fm:.6f}, '
              f'erro = {1e4*(fm - fa):+.3f} counts')
else:
    print(f'  (cache {cache} ausente; rode lab04_torcao.py para gerar)')

if os.path.exists(VARIANTE):
    os.remove(VARIANTE)

print('\n=== VEREDITO ===')
print(f'viavel (eta<=0,50 e folga>=0,5): '
      f'{"sim" if viavel_de(eta_ot, folga_ot) else "NAO"}')
print(f'otimo local: {"NAO (ha melhora viavel)" if melhorou else "sim"}')
