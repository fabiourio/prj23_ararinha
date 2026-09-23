'''
Metodo da secao critica.

O cl_norm de cada faixa da asa e linear em alfa no VLM (inclusive com o
profundor de compensacao, que tambem varia linearmente). Duas rodadas dao o
alfa em que cada faixa alcanca o clmax do seu perfil; a menor define o
estol. Uma terceira rodada, no alfa encontrado, confere o resultado.

clmax dos perfis do Lab 03 (XFoil corrigido, otimizacao_aerofolio/campanha/
resultados/revisao_plano_clmax.csv, plano "normal"), interpolado entre as
estacoes de referencia e constante fora delas.

Faixa da juncao com o winglet: a ultima faixa da asa, na juncao com o
winglet, tem um vertice de canto (o espacamento cosseno da malha a deixa
~1/ns^2 mais estreita que as vizinhas). Esse vertice e uma singularidade
de vórtice de canto: o cl_norm dessa faixa cresce sem limite conforme ela
estreita com o refino da malha (a area -> 0, entao CL/CDff nao sao
afetados), e um criterio de estol por secao nao vale num canto assim (a
aeronave real tem um filete na juncao). Por isso a faixa e as vizinhas
proximas o bastante do canto (a `DIST_JUNCAO` do winglet) sao excluidas da
busca da secao critica, tanto no ajuste linear quanto na conferencia.
DIST_JUNCAO = ct/4 (fracao da corda de ponta da asa, constante) e uma
distancia pequena mas suficiente para excluir so a(s) faixa(s) degeneradas
perto do canto, sem comer estacoes fisicamente relevantes; foi calibrada
comparando o alfa de estol entre malhas com asa_ns = 20 e 60 (ver teste em
tests/test_analises.py). Vale tambem para empenagem=False (asa isolada,
sem winglet): a exclusao e inofensiva ali (nao ha singularidade de canto
sem winglet, mas a ultima faixa continua estreita e perto da ponta).

Faixa da raiz: pelo mesmo motivo geometrico (espacamento cosseno deixa a
primeira faixa, colada no plano de simetria, muito estreita), a secante de
duas rodadas usada no ajuste linear (cl_norm x alfa) fica numericamente mal
condicionada ali: a inclinacao "b" quase zera e alfa_i = (clmax - a)/b sai
tanto de valores absurdos (milhares de graus, ou negativo) quanto de
valores pequenos e positivos que parecem plausiveis mas nao sao (ver malha
fina do teste de tests/test_analises.py, ou `--lc 4.8` com asa_ns=120: a
faixa da raiz sozinha da alfa_i ~= 4°, bem abaixo do estol real). Uma
checagem so de sinal/finitude em alfa_i nao pega esse caso. Por isso a
mesma logica geometrica do winglet vale no outro extremo da envergadura:
a(s) faixa(s) a ate `DIST_JUNCAO` do plano de simetria (y=0) tambem sao
excluidas da busca da secao critica.
'''

import numpy as np

from avl_run import caso

ETA_CLMAX = [0.1011, 0.398, 0.90]
CLMAX = [1.774, 1.7985, 1.7338]
MACH_BAIXO = 0.2
ALFAS_BASE = (8.0, 14.0)


def clmax_local(eta):
    return np.interp(eta, ETA_CLMAX, CLMAX)


def _excluidas(y, av):
    '''
    Mascara das faixas dentro de DIST_JUNCAO do winglet OU do plano de
    simetria (raiz) (True = excluida). Ver docstring do modulo.
    '''
    semi = av['asa']['yt']
    dist_juncao = av['asa']['ct']/4
    return ((semi - y) <= dist_juncao) | (y <= dist_juncao)


def modelo_estol(arquivo, it, trim, av, empenagem=True):
    '''alfa de estol de cada faixa da asa direita (modelo linear).'''
    semi = av['asa']['yt']
    a0, a1 = ALFAS_BASE
    r0, r1 = (caso(arquivo, MACH_BAIXO, alfa=a, it=it, trim=trim,
                   empenagem=empenagem, faixas=True) for a in ALFAS_BASE)
    eta = r0['faixas']['y']/semi
    c0, c1 = r0['faixas']['cl_norm'], r1['faixas']['cl_norm']
    b = (c1 - c0)/(a1 - a0)
    a = c0 - b*a0
    excluidas = _excluidas(r0['faixas']['y'], av)
    return {'eta': eta, 'alfa_i': (clmax_local(eta) - a)/b, 'excluidas': excluidas}


def estol(arquivo, it, trim, av, empenagem=True):
    m = modelo_estol(arquivo, it, trim, av, empenagem)
    candidatas = np.where(~m['excluidas'])[0]
    i = int(candidatas[np.argmin(m['alfa_i'][candidatas])])
    alfa = float(m['alfa_i'][i])
    r = caso(arquivo, MACH_BAIXO, alfa=alfa, it=it, trim=trim,
             empenagem=empenagem, faixas=True)
    semi = av['asa']['yt']
    eta = r['faixas']['y']/semi
    lim = clmax_local(eta)
    excluidas = _excluidas(r['faixas']['y'], av)
    excesso = r['faixas']['cl_norm'] - lim
    return {'alfa': alfa, 'CL': r['CL'], 'de': r['de'] if trim else 0.0,
            'eta_crit': float(m['eta'][i]), 'eta': eta,
            'cl_norm': r['faixas']['cl_norm'], 'clmax': lim,
            'excluidas': eta[excluidas],
            'excesso_max': float(np.max(excesso[~excluidas])),
            'alfa_i': m['alfa_i']}
