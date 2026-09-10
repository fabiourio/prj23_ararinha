'''
Testes das partes BARATAS do otimizador -- as que nao chamam o eulerblock.

Rodar antes de disparar a otimizacao: um erro no gradiente da restricao ou no
ponto de partida so apareceria depois de horas de maquina.

Rodar com:  python test_otimiza_secao.py
'''

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import descritores as dsc
import otimiza_secao as ot


def test_t01_bate_com_o_descritor():
    '''
    O t01 analitico do otimizador tem de dar EXATAMENTE o mesmo valor que o
    descritor usado no DOE. Se divergirem, o limiar foi calibrado numa funcao
    e aplicado noutra.
    '''
    for Au, Al in [(ot.AU_1411, ot.AL_1411),
                   ot.ponto_de_partida(0.1772)[::-1],
                   ot.ponto_de_partida(0.1082)[::-1]]:
        d = dsc.descritores(Au, Al)
        t01_ot = ot.t01_de(Al, Au)
        assert abs(t01_ot - d['t_01']) < 1e-12, \
            f'otimizador {t01_ot:.12f} != descritor {d["t_01"]:.12f}'


def test_gradiente_de_bluntez_bate_com_diferencas_finitas():
    '''
    A bluntez t01/sqrt(t_max) mistura um termo linear (t01) com um termo que
    vem da funcao KS (t_max). Conferir a regra do quociente contra diferencas
    finitas, usando o t_max do descritor como referencia.
    '''
    Al, Au = ot.ponto_de_partida(0.1772)
    h = 1e-6

    def bl_de(al, au):
        return ot.bluntez(ot.t01_de(al, au), dsc.descritores(au, al)['t_max'])

    # gradiente de t_max por diferencas finitas, para alimentar grad_bluntez
    xx = np.hstack([Al, Au, 0.05])
    dtmax = np.zeros_like(xx)
    for i in range(len(xx) - 1):
        xp, xm = xx.copy(), xx.copy()
        xp[i] += h
        xm[i] -= h
        alp, aup, _ = ot.desmonta(xp)
        alm, aum, _ = ot.desmonta(xm)
        dtmax[i] = (dsc.descritores(aup, alp)['t_max']
                    - dsc.descritores(aum, alm)['t_max']) / (2 * h)

    t01 = ot.t01_de(Al, Au)
    tmax = dsc.descritores(Au, Al)['t_max']
    g_analitico = ot.grad_bluntez(t01, tmax, dtmax)

    g_fd = np.zeros_like(xx)
    for i in range(len(xx)):
        xp, xm = xx.copy(), xx.copy()
        xp[i] += h
        xm[i] -= h
        alp, aup, _ = ot.desmonta(xp)
        alm, aum, _ = ot.desmonta(xm)
        g_fd[i] = (bl_de(alp, aup) - bl_de(alm, aum)) / (2 * h)

    erro = np.max(np.abs(g_analitico - g_fd)) / max(np.abs(g_fd).max(), 1e-12)
    assert erro < 1e-4, \
        f'gradiente da bluntez diverge: erro relativo {erro:.2e}\n' \
        f'  analitico: {g_analitico}\n  dif.finita: {g_fd}'


def test_gradiente_de_t01_bate_com_diferencas_finitas():
    '''
    O gradiente de t01 e analitico (a superficie CST e linear nos
    coeficientes para x fixo). Conferir contra diferenca central.
    '''
    xx = np.hstack([ot.AL_1411, ot.AU_1411, 0.05])
    g_analitico = ot.grad_t01()

    h = 1e-7
    g_fd = np.zeros_like(xx)
    for i in range(len(xx)):
        xp, xm = xx.copy(), xx.copy()
        xp[i] += h
        xm[i] -= h
        alp, aup, _ = ot.desmonta(xp)
        alm, aum, _ = ot.desmonta(xm)
        g_fd[i] = (ot.t01_de(alp, aup) - ot.t01_de(alm, aum)) / (2 * h)

    erro = np.max(np.abs(g_analitico - g_fd))
    assert erro < 1e-6, \
        f'gradiente de t01 diverge: erro maximo {erro:.2e}\n' \
        f'  analitico: {g_analitico}\n  dif.finita: {g_fd}'


def test_ponto_de_partida_atende_as_restricoes():
    '''
    O NACA 1411 puro tem t/c = 0,115 e violaria a restricao de espessura nas
    tres estacoes. O ponto de partida reescalado precisa atender espessura,
    bluntez e nao-cruzamento, senao o SLSQP comeca inviavel em varias
    restricoes ao mesmo tempo.
    '''
    for nome, est in ot.ESTACOES.items():
        Al, Au = ot.ponto_de_partida(est['tc_ref'])
        d = dsc.descritores(Au, Al)
        assert abs(d['t_max'] - est['tc_ref']) < 0.01, \
            f'{nome}: t/c de partida {d["t_max"]:.4f} != alvo {est["tc_ref"]:.4f}'
        bl = ot.bluntez(ot.t01_de(Al, Au), d['t_max'])
        assert bl >= ot.BLUNTEZ_MIN, \
            f'{nome}: bluntez de partida {bl:.5f} abaixo do limiar ' \
            f'{ot.BLUNTEZ_MIN} (t/c = {d["t_max"]:.4f})'
        assert d['t_min'] >= ot.MINT_MIN, \
            f'{nome}: superficies cruzadas na partida (t_min = {d["t_min"]:.5f})'


def test_ponto_de_partida_preserva_arqueamento():
    '''
    O reescalonamento mexe so na espessura; a linha media do NACA 1411 tem de
    sobreviver, senao mudamos o perfil de partida sem querer.
    '''
    d0 = dsc.descritores(ot.AU_1411, ot.AL_1411)
    Al, Au = ot.ponto_de_partida(0.1772)
    d1 = dsc.descritores(Au, Al)
    assert abs(d1['c_max'] - d0['c_max']) < 1e-6, \
        f'arqueamento mudou: {d0["c_max"]:.6f} -> {d1["c_max"]:.6f}'


def test_estacoes_dentro_dos_batentes():
    '''As tres estacoes tem de ser alcancaveis dentro dos batentes de CST.'''
    for nome, est in ot.ESTACOES.items():
        Al, Au = ot.ponto_de_partida(est['tc_ref'])
        assert np.all(Al >= np.array(ot.AL_LOWER)) and \
               np.all(Al <= np.array(ot.AL_UPPER)), \
            f'{nome}: Al de partida fora dos batentes: {Al}'
        assert np.all(Au >= np.array(ot.AU_LOWER)) and \
               np.all(Au <= np.array(ot.AU_UPPER)), \
            f'{nome}: Au de partida fora dos batentes: {Au}'


if __name__ == '__main__':
    falhas = 0
    for nome, fn in sorted(globals().items()):
        if nome.startswith('test_') and callable(fn):
            try:
                fn()
                print(f'PASS  {nome}')
            except AssertionError as exc:
                falhas += 1
                print(f'FALHA {nome}\n      {exc}')
            except Exception as exc:                          # noqa: BLE001
                falhas += 1
                print(f'ERRO  {nome}\n      {type(exc).__name__}: {exc}')
    print(f'\n{falhas} falha(s)')
    sys.exit(1 if falhas else 0)
