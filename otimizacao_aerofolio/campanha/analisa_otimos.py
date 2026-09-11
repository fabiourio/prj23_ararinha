'''
Pos-processamento das otimizacoes de secao -- Lab 03, PRJ-23.

Responde tres perguntas, nesta ordem de importancia:

  1. Quanto o otimizador ganhou, e a restricao de sustentacao maxima ficou
     ATIVA? Comparando a rodada com e sem ela, o custo em arrasto da
     restricao e a diferenca entre os dois otimos.

  2. O perfil otimo entrega mesmo cl_max >= 1,80? A substituta e um modelo
     com dispersao (R2 ~ 0,51); so o XFoil viscoso decide.

  3. A melhoria SOBREVIVE ao refinamento de malha? O nivel 1,0 superestima o
     arrasto em 70% (verificacao_malha.py), e a dissipacao numerica depende da
     intensidade do choque -- que e exatamente o que estamos minimizando.
     Se a melhoria sumir na malha fina, o resultado e artefato.

Rodar com:  python analisa_otimos.py [--malha] [--xfoil]
    --malha  reavalia partida e otimo nos niveis 1,5 e 2,0 (caro, ~30 min)
    --xfoil  verifica cl_max de cada otimo (rapido)
Sem opcoes, faz as duas.
'''

import os
import pickle
import shutil
import sys
import tempfile

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
PACOTE = os.path.normpath(os.path.join(AQUI, '..', 'eulerblock', 'package'))
sys.path.insert(0, PACOTE)
sys.path.insert(0, AQUI)

import descritores as dsc
import otimiza_secao as ot
import xfoil_runner as xr
from doe_clmax_xfoil import (ALPHA_SEQ, CLMAX_ALVO, MACH_DECOLAGEM,
                             RE_DECOLAGEM)

RES = os.path.join(AQUI, 'resultados')
NIVEIS_VERIFICACAO = [1.5, 2.0]


def carrega(pasta):
    cam = os.path.join(RES, pasta, 'historico.pickle')
    if not os.path.isfile(cam):
        return None
    with open(cam, 'rb') as fid:
        d = pickle.load(fid)
    h = d['hist']
    if not h['xx']:
        return None
    d['xx_ini'] = np.asarray(h['xx'][0])
    # o otimo e o ultimo ponto viavel de menor CD
    cd = np.array(h['CD'])
    cl = np.array(h['CL'])
    maxt = np.array(h['maxt'])
    viavel = (np.abs(cl - d['cl_ref']) < 5e-3) & (maxt >= d['tc_ref'] - 1e-4)
    idx = int(np.argmin(np.where(viavel, cd, np.inf))) if viavel.any() \
        else int(np.argmin(cd))
    d['i_otimo'] = idx
    d['xx_otimo'] = np.asarray(h['xx'][idx])
    d['n_aval'] = len(cd)
    d['tempo_min'] = h['tempo'][-1] / 60
    return d


def roda_euler(Al, Au, alpha, nivel, cfl=None):
    '''
    Uma avaliacao do Euler num nivel de malha, em diretorio isolado.

    O cfl precisa ser o da ESTACAO: com 0,20 o perfil grosso de partida da
    raiz diverge e a tabela sai com nan. Se nao for informado, tenta 0,20 e
    recua, igual ao avaliador da otimizacao.
    '''
    from eulerblock import euler_mod as eb
    nchord = int(30 * nivel) + 1
    nj = int(48 * nivel) + 1
    s0 = 48 / (nj - 1) * 0.5e-2
    cfls = [cfl] if cfl else [0.20, 0.10, 0.05]
    tmp = tempfile.mkdtemp(prefix='pos_')
    cwd = os.getcwd()
    try:
        os.chdir(tmp)
        for c in cfls:
            r = eb.run_cst(Al, Au, nchord, alpha, ot.MACH_N,
                           gamma=1.4, order=2, iter=20000, dt=0.001, CFL=c,
                           use_local_dt=1, res_NK=1e-5, res_tol=1e-8,
                           reinitialize=0, adj_funcs=[], plot=False,
                           NJ=nj, s0=s0)
            if np.isfinite(r['CL']) and np.isfinite(r['CD']):
                break
        return {'CL': float(r['CL']), 'CD': float(r['CD']),
                'CM': float(r['CM']), 'nivel': nivel,
                'Cp': np.asarray(r['distrib']['Cp']),
                'x': np.asarray(r['distrib']['xx']),
                'Mach': np.asarray(r['distrib']['Mach'])}
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)


def resumo_estacao(nome, com, sem):
    '''Tabela comparativa de uma estacao.'''
    print(f'\n{"="*70}')
    print(f'ESTACAO {nome.upper()}')
    print('=' * 70)
    if com is None:
        print('  (sem resultado)')
        return

    h = com['hist']
    cd0, cdf = h['CD'][0], h['CD'][com['i_otimo']]
    print(f'  avaliacoes: {com["n_aval"]}   tempo: {com["tempo_min"]:.0f} min')
    print(f'  CD: {cd0:.6f} -> {cdf:.6f}   ({(cdf/cd0-1)*100:+.1f}%)')
    print(f'  CL: {h["CL"][com["i_otimo"]]:.5f}  (alvo {com["cl_ref"]:.4f})')
    print(f'  t/c: {h["maxt"][com["i_otimo"]]:.5f}  (min {com["tc_ref"]:.4f})')

    Al, Au, alpha = ot.desmonta(com['xx_otimo'])
    d = dsc.descritores(Au, Al)
    bl = ot.bluntez(ot.t01_de(Al, Au), h['maxt'][com['i_otimo']])
    lim = ot.limiar_bluntez(nome)
    ativa = abs(bl - lim) < 2e-4
    print(f'  bluntez: {bl:.5f}  (limiar {lim:.4f})  '
          f'{"ATIVA" if ativa else "folgada"}')
    print(f'  alpha: {alpha*180/np.pi:.4f} deg')
    print(f'  espessura maxima em x/c = {d["x_tmax"]:.4f}  '
          f'(partida: NACA de 4 digitos, 0,30)')
    print(f'  arqueamento maximo {d["c_max"]:.5f} em x/c = {d["x_cmax"]:.4f}')

    if sem is not None:
        hs = sem['hist']
        cdf_s = hs['CD'][sem['i_otimo']]
        Als, Aus, _ = ot.desmonta(sem['xx_otimo'])
        bl_s = ot.bluntez(ot.t01_de(Als, Aus), hs['maxt'][sem['i_otimo']])
        print(f'\n  --- comparacao com a rodada SEM a restricao ---')
        print(f'  CD sem restricao: {cdf_s:.6f}   com: {cdf:.6f}')
        custo = (cdf / cdf_s - 1) * 100 if cdf_s else np.nan
        print(f'  custo da restricao: {custo:+.2f}% de arrasto')
        print(f'  bluntez do otimo irrestrito: {bl_s:.5f}  '
              f'({"atenderia o limiar sozinho" if bl_s >= lim else "VIOLA o limiar"})')
        if bl_s >= lim:
            print('  -> a restricao era INATIVA: nao custou nada.')
        else:
            print('  -> a restricao MORDEU: o custo acima e o preco de')
            print('     manter clmax_w = 1,80 no designTool.')


def verifica_xfoil(nome, dados):
    Al, Au, _ = ot.desmonta(dados['xx_otimo'])
    r = xr.clmax_cst(Au, Al, Re=RE_DECOLAGEM, Mach=MACH_DECOLAGEM,
                     alpha_seq=ALPHA_SEQ, timeout=900)
    ok = r['clmax'] >= CLMAX_ALVO
    print(f'  {nome:6s}: cl_max = {r["clmax"]:.4f} em alpha = '
          f'{r["alpha_clmax"]:.1f} deg  ({r["situacao"]})  '
          f'{"ATENDE" if ok else "ABAIXO do alvo " + str(CLMAX_ALVO)}')
    return r


def verifica_malha(nome, dados):
    print(f'\n  --- {nome}: a melhoria sobrevive ao refinamento? ---')
    Ali, Aui, ai = ot.desmonta(dados['xx_ini'])
    Alo, Auo, ao = ot.desmonta(dados['xx_otimo'])
    print(f'{"nivel":>6s} {"CD partida":>11s} {"CD otimo":>10s} {"ganho":>9s}')
    for nivel in [1.0] + NIVEIS_VERIFICACAO:
        ri = roda_euler(Ali, Aui, ai, nivel)
        ro = roda_euler(Alo, Auo, ao, nivel)
        ganho = (ro['CD'] / ri['CD'] - 1) * 100
        print(f'{nivel:6.2f} {ri["CD"]:11.6f} {ro["CD"]:10.6f} {ganho:8.1f}%',
              flush=True)


def main():
    fazer_malha = '--malha' in sys.argv or len(sys.argv) == 1
    fazer_xfoil = '--xfoil' in sys.argv or len(sys.argv) == 1

    resultados = {}
    for nome in ot.ESTACOES:
        com = carrega(f'otim_{nome}')
        sem = carrega(f'otim_{nome}_sem_bluntez')
        resultados[nome] = (com, sem)
        resumo_estacao(nome, com, sem)

    if fazer_xfoil:
        print(f'\n{"="*70}\nVERIFICACAO DE cl_max NO XFOIL '
              f'(Re = {RE_DECOLAGEM:.1e}, M = {MACH_DECOLAGEM})\n{"="*70}')
        for nome, (com, sem) in resultados.items():
            if com is not None:
                verifica_xfoil(nome, com)
            if sem is not None:
                verifica_xfoil(nome + '/sem', sem)

    if fazer_malha:
        print(f'\n{"="*70}\nVERIFICACAO DE MALHA DOS OTIMOS\n{"="*70}')
        for nome, (com, _) in resultados.items():
            if com is not None:
                verifica_malha(nome, com)


if __name__ == '__main__':
    main()
