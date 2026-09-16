'''
Dados numericos das tabelas do relatorio -- estacao MAC.

Escreve CSVs em relatorio/tables/, um por tabela do relatorio. O texto do
relatorio vive no Overleaf; este repositorio guarda as figuras e os numeros
que as alimentam, de forma que qualquer valor citado no relatorio possa ser
rastreado ate a rodada que o produziu.

Fontes de cada tabela:

  tab1_secao_referencia.csv   ENTREGA/01_ponto_de_projeto/tabelas_1_e_2.md
  tab2_alpha_naca1411.csv     resultados/ativ2_alpha_naca.json
  tab3_aerofolios.csv         ENTREGA/.../*_cst.txt + resultados/item5_*.json
  tab4_polar_ponto_projeto.csv  ENTREGA/07_polar_transonica/polar_*.csv
  tab5_subsonico.csv          resultados/item8_subsonico.json
  tab6_variantes_item9.csv    resultados/item9_variantes.json

Rodar com:  python tabelas_do_relatorio.py
'''

import csv
import json
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
ENTREGA = os.path.normpath(os.path.join(AQUI, '..', 'ENTREGA'))
RES = os.path.join(AQUI, 'resultados')
TABLES = os.path.normpath(os.path.join(AQUI, '..', '..', 'relatorio', 'tables'))

sys.path.insert(0, AQUI)


def _json(nome):
    cam = os.path.join(RES, nome)
    if not os.path.isfile(cam):
        print(f'  AVISO: {nome} nao existe -- rode o script que o gera')
        return None
    with open(cam) as fid:
        return json.load(fid)


def escreve(arq, cabecalho, linhas):
    os.makedirs(TABLES, exist_ok=True)
    cam = os.path.join(TABLES, arq)
    with open(cam, 'w', newline='', encoding='utf-8') as fid:
        w = csv.writer(fid)
        w.writerow(cabecalho)
        w.writerows(linhas)
    print(f'  {arq}: {len(linhas)} linhas')


# ---------------------------------------------------------------------------
def tab1():
    '''Ponto de projeto e secao de referencia (item 1).'''
    escreve('tab1_secao_referencia.csv',
            ['parametro', 'valor', 'unidade', 'descricao'], [
                ['h', 35000, 'ft', 'altitude do ponto de projeto'],
                ['M_inf', 0.850, '-', 'Mach do ponto de projeto'],
                ['W', 229669, 'kgf', 'peso da aeronave no ponto de projeto'],
                ['S_ref', 368.833, 'm2', 'area de referencia da aeronave'],
                ['rho_inf', 0.38046, 'kg/m3', 'densidade (ISA)'],
                ['a_inf', 296.587, 'm/s', 'velocidade do som (ISA)'],
                ['mu_inf', 1.4451e-05, 'Pa.s', 'viscosidade (Sutherland)'],
                ['V_inf', 252.1, 'm/s', 'velocidade de voo'],
                ['q_inf', 12089.7, 'Pa', 'pressao dinamica'],
                ['CL', 0.5053, '-', 'coeficiente de sustentacao da aeronave'],
                ['sweep_c2', 32.158, 'graus', 'enflechamento a 50% da corda'],
                ['M_n', 0.7196, '-', 'Mach normal ao enflechamento'],
                ['eta', 0.398, '-', 'posicao da estacao MAC'],
                ['c_ref', 6.899, 'm', 'corda da secao de referencia (MAC)'],
                ['c_n', 5.841, 'm', 'corda no plano normal'],
                ['cl_ref', 0.7319, '-', 'cl da secao no plano normal'],
                ['Re_ref', 3.282e+07, '-', 'Reynolds da secao de referencia'],
                ['tc_ref', 0.1772, '-', 'espessura relativa base'],
                ['Re_decolagem', 4.23e+07, '-', 'Reynolds na decolagem (item 8)'],
            ])


def tab2():
    '''Ajuste do alpha do NACA 1411 ao cl de projeto (item 2).'''
    d = _json('ativ2_alpha_naca.json')
    if d is None:
        return
    linhas = [[i + 1, round(h['alpha_deg'], 6), round(h['cl'], 6),
               round(h['cd'], 6)] for i, h in enumerate(d['hist'])]
    escreve('tab2_alpha_naca1411.csv',
            ['avaliacao', 'alpha_deg', 'cl', 'cd'], linhas)


def descritores(x, y):
    '''Espessura e arqueamento maximos por interpolacao nas duas superficies.

    Calculado aqui, e nao lido do cstfoil, porque o cstfoil devolve
    x_max_thickness errado -- filtra o vetor de espessuras e indexa o de
    abscissas completo (ver RESULTADOS.md, defeito 2).
    '''
    i = int(np.argmin(x))
    xl, yl = x[:i + 1][::-1], y[:i + 1][::-1]
    xu, yu = x[i:], y[i:]
    g = np.linspace(0.0, 1.0, 801)
    t = np.interp(g, xu, yu) - np.interp(g, xl, yl)
    c = 0.5 * (np.interp(g, xu, yu) + np.interp(g, xl, yl))
    it, ic = int(np.argmax(t)), int(np.argmax(np.abs(c)))
    return dict(tc=t[it], x_tc=g[it], camber=c[ic], x_camber=g[ic])


def _curva(arq):
    '''CSV longo -> {perfil: {coluna: array}}.'''
    from collections import defaultdict
    cam = os.path.join(TABLES, arq)
    if not os.path.isfile(cam):
        print(f'  AVISO: {arq} nao existe -- rode curvas_relatorio.py')
        return None
    bruto = defaultdict(lambda: defaultdict(list))
    with open(cam, encoding='utf-8') as fid:
        r = csv.reader(fid)
        cab = next(r)
        for lin in r:
            for col, v in zip(cab[1:], lin[1:]):
                bruto[lin[0]][col].append(float(v))
    return {k: {c: np.array(v) for c, v in d.items()}
            for k, d in bruto.items()}


def tab3():
    '''Coeficientes CST e descritores dos dois perfis (itens 3 e 5).'''
    from otimiza_secao import ponto_de_partida
    from curvas_relatorio import le_cst, TC_REF

    geo = _curva('curva_geometria.csv')
    esc = _curva('curva_escalares.csv')
    if geo is None or esc is None:
        return

    Al_i, Au_i = ponto_de_partida(TC_REF)
    Al_o, Au_o, _ = le_cst('meio_roteiro')
    di = descritores(geo['partida']['x'], geo['partida']['y'])
    do = descritores(geo['otimizado']['x'], geo['otimizado']['y'])

    def col(Al, Au, chave, desc):
        a = esc[chave]
        v = [f'{x:.8f}' for x in list(Al) + list(Au)]
        v.append(f'{a["alpha_deg"][0]:.6f}')
        v += [f'{a[k][0]:.6f}' for k in ('cl', 'cd', 'cm')]
        v += [f'{desc[k]:.6f}'
              for k in ('tc', 'x_tc', 'camber', 'x_camber')]
        return v

    nomes = (['Al1', 'Al2', 'Al3', 'Al4', 'Au1', 'Au2', 'Au3', 'Au4', 'alpha_deg',
              'cl', 'cd', 'cm', 'tc_max', 'x_tc_max', 'camber_max',
              'x_camber_max'])
    ci = col(Al_i, Au_i, 'partida', di)
    co = col(Al_o, Au_o, 'otimizado', do)
    escreve('tab3_aerofolios.csv',
            ['parametro', 'naca1411_reescalado', 'otimizado'],
            [[n, a, b] for n, a, b in zip(nomes, ci, co)])


def tab4():
    '''Polar transonica na vizinhanca do ponto de projeto (item 7).'''
    cam = os.path.join(ENTREGA, '07_polar_transonica',
                       'polar_transonica_meio.csv')
    if not os.path.isfile(cam):
        print('  AVISO: polar_transonica_meio.csv nao encontrada')
        return
    with open(cam) as fid:
        linhas = [r for r in csv.DictReader(fid)]
    escreve('tab4_polar_transonica.csv',
            ['perfil', 'alpha_deg', 'cl', 'cd', 'cm'],
            [[r['perfil'], r['alpha_deg'], r['CL'], r['CD'], r['CM']]
             for r in linhas])


VIES_XFOIL = 0.24        # o XFoil le cl_max 0,24 acima do experimental
                         # (afere_xfoil.py, 5 perfis NACA de Abbott)

# c_d transonico de cada variante, medido no Euler (ENTREGA/RESULTADOS.md).
# So o do perfil do roteiro sai de curva_escalares.csv; as variantes do item 9
# nao sao re-rodadas no Euler por este pipeline.
CD_TRANSONICO = {'otimizado': None,              # lido de curva_escalares
                 'cusp_liberado': 0.0082152,
                 'sem_restricao_clmax': 0.0106635}


def _metricas(p):
    '''Metricas subsonicas de uma polar do XFoil.'''
    a, cl, cd, cm = (p[k] for k in ('alpha_deg', 'cl', 'cd', 'cm'))
    m = (a >= -2) & (a <= 4)
    i = int(np.argmin(cd))
    dentro = cl[cd <= 1.10 * cd[i]]
    return dict(
        clmax=float(np.max(cl)),
        cl_alpha=float(np.polyfit(a[m], cl[m], 1)[0]),
        cd_min=float(cd[i]),
        alpha_estol=float(a[int(np.argmax(cl))]),
        cm_medio=float(np.mean(cm[m])),
        largura_bacia=(float(dentro.max() - dentro.min())
                       if len(dentro) else float('nan')))


def tab5():
    '''Comparacao subsonica dos dois perfis (item 8).'''
    d = _curva('curva_subsonico.csv')
    if d is None:
        return
    mi, mo = _metricas(d['partida']), _metricas(d['otimizado'])
    chaves = [('cl_alpha', 'cl_alpha [1/grau]'),
              ('clmax', 'cl_max (XFoil)'),
              ('alpha_estol', 'alpha_estol [graus]'),
              ('cd_min', 'cd_min'),
              ('cm_medio', 'cm medio (-2 a 4 graus)')]
    linhas = [[rot, f'{mi[k]:.6f}', f'{mo[k]:.6f}'] for k, rot in chaves]
    linhas.append([f'cl_max (corrigido, vies -{VIES_XFOIL:.2f})',
                   f'{mi["clmax"] - VIES_XFOIL:.6f}',
                   f'{mo["clmax"] - VIES_XFOIL:.6f}'])
    escreve('tab5_subsonico.csv',
            ['grandeza', 'naca1411_reescalado', 'otimizado'], linhas)


def tab6():
    '''Variantes de formulacao: efeito transonico e subsonico (item 9).'''
    d = _curva('curva_subsonico.csv')
    esc = _curva('curva_escalares.csv')
    if d is None:
        return
    rotulo = {'otimizado': 'batentes do roteiro',
              'cusp_liberado': 'cusp liberado',
              'sem_restricao_clmax': 'sem restricao de cl_max'}
    linhas = []
    for chave, rot in rotulo.items():
        if chave not in d:
            continue
        m = _metricas(d[chave])
        cdt = CD_TRANSONICO[chave]
        if cdt is None and esc and chave in esc:
            cdt = float(esc[chave]['cd'][0])
        linhas.append([rot, '' if cdt is None else f'{cdt:.6f}',
                       f'{m["clmax"]:.4f}', f'{m["cd_min"]:.5f}',
                       f'{m["alpha_estol"]:.1f}', f'{m["cm_medio"]:+.4f}',
                       f'{m["largura_bacia"]:.3f}'])
    escreve('tab6_variantes_item9.csv',
            ['formulacao', 'cd_transonico', 'clmax_xfoil', 'cd_min',
             'alpha_estol_deg', 'cm_medio', 'largura_bacia_cl'], linhas)


def main():
    print(f'escrevendo em {TABLES}')
    for f in (tab1, tab2, tab3, tab4, tab5, tab6):
        f()
    return 0


if __name__ == '__main__':
    sys.exit(main())
