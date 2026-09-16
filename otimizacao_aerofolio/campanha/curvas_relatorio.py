'''
Roda os solvers UMA VEZ e despeja as curvas do relatorio em CSV.

Separado de figuras_relatorio.py de proposito: desenhar bem exige iterar no
estilo, e antes disso cada tentativa custava ~6 min de solver porque as curvas
eram recomputadas a cada plot. Com as curvas em disco, replotar e instantaneo
-- e elas passam a ser versionadas junto com os numeros, em vez de viverem
so no historico.pickle, que o .gitignore descarta.

Escreve em relatorio/tables/:

  curva_geometria.csv      perfil, x, y            (3 perfis)
  curva_cp_mach.csv        perfil, x, cp, mach     (Euler, condicao de projeto)
  curva_subsonico.csv      perfil, alpha, cl, cd, cm  (XFoil viscoso)
  curva_convergencia.csv   avaliacao, cd, cl, tc_max, bluntez, tempo_min

A polar transonica ja vem pronta de ENTREGA/07_polar_transonica e vira
tab4_polar_transonica.csv em tabelas_do_relatorio.py.

Rodar com:  python curvas_relatorio.py [--geometria] [--cpmach]
                                        [--subsonico] [--convergencia]
            (sem argumento faz tudo; ~6 min)
'''

import csv
import os
import pickle
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
PACOTE = os.path.normpath(os.path.join(AQUI, '..', 'eulerblock', 'package'))
ENTREGA = os.path.normpath(os.path.join(AQUI, '..', 'ENTREGA'))
TABLES = os.path.normpath(os.path.join(AQUI, '..', '..', 'relatorio', 'tables'))
RES = os.path.join(AQUI, 'resultados')
sys.path.insert(0, PACOTE)
sys.path.insert(0, AQUI)

TC_REF, CL_REF, MACH_N = 0.1772, 0.7319, 0.7196
RE_DECOLAGEM = 4.23e7
ALPHA_PARTIDA_PADRAO = 3.567408

ROT_PART, ROT_OTI, ROT_CUSP = 'partida', 'otimizado', 'cusp_liberado'
ROT_SEM = 'sem_restricao_clmax'

# A variante sem restricao de cl_max e o CONTROLE do item 9: mostra o que a
# substituta geometrica esta comprando. Nao tem .dat exportado -- nao e um
# perfil que a equipe entrega -- entao os coeficientes vem da linha
# "otimizado sem restr. de cl_max" da tabela da MAC em
# ENTREGA/01_ponto_de_projeto/tabelas_1_e_2.md.
SEM_CLMAX_AL = np.array([-0.0724465, -0.445839, -0.176163, -0.05])
SEM_CLMAX_AU = np.array([0.190881, 0.0759592, 0.413892, 0.138493])


def _escreve(arq, cabecalho, linhas):
    os.makedirs(TABLES, exist_ok=True)
    cam = os.path.join(TABLES, arq)
    with open(cam, 'w', newline='', encoding='utf-8') as fid:
        w = csv.writer(fid)
        w.writerow(cabecalho)
        w.writerows(linhas)
    print(f'  {arq}: {len(linhas)} linhas')


def le_cst(nome):
    '''Le Al, Au e alpha de um *_cst.txt da ENTREGA (fonte canonica).'''
    import re
    txt = open(os.path.join(ENTREGA, '02_perfis_otimizados',
                            f'{nome}_cst.txt')).read()

    def vetor(chave):
        m = re.search(chave + r'\s*=\s*\[([^\]]*)\]', txt)
        return np.array([float(v) for v in m.group(1).split(',')])

    return (vetor('Al'), vetor('Au'),
            float(re.search(r'alpha\s*=\s*([-\d.eE+]+)', txt).group(1)))


def alpha_partida():
    import json
    cam = os.path.join(RES, 'ativ2_alpha_naca.json')
    if os.path.isfile(cam):
        with open(cam) as fid:
            return float(json.load(fid)['alpha_final_deg'])
    print(f'  aviso: ativ2_alpha_naca.json ausente; usando '
          f'{ALPHA_PARTIDA_PADRAO:.6f} graus')
    return ALPHA_PARTIDA_PADRAO


def _le_dat(nome):
    d = np.loadtxt(os.path.join(ENTREGA, '02_perfis_otimizados', f'{nome}.dat'),
                   skiprows=1)
    return d[:, 0], d[:, 1]


# ---------------------------------------------------------------------------
def dump_geometria():
    from otimiza_secao import ponto_de_partida
    from airfoil_mod import cstfoil

    Al_i, Au_i = ponto_de_partida(TC_REF)
    d0 = cstfoil(Au_i, Al_i, 0.5 * (1 - np.cos(np.linspace(0, np.pi, 200))))
    perfis = [(ROT_PART, np.real(d0['x_coord']), np.real(d0['y_coord'])),
              (ROT_OTI, *_le_dat('meio_roteiro')),
              (ROT_CUSP, *_le_dat('meio_cusp_liberado'))]

    linhas = [[rot, f'{xi:.8f}', f'{yi:.8f}']
              for rot, x, y in perfis for xi, yi in zip(x, y)]
    _escreve('curva_geometria.csv', ['perfil', 'x', 'y'], linhas)


def dump_cp_mach():
    from otimiza_secao import ponto_de_partida
    from analisa_otimos import roda_euler

    Al_i, Au_i = ponto_de_partida(TC_REF)
    Al_o, Au_o, a_oti = le_cst('meio_roteiro')
    a_ini = alpha_partida()

    linhas, escalares = [], []
    for rot, Al, Au, a in ((ROT_PART, Al_i, Au_i, a_ini),
                           (ROT_OTI, Al_o, Au_o, a_oti)):
        print(f'  Euler: {rot} (alpha = {a:.4f} deg) ...', flush=True)
        r = roda_euler(Al, Au, a * np.pi / 180, 1.0)
        print(f"    cl={r['CL']:.5f} cd={r['CD']:.6f} cm={r['CM']:.5f} "
              f"Mach_max={np.max(r['Mach']):.4f}")
        linhas += [[rot, f'{x:.8f}', f'{cp:.8f}', f'{m:.8f}']
                   for x, cp, m in zip(r['x'], r['Cp'], r['Mach'])]
        escalares.append([rot, f'{a:.6f}', f'{r["CL"]:.6f}', f'{r["CD"]:.8f}',
                          f'{r["CM"]:.6f}', f'{np.max(r["Mach"]):.6f}'])
    _escreve('curva_cp_mach.csv', ['perfil', 'x', 'cp', 'mach'], linhas)
    # os escalares do Euler nao sao derivaveis das curvas, entao ficam num
    # arquivo proprio -- e o que tabelas_do_relatorio.py le para a Tab. 2
    _escreve('curva_escalares.csv',
             ['perfil', 'alpha_deg', 'cl', 'cd', 'cm', 'mach_max'], escalares)


def dump_subsonico():
    from otimiza_secao import ponto_de_partida
    import xfoil_runner as xr

    Al_i, Au_i = ponto_de_partida(TC_REF)
    Al_o, Au_o, _ = le_cst('meio_roteiro')
    Al_c, Au_c, _ = le_cst('meio_cusp_liberado')

    linhas = []
    for rot, Au, Al in ((ROT_PART, Au_i, Al_i), (ROT_OTI, Au_o, Al_o),
                        (ROT_CUSP, Au_c, Al_c),
                        (ROT_SEM, SEM_CLMAX_AU, SEM_CLMAX_AL)):
        print(f'  XFoil: {rot} em Re = {RE_DECOLAGEM:.2e} ...', flush=True)
        # ate 28 graus: com 20 a varredura parava antes do estol
        r = xr.clmax_cst(Au, Al, Re=RE_DECOLAGEM, Mach=0.0,
                         alpha_seq=(-6.0, 28.0, 0.5), timeout=900)
        p = r['polar']
        print(f"    {r['n_pontos']} pontos, cl_max = {r['clmax']:.4f} "
              f"({r['situacao']})")
        linhas += [[rot, f'{a:.4f}', f'{cl:.6f}', f'{cd:.6f}', f'{cm:.6f}']
                   for a, cl, cd, cm in zip(p['alpha'], p['CL'],
                                            p['CD'], p['CM'])]
    _escreve('curva_subsonico.csv',
             ['perfil', 'alpha_deg', 'cl', 'cd', 'cm'], linhas)


def dump_convergencia():
    '''
    Historico da otimizacao. O pickle e regenerado a cada rodada e fica fora
    do git; passa-lo para CSV e o que torna a figura reproduzivel sem gastar
    os 17,6 min da otimizacao de novo.
    '''
    import otimiza_secao as ot

    cam = os.path.join(RES, 'otim_meio', 'historico.pickle')
    if not os.path.isfile(cam):
        print('  AVISO: otim_meio/historico.pickle ausente -- '
              'rode `python otimiza_secao.py meio` antes')
        return
    with open(cam, 'rb') as fid:
        h = pickle.load(fid)['hist']

    limiar = ot.limiar_bluntez('meio')
    linhas = []
    for i, (cd, cl, mt, t01, t) in enumerate(zip(h['CD'], h['CL'], h['maxt'],
                                                 h['t01'], h['tempo']), 1):
        linhas.append([i, f'{cd:.8f}', f'{cl:.6f}', f'{mt:.6f}',
                       f'{ot.bluntez(t01, mt):.6f}', f'{t / 60:.3f}'])
    _escreve('curva_convergencia.csv',
             ['avaliacao', 'cd', 'cl', 'tc_max', 'bluntez', 'tempo_min'],
             linhas)
    print(f'    alvos: cl = {CL_REF}, tc_min = {TC_REF}, '
          f'limiar_bluntez = {limiar:.4f}')


def main():
    flags = [a for a in sys.argv[1:] if a.startswith('--')]
    tudo = not flags
    print(f'despejando curvas em {TABLES}')
    if tudo or '--geometria' in flags:
        print('=== geometria ===')
        dump_geometria()
    if tudo or '--cpmach' in flags:
        print('=== Cp e Mach (Euler) ===')
        dump_cp_mach()
    if tudo or '--subsonico' in flags:
        print('=== subsonico (XFoil) ===')
        dump_subsonico()
    if tudo or '--convergencia' in flags:
        print('=== convergencia ===')
        dump_convergencia()
    return 0


if __name__ == '__main__':
    sys.exit(main())
