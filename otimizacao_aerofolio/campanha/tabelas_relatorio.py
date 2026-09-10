'''
Tabelas 1 e 2 do roteiro, prontas para o relatorio.

  Tab. 1 -- dados da secao de referencia (uma por estacao)
  Tab. 2 -- coeficientes CST e propriedades do perfil de partida e do otimizado

Sai em Markdown (para conferir aqui) e em LaTeX (para colar no Overleaf).

ATENCAO a duas armadilhas ja documentadas, ambas embutidas aqui:
  - x_t/c,max NAO pode vir do cstfoil: ele tem um erro de indexacao e devolve
    0,117 no NACA 1411, quando o valor correto e 0,30. Usamos descritores.py.
  - o c_d do nivel 1,0 de malha superestima em ~70%; a coluna de c_d traz o
    valor do nivel 1,0 (que e onde a otimizacao rodou) E o extrapolado.

Rodar com:  python tabelas_relatorio.py
'''

import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
REPO = os.path.normpath(os.path.join(AQUI, '..', '..'))
sys.path.insert(0, REPO)

import descritores as dsc
import otimiza_secao as ot
from analisa_otimos import RES, carrega, roda_euler
from compara_korn import FATOR_MALHA
from doe_clmax_xfoil import MACH_DECOLAGEM, RE_DECOLAGEM

# --- ponto de projeto da aeronave B (documento de projeto, secao 3) ---------
PONTO = {
    'h_ft': 35000.0, 'h_m': 10668.0,
    'M_inf': 0.850, 'a_inf': 296.587, 'V': 252.10,
    'rho_inf': 0.38046, 'mu_inf': 1.4451e-5, 'q_inf': 12089.69,
    'S_ref': 368.833, 'W_kgf': 229669.3,       # peso medio de cruzeiro
    'CL_aeronave': 0.5053,
    'sweep_c2_deg': 32.158, 'cos_sweep': 0.84659,
}

# corda e Reynolds por estacao, no plano normal (documento, secao 4)
GEOM = {
    'raiz':  dict(c_m=9.130, c_n=7.729, Re_n=4.343e7),
    'meio':  dict(c_m=6.899, c_n=5.841, Re_n=3.282e7),
    'ponta': dict(c_m=3.125, c_n=2.646, Re_n=1.487e7),
}


def tabela1():
    linhas = []
    for nome, est in ot.ESTACOES.items():
        g = GEOM[nome]
        linhas.append({
            'estacao': nome, 'eta': est['eta'],
            'c_ref [m]': g['c_m'], 'c_normal [m]': g['c_n'],
            'cl_ref': est['cl_ref'], '(t/c)_ref': est['tc_ref'],
            'Re_ref': g['Re_n'], 'M_n': ot.MACH_N,
        })
    return linhas


def props(Al, Au, alpha, nome_estacao, rodar_euler=True):
    d = dsc.descritores(Au, Al)
    saida = {
        'alpha [deg]': alpha * 180 / np.pi,
        '(t/c)max': d['t_max'], 'x_t/c,max': d['x_tmax'],
        '(h/c)max': d['c_max'], 'x_h/c,max': d['x_cmax'],
    }
    for i in range(len(Al)):
        saida[f'Al{i+1}'] = Al[i]
    for i in range(len(Au)):
        saida[f'Au{i+1}'] = Au[i]
    if rodar_euler:
        r = roda_euler(Al, Au, alpha, 1.0)
        saida['cl'] = r['CL']
        saida['cd'] = r['CD']
        saida['cm'] = r['CM']
        saida['cd_extrap'] = r['CD'] * FATOR_MALHA.get(nome_estacao, 1.0)
    return saida


def md(linhas, titulo):
    if not linhas:
        return ''
    cols = list(linhas[0])
    out = [f'\n### {titulo}\n', '| ' + ' | '.join(cols) + ' |',
           '|' + '---|' * len(cols)]
    for l in linhas:
        out.append('| ' + ' | '.join(
            f'{l[c]:.6g}'.replace('.', ',') if isinstance(l[c], float)
            else str(l[c]) for c in cols) + ' |')
    return '\n'.join(out)


def latex(linhas, titulo, rotulo):
    if not linhas:
        return ''
    cols = list(linhas[0])
    out = ['\\begin{table}[htb]', '\\centering',
           f'\\caption{{{titulo}}}', f'\\label{{tab:{rotulo}}}',
           '\\begin{tabular}{l' + 'r' * (len(cols) - 1) + '}', '\\hline',
           ' & '.join(c.replace('_', '\\_') for c in cols) + ' \\\\', '\\hline']
    for l in linhas:
        out.append(' & '.join(
            f'{l[c]:.5g}'.replace('.', ',') if isinstance(l[c], float)
            else str(l[c]) for c in cols) + ' \\\\')
    out += ['\\hline', '\\end{tabular}', '\\end{table}']
    return '\n'.join(out)


def main():
    partes_md = ['# Tabelas para o relatorio -- Lab 03\n',
                 '## Ponto de projeto da aeronave\n']
    for k, v in PONTO.items():
        partes_md.append(f'- **{k}**: {v}')

    t1 = tabela1()
    partes_md.append(md(t1, 'Tab. 1 -- secoes de referencia'))

    print('Rodando o Euler no perfil de partida e no otimo de cada '
          'estacao...', flush=True)
    for nome in ot.ESTACOES:
        com = carrega(f'otim_{nome}')
        if com is None:
            print(f'  {nome}: sem resultado')
            continue
        print(f'  {nome}...', flush=True)
        Ali, Aui, ai = ot.desmonta(com['xx_ini'])
        Alo, Auo, ao = ot.desmonta(com['xx_otimo'])
        linhas = [{'perfil': 'partida (NACA 1411 reesc.)',
                   **props(Ali, Aui, ai, nome)},
                  {'perfil': 'otimizado', **props(Alo, Auo, ao, nome)}]
        sem = carrega(f'otim_{nome}_sem_bluntez')
        if sem is not None:
            Als, Aus, asx = ot.desmonta(sem['xx_otimo'])
            linhas.append({'perfil': 'otimizado sem restr. de cl_max',
                           **props(Als, Aus, asx, nome)})
        partes_md.append(md(linhas, f'Tab. 2 -- estacao {nome}'))

    texto = '\n'.join(partes_md)
    cam = os.path.join(RES, 'tabelas_relatorio.md')
    with open(cam, 'w', encoding='utf-8') as fid:
        fid.write(texto + '\n')
    print(texto)
    print(f'\ngravado em {cam}')

    cam_tex = os.path.join(RES, 'tabelas_relatorio.tex')
    with open(cam_tex, 'w', encoding='utf-8') as fid:
        fid.write(latex(t1, 'Secoes de referencia', 'secoes') + '\n')
    print(f'LaTeX da Tab. 1 em {cam_tex}')


if __name__ == '__main__':
    main()
