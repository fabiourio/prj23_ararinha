'''
Figuras e resumo da investigacao da otimizacao de torcao.

Le os resultados de runs/*.json (gerados por roda_otim.py) e produz:
  convergencia_otimizacao.png  -- CDff e restricao de estol ao longo das
                                  iteracoes do SLSQP (caso de producao)
  multipartida_torcao.png      -- torcao otima de varias partidas e passos
                                  de diferenca finita, e a solucao sem a
                                  restricao de estol
  resumo.csv                   -- uma linha por rodada

Rodar de dentro de avl/investigacao_torcao/:   python figuras_investigacao.py
'''

import csv
import glob
import json
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(AQUI))
from estilo import plt, AZUL, LARANJA, VERDE, CINZA, TINTA2  # noqa: E402

ETAS = [0.0, 0.1011, 0.22, 0.32, 0.398, 0.48, 0.56, 0.70, 0.82, 0.90, 1.0]

DESCRICAO = {
    'A_prod': 'producao reproduzida (partida sem torcao, eps 0,2 grau)',
    'eps01': 'passo de diferenca finita 0,1 grau',
    'eps05': 'passo de diferenca finita 0,5 grau',
    'ms_alt': 'partida alternativa',
    'ms_washout': 'partida com washout linear',
    'full_semestol': 'sem restricao de estol',
    'full_m10': 'margem de estol de 1 grau sobre o aileron',
    'full_m20': 'margem de estol de 2 graus sobre o aileron',
    'pchip4_m10': 'spline PCHIP de 4 nos, margem 1 grau',
    'pchip4_m20': 'spline PCHIP de 4 nos, margem 2 graus',
    'pl2_m02': 'linear por partes, 2 graus de liberdade',
    'pl3_m02': 'linear por partes, 3 graus de liberdade',
    'E_asafus': 'asa + fuselagem, sem EH (sem compensacao)',
    'E_fwd': 'aeronave completa, CG dianteiro',
}


ROTULO = {'A_prod': 'produção (partida sem torção, passo 0,2°)',
          'eps01': 'passo de diferença finita 0,1°',
          'ms_alt': 'partida alternativa',
          'ms_washout': 'partida com washout linear'}


def carrega():
    runs = {}
    for f in sorted(glob.glob(os.path.join(AQUI, 'runs', '*.json'))):
        with open(f, encoding='utf-8') as fh:
            d = json.load(fh)
        if 'twist_opt' in d or 'iters' in d:
            runs[os.path.splitext(os.path.basename(f))[0]] = d
    return runs


def figura_convergencia(d):
    it = [0] + [r['it'] for r in d['iters']]
    cd = [d['CDff_base']] + [r['CDff'] for r in d['iters']]
    g = [d['evals'][0]['gmin']] + [r['gmin'] for r in d['iters']]
    ev_k = [e['k'] for e in d['evals']]
    ev_cd = [e['CDff'] for e in d['evals']]
    it_nev = [0] + [r['nev'] for r in d['iters']]

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(12.5, 3.8))
    a1.plot(ev_k, ev_cd, '.', color=CINZA, ms=3, label='cada avaliação')
    a1.plot(it_nev, cd, '-o', color=AZUL, label='fim de cada iteração')
    a1.set_xlabel('número de avaliações')
    a1.set_ylabel('$C_{D_{ff}}$ [count]')
    a1.set_title('Objetivo')
    a1.legend(fontsize=7)

    a2.axhline(0, color=TINTA2, lw=0.8, ls='--')
    a2.plot(it, g, '-o', color=LARANJA)
    a2.set_xlabel('iteração do SLSQP')
    a2.set_ylabel('menor folga da restrição [°]')
    a2.set_title('Restrição de estol (≥ 0: atendida)')
    a2.set_ylim(min(-0.5, min(g) - 0.1), max(g) + 0.2)

    snaps = [0, 1, 2, 4]
    cores = plt.cm.Blues(np.linspace(0.35, 0.8, len(snaps)))
    for s, c in zip(snaps, cores):
        tw = d['evals'][0]['twist'] if s == 0 else d['iters'][s - 1]['twist']
        a3.plot(ETAS, tw, '-o', ms=3, color=c,
                label='partida' if s == 0 else f'iteração {s}')
    a3.plot(ETAS, d['twist_opt'], '-o', ms=3, color='#08306b', label='final')
    a3.set_xlabel('η = 2y/b')
    a3.set_ylabel('torção [°]')
    a3.set_title('Evolução da torção')
    a3.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, 'convergencia_otimizacao.png'))


def figura_multipartida(runs):
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    estilos = [('A_prod', AZUL, '-'), ('eps01', LARANJA, '-'), ('ms_alt', VERDE, '-'),
               ('ms_washout', TINTA2, '-')]
    for tag, cor, ls in estilos:
        if tag in runs:
            r = runs[tag]
            ax.plot(ETAS, r['twist_opt'], ls, marker='o', ms=3, color=cor,
                    label=f'{ROTULO[tag]} ({r["CDff_opt"]:.2f} count)')
    if 'full_semestol' in runs:
        r = runs['full_semestol']
        ax.plot(ETAS, r['twist_opt'], '--', color=CINZA, marker='o', ms=3,
                label=f'sem restrição de estol ({r["CDff_opt"]:.2f} count)')
    ax.axvline(0.56, color=TINTA2, lw=0.8, ls=':')
    ax.set_xlabel('η = 2y/b')
    ax.set_ylabel('torção [°]')
    ax.set_title('Torção ótima: partidas e passos diferentes chegam ao mesmo resultado')
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(AQUI, 'multipartida_torcao.png'))


def resumo(runs):
    linhas = []
    for tag, r in runs.items():
        v = r.get('verif') or {}
        if r.get('twist_opt') is not None:
            status = 'convergiu' if r.get('sucesso') else f'nao convergiu: {r.get("mensagem")}'
            cd, tw = r.get('CDff_opt'), r.get('twist_opt')
        else:
            ult = r['iters'][-1]
            status = 'interrompida (sem convergir; restricao nao atendida)'
            cd, tw = ult['CDff'], ult['twist']
        linhas.append({
            'rodada': tag, 'descricao': DESCRICAO.get(tag, ''),
            'CDff_count': cd, 'e': v.get('e', r.get('e_opt')),
            'eta_inicio_estol': v.get('eta', r.get('eta_on_opt')),
            'alfa_estol_graus': v.get('alfa_st', r.get('alfa_on_opt')),
            'margem_aileron_graus': r.get('margem_ext_opt'),
            'folga_restricao_graus': r.get('gmin_opt'),
            'status': status,
            'torcao_graus': ' '.join(f'{t:+.2f}' for t in tw),
        })
    with open(os.path.join(AQUI, 'resumo.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)
    for l in linhas:
        print(f'{l["rodada"]:14s} CDff {l["CDff_count"]:8.2f}  {l["status"]}')


if __name__ == '__main__':
    runs = carrega()
    figura_convergencia(runs['A_prod'])
    figura_multipartida(runs)
    resumo(runs)
