'''
Item 1 da Atividade 06: diagramas V-n no inicio e no fim do cruzeiro.

Entradas aerodinamicas, todas do modelo do Lab 04 (fwd.avl, CG dianteiro):
  - CLmax: metodo da secao critica em M 0,2, compensado (tab_estol.csv);
  - CLmin: mesmo metodo com o cl minimo dos perfis (clmin_perfis.py);
  - CLa: derivada do AVL no Mach de cruzeiro, no CL de cruzeiro de cada peso,
    compensado (e a inclinacao da aeronave completa, usada na rajada).

O V-n fica na altitude de cruzeiro (velocidades em EAS); o ponto simulado
no AVL (V_A, n limite, nivel do mar) e marcado em cada diagrama. Em EAS a
curva de estol nao depende da altitude, entao o ponto e o mesmo.

Saidas em resultados/01_envelope/: envelope.json, tab_envelope.csv,
vn_cruzeiro.png.
'''

import os

import matplotlib.pyplot as plt
import numpy as np

import comum as cm
import vn
from avl_run import caso
from estilo import AZUL, CINZA, LARANJA, TINTA, TINTA2
from estol import MACH_BAIXO, _excluidas

ALFAS_NEG = (-6.0, -12.0)
ARQ_CLMIN = os.path.join(cm.DIR_ENVELOPE, 'clmin_perfis.csv')


def clmin_local(eta):
    import csv
    with open(ARQ_CLMIN, encoding='utf-8') as f:
        d = list(csv.DictReader(f))
    etas = [float(l['eta']) for l in d]
    cls = [float(l['clmin_corrigido']) for l in d]
    return np.interp(eta, etas, cls)


def clmin_aviao(av):
    '''Secao critica no estol negativo (espelho de avl/estol.py).'''
    asa = {'asa': av['asa']}
    semi = av['asa']['yt']
    r0, r1 = (caso(cm.ARQ_AVL, MACH_BAIXO, alfa=a, it=cm.IT_FWD, trim=True, faixas=True)
              for a in ALFAS_NEG)
    y = r0['faixas']['y']
    eta = y/semi
    c0, c1 = r0['faixas']['cl_norm'], r1['faixas']['cl_norm']
    b = (c1 - c0)/(ALFAS_NEG[1] - ALFAS_NEG[0])
    a = c0 - b*ALFAS_NEG[0]
    alfa_i = (clmin_local(eta) - a)/b
    cand = np.where(~_excluidas(y, asa))[0]
    i = int(cand[np.argmax(alfa_i[cand])])       # a primeira faixa a estolar
    r = caso(cm.ARQ_AVL, MACH_BAIXO, alfa=float(alfa_i[i]), it=cm.IT_FWD,
             trim=True, faixas=True)
    folga = r['faixas']['cl_norm'][cand] - clmin_local(eta[cand])
    return {'CLmin': r['CL'], 'alfa': float(alfa_i[i]), 'de': r['de'],
            'eta_crit': float(eta[i]), 'y_crit': float(y[i]),
            'folga_min': float(np.min(folga))}


def cla_cruzeiro(av, W):
    q = 0.5*av['rho_cruzeiro']*(av['M_cruzeiro']*av['a_cruzeiro'])**2
    CL = W/(q*av['S'])
    r = caso(cm.ARQ_AVL, av['M_cruzeiro'], cl=CL, it=cm.IT_FWD, trim=True, derivadas=True)
    return r['CLa'], CL


def calcula():
    av = cm.aeronave()
    neg = clmin_aviao(av)
    out = {'CLmax': cm.CLMAX, 'CLmin': neg, 'condicoes': {}}
    for k, cond in av['condicoes'].items():
        CLa, CL_cru = cla_cruzeiro(av, cond['W'])
        e = vn.envelope(cond['W'], av['S'], av['c_ref'], av['h_cruzeiro'],
                        cm.CLMAX, neg['CLmin'], CLa, av['M_cruzeiro'])
        e['CL_cruzeiro'] = CL_cru
        out['condicoes'][k] = e
    return av, out


def _kgf(W):
    return W/cm.GRAVIDADE


virgula = cm.virgula


def tabela(av, out):
    ini, fim = out['condicoes']['inicio'], out['condicoes']['fim']
    linhas = []

    def add(par, un, vi, vf, obs=''):
        linhas.append({'parametro': par, 'unidade': un, 'inicio_cruzeiro': vi,
                       'fim_cruzeiro': vf, 'observacao': obs})

    add('MTOW', 'kgf', _kgf(av['W0']), _kgf(av['W0']), 'designTool')
    add('peso da condicao', 'kgf', _kgf(ini['W']), _kgf(fim['W']),
        'MTOW x fracoes de partida, taxi, decolagem e subida (x cruzeiro no fim)')
    add('area de referencia', 'm2', av['S'], av['S'])
    add('corda de referencia (CMA)', 'm', av['c_ref'], av['c_ref'])
    add('altitude de cruzeiro', 'm', av['h_cruzeiro'], av['h_cruzeiro'])
    add('Mach de cruzeiro (M_C)', '-', av['M_cruzeiro'], av['M_cruzeiro'])
    add('Mach de mergulho (M_D)', '-', ini['M_D'], fim['M_D'], 'M_C + 0,07 (25.335(b)(2))')
    add('CLmax', '-', cm.CLMAX, cm.CLMAX, 'secao critica, M 0,2, compensado (Lab 04)')
    add('CLmin', '-', out['CLmin']['CLmin'], out['CLmin']['CLmin'],
        'secao critica negativa, M 0,2, compensado')
    add('CLa', '1/rad', ini['CLa'], fim['CLa'], 'AVL, M 0,85, CL de cruzeiro do peso')
    add('U_ref em V_C (EAS)', 'm/s', ini['U_ref'], fim['U_ref'], '25.341(a)(5)(i), na altitude de cruzeiro')
    add('U_ref em V_C (EAS)', 'ft/s', ini['U_ref']/0.3048, fim['U_ref']/0.3048)
    add('n limite positivo', '-', ini['n_max'], fim['n_max'], '25.337(b)')
    add('n limite negativo', '-', ini['n_min'], fim['n_min'], '25.337(c)')
    add('razao de massa mu', '-', ini['mu'], fim['mu'])
    add('fator de alivio Kg', '-', ini['Kg'], fim['Kg'])
    for v in ('V_S1', 'V_Sneg', 'V_A', 'V_B', 'V_C', 'V_D'):
        add(f'{v} (EAS)', 'm/s', ini[v], fim[v])
        add(f'{v} (Mach em FL350)', '-', ini['M'][v], fim['M'][v])
    for v in ('V_B', 'V_C', 'V_D'):
        add(f'n de rajada positiva em {v}', '-', ini['n_rajada_pos'][v], fim['n_rajada_pos'][v])
        add(f'n de rajada negativa em {v}', '-', ini['n_rajada_neg'][v], fim['n_rajada_neg'][v])
    add('n de estol em V_C', '-', ini['n_estol_VC'], fim['n_estol_VC'])
    add('n de estol em V_D', '-', ini['n_estol_VD'], fim['n_estol_VD'])
    add('folga V_C - (V_B + 1,32 U_ref)', 'm/s', ini['folga_VC_VB'], fim['folga_VC_VB'], '25.335(a)(1)')
    return linhas


def figura(av, out, caminho):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, k in zip(axs, ('inicio', 'fim')):
        e = out['condicoes'][k]
        cond = av['condicoes'][k]
        V, n = vn.contorno_manobra(e)
        ax.fill(V, n, color=AZUL, alpha=0.10, lw=0)
        ax.plot(V, n, color=AZUL, lw=2, label='manobra (25.333/337)')
        Vg, ng = vn.contorno_rajada(e)
        ax.plot(Vg, ng, color=LARANJA, lw=1.6, ls='--', label='rajada (25.341)')
        # curva de estol alem do envelope (ate V_A), para localizar o ponto A
        if e['V_A'] > e['V_D']:
            Vx = np.linspace(e['V_D'], e['V_A'], 50)
            ax.plot(Vx, (Vx/e['V_S1'])**2, color=AZUL, lw=1.2, ls=':')
        ax.plot([e['V_A']], [e['n_max']], marker='*', ms=13, color=TINTA, ls='none',
                label='caso simulado: $V_A$, $n$ = 2,5 (nível do mar)')
        ax.axhline(0, color=CINZA, lw=0.8)
        ax.axhline(1, color=CINZA, lw=0.6, ls=':')
        topo = 3.15
        nomes = {'V_S1': '$V_{S1}$', 'V_A': '$V_A$', 'V_B': '$V_B$', 'V_C': '$V_C$', 'V_D': '$V_D$'}
        ordem = sorted(nomes, key=lambda v: e[v])
        alin = {v: 'center' for v in ordem}
        for v0, v1 in zip(ordem, ordem[1:]):      # rotulos proximos: afasta para os lados
            if e[v1] - e[v0] < 9.0:
                alin[v0], alin[v1] = 'right', 'left'
        for v in ordem:
            ax.axvline(e[v], color=CINZA, lw=0.7, ls=':')
            ax.text(e[v], topo, f'{nomes[v]}\n{e[v]:.0f}', ha=alin[v], va='bottom',
                    fontsize=7.5, color=TINTA2)
        ax.set_title(f'{cond["nome"]}: W = {virgula(cond["W"]/cm.GRAVIDADE/1000, 1)} t, '
                     f'h = {av["h_cruzeiro"]:.0f} m', fontsize=9, pad=26)
        ax.set_xlabel('velocidade equivalente [m/s]')
        ax.set_xlim(0, 1.06*max(e['V_A'], e['V_D']))
        ax.set_ylim(-1.6, 3.1)
        if e['V_A'] > e['V_D']:
            nota = (f'em FL350, $V_A$ cai em M {virgula(e["M"]["V_A"], 2)}, além de $V_D$:\n'
                    f'o envelope é limitado pelo estol ($n$ = {virgula(e["n_estol_VD"], 2)} em $V_D$)')
        else:
            nota = f'em FL350, $V_A$ cai em M {virgula(e["M"]["V_A"], 2)}'
        ax.text(0.02, 0.04, nota, transform=ax.transAxes, fontsize=8, color=TINTA2)
    axs[0].set_ylabel('fator de carga $n$')
    axs[1].legend(loc='lower right', fontsize=7.5)
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def main():
    av, out = calcula()
    cm.salva_json(out, os.path.join(cm.DIR_ENVELOPE, 'envelope.json'))
    cm.salva_csv(tabela(av, out), os.path.join(cm.DIR_ENVELOPE, 'tab_envelope.csv'))
    figura(av, out, os.path.join(cm.DIR_ENVELOPE, 'vn_cruzeiro.png'))
    neg = out['CLmin']
    print(f'CLmax = {cm.CLMAX:.4f}  CLmin = {neg["CLmin"]:.4f} (alfa {neg["alfa"]:.2f} graus, '
          f'eta crit {neg["eta_crit"]:.3f}, de {neg["de"]:.2f})')
    for k, e in out['condicoes'].items():
        print(f'{k:6s} W = {e["W"]/cm.GRAVIDADE:9.1f} kgf  CLa = {e["CLa"]:.3f}  '
              f'VS1 = {e["V_S1"]:.1f}  VA = {e["V_A"]:.1f} (M {e["M"]["V_A"]:.3f})  '
              f'VB = {e["V_B"]:.1f}  VC = {e["V_C"]:.1f}  VD = {e["V_D"]:.1f}  '
              f'n_raj(VC) = {e["n_rajada_pos"]["V_C"]:.2f}  n_estol(VD) = {e["n_estol_VD"]:.2f}')


if __name__ == '__main__':
    main()
