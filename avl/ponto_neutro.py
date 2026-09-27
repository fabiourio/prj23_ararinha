'''
Comparacao do ponto neutro (NP) entre o designTool e o AVL, por etapa de
"build-up" da configuracao (asa -> +winglet -> +fuselagem -> +nacele ->
completo), em M 0,2 e M 0,85 (cruzeiro). Suporte ao item 8 do relatorio.

A margem estatica (MS_aft) nao entra na escolha de Lc_h em varredura_eh.py
(ver docstring la): o produto eta_h*S_h e mantido fixo pelo designTool, e a
contribuicao da EH ao NP praticamente nao muda com Lc_h. Este script
investiga POR QUE o NP do AVL discorda do NP do designTool (achado
reportado no item 8): a diferenca vem principalmente

  - da fuselagem: o momento de Munk (corpo esbelto) calculado pelo AVL e
    2-3x o valor empirico de Raymer/Gilruth usado no designTool (balance.py,
    Eq. 16.25), e essa diferenca e amplificada pelo fator de
    Prandtl-Glauert em M 0,85;
  - do anel da nacele (-6/-7% CMA), que o designTool nao modela;
  - o a.c. isolado da asa no AVL fica a re de c/4 (34-39% CMA), tambem
    diferente da hipotese do designTool (exatamente c/4).

A contribuicao da propria EH concorda bem entre os dois metodos (+37/+31%
CMA no AVL vs +38/+34% no designTool, M 0,85/M 0,2).

`np_designtool(Lc_h, M)` replica as formulas de designTool/balance.py (so a
parte do ponto neutro) para poder variar o Mach livremente -- balance.py
roda so no Mach_cruise do designTool. A funcao e testada contra
aeronave(Lc_h)['xnp_dt'] em tests/test_ponto_neutro.py.

Modelos parciais do AVL: escreve o .avl completo (escreve_avl) e recorta os
blocos SURFACE/BODY por nome, remontando um arquivo por etapa (arquivos
temporarios em resultados/_tmp/pn_*.avl, nomes distintos dos demais
scripts do Lab 04).

Rodar de dentro de avl/:   python ponto_neutro.py [Lc_h]
(default: Lc_h escolhido em resultados/lc_h_escolhido.json, senao 4.6)
'''

import csv
import json
import os
import re
import sys

import numpy as np

from estilo import plt, AZUL, LARANJA, VERDE, TINTA2
from aeronave import _analisa, aeronave
from analises import it_para_de_zero
from avl_run import caso
from designTool.geometry import change_sweep
from gera_avl import escreve_avl, AQUI

SAIDA = os.path.join(AQUI, 'resultados')
LC_H_FALLBACK = 4.6
MACHS = (0.2, 0.85)

# Etapas do build-up e os blocos SURFACE/BODY (nomeados em gera_avl.py) que
# cada uma acumula. O designTool nao modela winglet nem nacele: nessas
# etapas o valor do designTool fica igual ao da etapa anterior (patamar).
ETAPAS = [('asa', ['Wing']),
          ('+winglet', ['Wing', 'Winglet']),
          ('+fuselagem', ['Wing', 'Winglet', 'Fuselage']),
          ('+nacele', ['Wing', 'Winglet', 'Fuselage', 'Nacelle']),
          ('completo', ['Wing', 'Winglet', 'Fuselage', 'Nacelle',
                        'Horizontal tail', 'Vertical tail'])]


def np_designtool(Lc_h, M):
    '''
    Replica o calculo do ponto neutro de designTool/balance.py (Eq. 16.9 e
    16.23 do Raymer), mas com o Mach como parametro livre (balance.py usa
    sempre Mach_cruise). Devolve xac_w (a.c. da asa isolada), xac_wf (com a
    fuselagem) e xnp (ponto neutro completo, com a EH), tudo em metros
    (posicao absoluta em x), alem de CLa_w, CLa_h e deda.
    '''
    ap = _analisa(Lc_h)
    inp, g = ap['inputs'], ap['geometry']
    AR = ap['aerodynamics']['AR_eff']
    S_w, c = inp['S_w'], g['cm_w']
    xm = g['xm_w']

    sweep_maxt_w = change_sweep(0.25, 0.40, inp['sweep_w'], g['b_w']/2,
                                g['cr_w'], g['ct_w'])
    beta2 = 1 - M**2
    CLa_w = 2*np.pi*AR/(2 + np.sqrt(4 + AR**2*beta2/0.95**2 *
                                     (1 + np.tan(sweep_maxt_w)**2/beta2)))
    xac_w = xm + 0.25*c

    Kf = 0.1462*np.exp(4.8753*(inp['xr_w'] + 0.25*g['cr_w'])/inp['L_f'])
    CMa_f = Kf*inp['D_f']**2*inp['L_f']/c/S_w
    CLa_wf = 0.98*CLa_w
    xac_wf = xac_w - CMa_f/CLa_wf*c

    sweep_maxt_h = change_sweep(0.25, 0.40, inp['sweep_h'], g['b_h']/2,
                                g['cr_h'], g['ct_h'])
    AR_h = inp['AR_h']
    CLa_h = 2*np.pi*AR_h/(2 + np.sqrt(4 + AR_h**2*beta2/0.95**2 *
                                        (1 + np.tan(sweep_maxt_h)**2/beta2)))*0.98
    xac_h = g['xm_h'] + 0.25*g['cm_h']

    K_A = 1/AR - 1/(1 + AR**1.7)
    K_lambda = (10 - 3*inp['taper_w'])/7
    h_H = abs(g['zm_h'] - inp['zr_w'])
    L_H = Lc_h*c
    b_w = g['b_w']
    K_H = (1 - h_H/b_w)/(2*L_H/b_w)**(1/3)
    CLa_w0 = 2*np.pi*AR/(2 + np.sqrt(4 + AR**2/0.95**2*(1 + np.tan(sweep_maxt_w)**2)))
    deda = 4.44*(K_A*K_lambda*K_H*np.sqrt(np.cos(inp['sweep_w'])))**1.19*CLa_w/CLa_w0

    t = inp['eta_h']*g['S_h']/S_w*CLa_h*(1 - deda)
    xnp = (CLa_wf*xac_wf + t*xac_h)/(CLa_wf + t)
    return dict(xac_w=xac_w, xac_wf=xac_wf, xnp=xnp, CLa_w=CLa_w, CLa_h=CLa_h,
                deda=deda, AR=AR)


def _dt_por_etapa(Lc_h, M):
    '''Valor do designTool em cada etapa do build-up (patamar nas etapas
    que o designTool nao modela: winglet e nacele).'''
    d = np_designtool(Lc_h, M)
    return {'asa': d['xac_w'], '+winglet': d['xac_w'],
            '+fuselagem': d['xac_wf'], '+nacele': d['xac_wf'],
            'completo': d['xnp']}


def constroi_modelos(av):
    '''Escreve o .avl completo e recorta um arquivo por etapa do build-up
    (temporarios em resultados/_tmp/pn_*.avl). Devolve (arquivo_completo,
    {etapa: arquivo}).'''
    base = escreve_avl(av, 'resultados/_tmp/pn_full.avl', cg='aft')
    with open(os.path.join(AQUI, base)) as f:
        txt = f.read()
    blocos = re.split(r'(?=#-{20,}\n)', txt)

    def nome(bl):
        m = re.search(r'(SURFACE|BODY)\n(.+)\n', bl)
        return m.group(2).strip() if m else 'cab'

    B = {nome(bl): bl for bl in blocos}
    arquivos = {}
    for tag, blocos_incluidos in ETAPAS:
        slug = tag.replace('+', 'p').replace(' ', '_')
        p = f'resultados/_tmp/pn_{slug}.avl'
        with open(os.path.join(AQUI, p), 'w') as f:
            f.write(''.join(B[k] for k in ['cab'] + blocos_incluidos))
        arquivos[tag] = p
    return base, arquivos


def tabela(Lc_h):
    '''Monta a tabela do build-up do NP (uma linha por Mach x etapa x
    ferramenta), com a margem estatica (fwd e aft) calculada a partir do
    xnp de cada linha. As linhas com etapa == 'completo' sao as margens
    estaticas "oficiais" das tres formas descritas na docstring do
    modulo.'''
    av = aeronave(Lc_h)
    c = av['Cref']
    ap = _analisa(Lc_h)
    xm_w = ap['geometry']['xm_w']

    base, arquivos = constroi_modelos(av)
    it_aft, _ = it_para_de_zero(base, av['M'], av['CL'])

    def pc(x):
        return (x - xm_w)/c*100

    linhas = []
    for M in MACHS:
        dt = _dt_por_etapa(Lc_h, M)
        for tag, blocos_incluidos in ETAPAS:
            emp = 'Horizontal tail' in blocos_incluidos
            r = caso(arquivos[tag], M, alfa=2.0, it=it_aft if emp else 0.0,
                      de=0.0, empenagem=emp, derivadas=True)
            xnp_avl = r['xnp']
            xnp_dt = dt[tag]
            for ferramenta, xnp, CLa in (('AVL', xnp_avl, r['CLa']),
                                          ('designTool', xnp_dt, float('nan'))):
                linhas.append({
                    'M': M, 'etapa': tag, 'ferramenta': ferramenta,
                    'xnp_m': xnp, 'xnp_pctMAC': pc(xnp), 'CLa': CLa,
                    'MS_fwd': (xnp - av['xcg_fwd'])/c,
                    'MS_aft': (xnp - av['xcg_aft'])/c,
                })
    return av, it_aft, linhas


def margens(linhas):
    '''Extrai as margens estaticas (fwd, aft) das linhas 'completo': AVL em
    M 0,2 e M 0,85, e designTool (independe do Mach na tabela, mas listado
    nos dois Machs por simetria).'''
    out = {}
    for l in linhas:
        if l['etapa'] != 'completo':
            continue
        chave = f"{l['ferramenta']}_M{l['M']:g}" if l['ferramenta'] == 'AVL' else l['ferramenta']
        if chave not in out:
            out[chave] = {'MS_fwd': l['MS_fwd'], 'MS_aft': l['MS_aft']}
    return out


def imprime_tabela(linhas):
    print(f"{'M':>5} {'etapa':14} {'ferramenta':11} {'NP [% CMA]':>11} "
          f"{'MS_fwd':>8} {'MS_aft':>8}")
    for l in linhas:
        print(f"{l['M']:5.2f} {l['etapa']:14} {l['ferramenta']:11} "
              f"{l['xnp_pctMAC']:11.1f} {l['MS_fwd']:8.3f} {l['MS_aft']:8.3f}")


def figura(av, linhas):
    etapas = [e[0] for e in ETAPAS]
    n = len(etapas)
    y = np.arange(n)
    # Segmentos que o designTool nao modela (winglet, nacele): traco
    # pontilhado, para deixar claro que sao patamares (o valor nao muda),
    # nao uma medida real da etapa.
    NAO_MODELADO_DT = {('asa', '+winglet'), ('+fuselagem', '+nacele')}
    series = [('AVL', 0.85, AZUL, 'AVL, M 0,85'),
              ('AVL', 0.2, LARANJA, 'AVL, M 0,2'),
              ('designTool', 0.85, VERDE, 'designTool')]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for ferramenta, M, cor, rot in series:
        vals = [next(l['xnp_pctMAC'] for l in linhas
                     if l['etapa'] == et and l['ferramenta'] == ferramenta and l['M'] == M)
                for et in etapas]
        if ferramenta == 'designTool':
            for i in range(n - 1):
                tracejado = (etapas[i], etapas[i + 1]) in NAO_MODELADO_DT
                ax.plot(vals[i:i + 2], y[i:i + 2], color=cor, lw=1.4,
                        ls=':' if tracejado else '-')
            ax.plot(vals, y, 'o', color=cor, ms=4, label=rot)
        else:
            ax.plot(vals, y, '-o', color=cor, label=rot, lw=1.4, ms=4)

    c = av['Cref']
    # posicoes do CG em %CMA usam a mesma referencia (xm da asa) da tabela
    xnp_ref = next(l for l in linhas if l['etapa'] == 'completo'
                   and l['ferramenta'] == 'AVL' and l['M'] == 0.85)
    xm_w = xnp_ref['xnp_m'] - xnp_ref['xnp_pctMAC']/100*c
    xcg_fwd_pct = (av['xcg_fwd'] - xm_w)/c*100
    xcg_aft_pct = (av['xcg_aft'] - xm_w)/c*100
    ax.axvline(xcg_fwd_pct, color=TINTA2, lw=0.9, ls='--', zorder=0)
    ax.axvline(xcg_aft_pct, color=TINTA2, lw=0.9, ls='--', zorder=0)

    pad_topo, pad_baixo = 2.4, 0.5
    ax.set_ylim(n - 1 + pad_baixo, -pad_topo)  # asa no topo, completo embaixo
    ytxt = -pad_topo + 0.15
    caixa = dict(facecolor='white', edgecolor='none', pad=1.0)
    ax.text(xcg_fwd_pct, ytxt, 'CG dianteiro', color=TINTA2, fontsize=7,
            ha='center', va='top', rotation=90, bbox=caixa)
    ax.text(xcg_aft_pct, ytxt, 'CG traseiro', color=TINTA2, fontsize=7,
            ha='center', va='top', rotation=90, bbox=caixa)

    ax.set_yticks(y)
    ax.set_yticklabels(etapas)
    ax.set_xlabel('ponto neutro [% CMA]')
    ax.set_title(f'Ponto neutro por etapa do build-up (Lc_h = {av["Lc_h"]:.4f})')
    ax.legend(fontsize=8, loc='lower center', bbox_to_anchor=(0.5, -0.24),
              ncol=3, frameon=False)
    fig.tight_layout(rect=(0.0, 0.06, 1.0, 1.0))
    fig.savefig(os.path.join(SAIDA, 'ponto_neutro.png'))


def lc_h_default():
    p = os.path.join(SAIDA, 'lc_h_escolhido.json')
    if os.path.exists(p):
        with open(p) as f:
            d = json.load(f)
        if d.get('Lc_h') is not None:
            return float(d['Lc_h'])
    return LC_H_FALLBACK


def main():
    Lc_h = float(sys.argv[1]) if len(sys.argv) > 1 else lc_h_default()
    print(f'Lc_h = {Lc_h:.4f}')
    av, it_aft, linhas = tabela(Lc_h)
    print(f"it (compensacao no cruzeiro, modelo completo, CG traseiro): {it_aft:.3f}")
    imprime_tabela(linhas)

    m = margens(linhas)
    print('\nMargens estaticas (etapa completo):')
    for chave, v in m.items():
        print(f"  {chave:14} MS_fwd = {v['MS_fwd']:7.3f}   MS_aft = {v['MS_aft']:7.3f}")

    os.makedirs(SAIDA, exist_ok=True)
    with open(os.path.join(SAIDA, 'ponto_neutro.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)

    figura(av, linhas)


if __name__ == '__main__':
    main()
