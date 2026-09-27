'''
Lab 04 -- roteiro da secao 2 (itens 1 a 8), sobre os arquivos fwd.avl e
aft.avl.

Saidas em resultados/roteiro/ (tabelas .csv, figuras .png e notas_roteiro.txt
com as escolhas e as discussoes geradas automaticamente):

  1. tab_modelo.csv            referencias, CDp, malha, superficies, controles
  2. tab1_ponto_projeto.csv    Tabela 1 (ponto de projeto)
  3. tab_incidencia.csv        it que zera o profundor no cruzeiro (fwd, aft)
     trefftz_fwd.png / trefftz_aft.png   distribuicoes no plano de Trefftz
     (faixas do fs de todas as superficies); se o Ghostscript do MiKTeX
     estiver disponivel, tambem a copia do proprio grafico do AVL
     (trefftz_*_avl.png, a partir do plot.ps do AVL)
  4. tab_estol.csv, estol_cl_y.png     secao critica em M 0,2 (4 casos)
  5. polares.csv, polar_cd_cl.png, tab_polar_ponto_projeto.csv
  6. cl_alpha.png
  7. cl_delta_e.png
  8. tab_ponto_neutro.csv

Escolhas (tambem gravadas em notas_roteiro.txt):
  - Polares, CL x alfa e CL x delta_e no Mach de cruzeiro (0,85), na
    altitude de projeto; it de cada CG = o do item 3.
  - CLmax de cada caso pelo metodo da secao critica em M 0,2 (item 4), com
    o clmax dos perfis do Lab 03 (estol.py). Cada polar vai do alfa que da
    CL = -0,5 ate o alfa que da CL = CLmax do caso correspondente do item 4
    (o enunciado pede que "cada polar termine no CLmax"): como a polar e
    calculada em M 0,85, o alfa final e menor que o alfa_max de M 0,2
    (efeito de compressibilidade de Prandtl-Glauert no CL_alfa).
    Sem compensacao -> CLmax sem trimagem; compensada -> CLmax com trimagem.
  - Passo de alfa de 0,5 grau (mais os extremos exatos).

Uso (de dentro de avl/):
  python lab04_roteiro.py                       # fwd.avl e aft.avl
  python lab04_roteiro.py --gera-dev            # escreve e usa copias em
                                                # resultados/_tmp/rot_*.avl
  python lab04_roteiro.py --fwd X.avl --aft Y.avl --saida resultados/outra
'''

import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import sys

import numpy as np

from estilo import plt, AZUL, LARANJA, VERDE, TINTA, TINTA2, CINZA
from aeronave import aeronave, _analisa
from analises import it_para_de_zero
from avl_run import caso, roda
from avl_saida import (casos_ft, superficies, faixas_todas, cabecalho_avl,
                       secoes_superficie, sessao, derivadas)
from estol import estol, clmax_local, MACH_BAIXO
from gera_avl import (escreve_avl, malha_adotada, AQUI, ETA_AILERON,
                      TAPER_WINGLET)
from ponto_neutro import np_designtool
from designTool.aerodynamics import aerodynamics

CL_MIN = -0.5
PASSO_ALFA = 0.5
DE_LIMITE = 25.0            # [graus] autoridade tipica do profundor
TOL_DE = 0.01               # [graus] profundor residual aceito no item 3

CASOS = [  # (chave, cg, trim, rotulo, cor, estilo de linha)
    ('a', 'fwd', False, '(a) CG dianteiro, δe = 0', AZUL, '--'),
    ('b', 'aft', False, '(b) CG traseiro, δe = 0', LARANJA, '--'),
    ('c', 'fwd', True, '(c) CG dianteiro, compensado (Cm = 0)', AZUL, '-'),
    ('d', 'aft', True, '(d) CG traseiro, compensado (Cm = 0)', LARANJA, '-'),
]
NOME_CG = {'fwd': 'dianteiro', 'aft': 'traseiro'}


# --------------------------------------------------------------------------
# utilidades
# --------------------------------------------------------------------------

def lc_h_escolhido():
    with open(os.path.join(AQUI, 'resultados', 'lc_h_escolhido.json')) as f:
        return float(json.load(f)['Lc_h'])


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        if v == 0:
            v = 0.0
        return f'{v:.6g}' if abs(v) < 1e5 else f'{v:.1f}'
    return v


def _v(x, n=2):
    '''Numero com virgula decimal (rotulos das figuras).'''
    return f'{x:.{n}f}'.replace('.', ',')


def grava_csv(caminho, linhas, campos=None):
    campos = campos or list(linhas[0].keys())
    with open(caminho, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for ln in linhas:
            w.writerow({k: _fmt(ln.get(k, '')) for k in campos})


def imprime(titulo, linhas):
    print(f'\n== {titulo}')
    campos = list(linhas[0].keys())
    larg = {k: max(len(k), *(len(str(_fmt(l.get(k, '')))) for l in linhas)) for k in campos}
    print('  '.join(k.ljust(larg[k]) for k in campos))
    for l in linhas:
        print('  '.join(str(_fmt(l.get(k, ''))).ljust(larg[k]) for k in campos))


# --------------------------------------------------------------------------
# item 1 -- modelo
# --------------------------------------------------------------------------

def item1(arqs, saida):
    cab = {cg: cabecalho_avl(os.path.join(AQUI, a)) for cg, a in arqs.items()}
    s = sessao(arqs['aft'], 0.85, 0.0, ['a a 2', 'x', 'fs', ''])
    sup = superficies(s)
    tw = secoes_superficie(os.path.join(AQUI, arqs['aft']), 'Wing')
    tw_fwd = secoes_superficie(os.path.join(AQUI, arqs['fwd']), 'Wing')
    malha = malha_adotada()
    semi = tw[-1, 1]

    linhas = []

    def add(par, fwd, aft, unid='', obs=''):
        linhas.append({'parametro': par, 'fwd': fwd, 'aft': aft,
                       'unidade': unid, 'observacao': obs})

    add('arquivo', arqs['fwd'], arqs['aft'])
    for k, u in (('Mach', '-'), ('Sref', 'm2'), ('Cref', 'm'), ('Bref', 'm'),
                 ('Xref', 'm'), ('Zref', 'm'), ('CDp', '-')):
        obs = {'CDp': 'CD0 do designTool no ponto de projeto',
               'Xref': 'CG dianteiro / traseiro do designTool (origem no nariz)',
               'Sref': 'area de referencia do designTool (asa trapezoidal)'}.get(k, '')
        add(k, cab['fwd'][k], cab['aft'][k], u, obs)
    for sp in sup:
        add(f'superficie {sp["numero"]}: {sp["nome"]}', sp['area'], sp['area'], 'm2',
            f'area integrada pelo AVL; malha {sp["n_corda"]} (corda) x {sp["n_env"]} (envergadura)')
    add('fuselagem (BODY)', f'nbody = {malha["nbody"]}', f'nbody = {malha["nbody"]}', '',
        'fuselage_reamostrada.dat (BFILE)')
    add('torcao da asa (Ainc por secao)',
        ' '.join(f'{v:.3f}' for v in tw_fwd[:, 4]), ' '.join(f'{v:.3f}' for v in tw[:, 4]),
        'graus', 'eta = ' + ' '.join(f'{y/semi:.3f}' for y in tw[:, 1]))
    add('d1', 'aileron', 'aileron', '',
        f'asa, eta {ETA_AILERON[0]:.2f}-{ETA_AILERON[1]:.2f}, charneira em 73% da corda, SgnDup -1')
    add('d2', 'elevator', 'elevator', '', 'EH inteira, charneira em 70% da corda, SgnDup +1')
    add('d3', 'rudder', 'rudder', '', 'EV inteira, charneira em 70% da corda')
    add('variavel de projeto 1', 'it', 'it', 'graus', 'incidencia da EH (DESIGN it 1.0)')
    add('perfis da asa', 'Lab 03 (airfoils/)', 'Lab 03 (airfoils/)', '', 'EH e EV: NACA 0010')
    grava_csv(os.path.join(saida, 'tab_modelo.csv'), linhas)
    imprime('Item 1 -- modelo', linhas)
    return {'semi': semi, 'ct': tw[-1, 3]}


# --------------------------------------------------------------------------
# item 2 -- ponto de projeto
# --------------------------------------------------------------------------

def item2(av, saida):
    linhas = [
        {'parametro': 'W0', 'explicacao': 'peso maximo de decolagem [N]', 'valor': av['W0']},
        {'parametro': 'W', 'explicacao': 'peso no ponto de projeto [N]', 'valor': av['W']},
        {'parametro': 'h', 'explicacao': 'altitude [m]', 'valor': av['h']},
        {'parametro': 'rho', 'explicacao': 'densidade do ar [kg/m3]', 'valor': av['rho']},
        {'parametro': 'a', 'explicacao': 'velocidade do som [m/s]', 'valor': av['a']},
        {'parametro': 'M', 'explicacao': 'Mach [-]', 'valor': av['M']},
        {'parametro': 'V', 'explicacao': 'velocidade [m/s]', 'valor': av['V']},
        {'parametro': 'CL', 'explicacao': 'CL de projeto [-]', 'valor': av['CL']},
        {'parametro': 'Sref', 'explicacao': 'area de referencia [m2]', 'valor': av['Sref']},
    ]
    for l in linhas:
        l['observacao'] = ''
    linhas[1]['observacao'] = (f'fracao de combustivel {av["fuel_frac"]:.4f} (mesmo ponto do Lab 03), '
                               f'100% de carga paga; W = {av["W_kgf"]:.1f} kgf')
    linhas[0]['observacao'] = f'{av["W0_kgf"]:.1f} kgf; Lc_h = {av["Lc_h"]:.4f}'
    grava_csv(os.path.join(saida, 'tab1_ponto_projeto.csv'), linhas)
    imprime('Item 2 -- Tabela 1', linhas)


# --------------------------------------------------------------------------
# item 3 -- incidencia da EH e plano de Trefftz
# --------------------------------------------------------------------------

def _mgs():
    exe = shutil.which('mgs') or shutil.which('gswin64c') or shutil.which('gs')
    if exe:
        return exe
    cand = glob.glob(os.path.expanduser(r'~\AppData\Local\Programs\MiKTeX\miktex\bin\x64\mgs.exe'))
    return cand[0] if cand else None


def trefftz_avl(arquivo, it, av, destino_base):
    '''
    Grafico de Trefftz do proprio AVL (hardcopy PostScript), sem janela
    (PLOP G desliga a tela). O AVL grava plot.ps na pasta avl/; o arquivo e
    movido para destino_base.ps e, se houver Ghostscript, convertido em PNG
    (girado para a horizontal). Opcional: falhas so geram aviso.
    '''
    ps_avl = os.path.join(AQUI, 'plot.ps')
    if os.path.exists(ps_avl):
        os.remove(ps_avl)
    cmds = ['plop', 'g', '', f'load {arquivo}', 'oper', 'm', f'mn {av["M"]}', '',
            'de', f'1 {it}', '', 'd2 pm 0', f'a c {av["CL"]}', 'x', 't', 'h', '', '', 'quit']
    try:
        roda(cmds, timeout=120)
    except Exception as e:  # noqa: BLE001
        print(f'  aviso: hardcopy do Trefftz falhou ({e})')
        return None
    if not os.path.exists(ps_avl):
        print('  aviso: o AVL nao gravou plot.ps')
        return None
    ps = destino_base + '.ps'
    shutil.move(ps_avl, ps)
    exe = _mgs()
    if not exe:
        return ps
    png = destino_base + '.png'
    subprocess.run([exe, '-q', '-dNOPAUSE', '-dBATCH', '-sDEVICE=png16m', '-r150',
                    f'-sOutputFile={png}', ps], capture_output=True, timeout=120)
    if os.path.exists(png):
        from PIL import Image
        with Image.open(png) as im:
            im.rotate(-90, expand=True).save(png)
        os.remove(ps)
        return png
    return ps


def _winglet_desdobrado(f, y_ponta, ct):
    '''Winglet vertical: todas as faixas tem o mesmo Yle. Para o grafico,
    usa a coordenada "desdobrada" s = y_ponta + altura, com a altura tirada
    da corda local (corda linear de ct a TAPER*ct ao longo da altura ct).'''
    altura = (ct - f['corda'])/(1.0 - TAPER_WINGLET)
    return y_ponta + altura


def figura_trefftz(av, geo, r, it, cg, destino):
    fx = r['faixas_todas']
    asa, wl, eh = fx['Wing'], fx['Winglet'], fx['Horizontal tail']
    cref = av['Cref']
    y_ponta = geo['semi']
    s_wl = _winglet_desdobrado(wl, y_ponta, geo['ct'])

    fig, axs = plt.subplots(3, 1, figsize=(7.2, 8.0), sharex=True)
    series = [('Asa', asa['y'], asa, AZUL), ('Winglet (desdobrado)', s_wl, wl, VERDE),
              ('EH', eh['y'], eh, LARANJA)]
    campos = [('ccl', 'c·cl / cref', True), ('cl_norm', 'cl⊥ (cl_norm)', False),
              ('cl', 'cl', False)]
    for ax, (campo, rot, divide) in zip(axs, campos):
        for nome, y, f, cor in series:
            v = f[campo]/cref if divide else f[campo]
            ax.plot(y, v, color=cor, lw=1.4, label=nome)
        ax.axvline(y_ponta, color=CINZA, lw=0.8, ls=':')
        ax.axhline(0.0, color=CINZA, lw=0.6)
        ax.set_ylabel(rot)
    axs[0].legend(fontsize=8, loc='upper right')
    axs[0].text(y_ponta, axs[0].get_ylim()[1], ' junção asa-winglet', color=TINTA2,
                fontsize=7, va='top', ha='left')
    axs[-1].set_xlabel('y [m]  (semiasa direita; winglet: y da ponta + altura)')
    txt = (f'M = {av["M"]:.2f}   α = {r["alfa"]:.3f}°   it = {it:.3f}°   '
           f'δe = {r["de"]:.3f}°\n'
           f'CL = {r["CL"]:.4f}   CDi (Trefftz) = {1e4*r["CDff"]:.2f} counts   '
           f'e = {r["e"]:.4f}   Cm = {r["Cm"]:.5f}')
    axs[0].set_title(f'Plano de Trefftz, CG {NOME_CG[cg]}, cruzeiro\n{txt}', fontsize=8.5)
    fig.tight_layout()
    fig.savefig(destino)
    plt.close(fig)


def item3(arqs, av, geo, saida):
    linhas, its = [], {}
    for cg in ('fwd', 'aft'):
        it, _ = it_para_de_zero(arqs[cg], av['M'], av['CL'], tol=TOL_DE)
        r = caso(arqs[cg], av['M'], cl=av['CL'], it=it, trim=True, faixas=True)
        r['faixas_todas'] = faixas_todas(r['saida'])
        its[cg] = it
        linhas.append({'cg': cg, 'arquivo': arqs[cg], 'it_graus': it, 'de_residual_graus': r['de'],
                       'alfa_graus': r['alfa'], 'CL': r['CL'], 'CD_trefftz': r['CDvis'] + r['CDff'],
                       'CDp': r['CDvis'], 'CDff': r['CDff'], 'e': r['e'],
                       'CD_campo_proximo': r['CD'], 'CDind': r['CDind'], 'Cm': r['Cm'],
                       'M': av['M']})
        figura_trefftz(av, geo, r, it, cg, os.path.join(saida, f'trefftz_{cg}.png'))
        p = trefftz_avl(arqs[cg], it, av, os.path.join(saida, f'trefftz_{cg}_avl'))
        if p:
            print(f'  hardcopy do AVL: {os.path.relpath(p, AQUI)}')
    grava_csv(os.path.join(saida, 'tab_incidencia.csv'), linhas)
    imprime('Item 3 -- incidencia da EH que zera o profundor', linhas)
    return its


# --------------------------------------------------------------------------
# item 4 -- estol (secao critica, M 0,2)
# --------------------------------------------------------------------------

def item4(arqs, av, its, saida):
    semi = av['asa']['yt']
    res, linhas = {}, []
    for cg in ('fwd', 'aft'):
        for trim in (False, True):
            e = estol(arqs[cg], it=its[cg], trim=trim, av=av)
            res[(cg, trim)] = e
            dentro = e['eta_crit'] < ETA_AILERON[0]
            linhas.append({'cg': cg, 'trimagem': 'com' if trim else 'sem',
                           'it_graus': its[cg], 'M': MACH_BAIXO, 'alfa_max_graus': e['alfa'],
                           'CLmax': e['CL'], 'de_graus': e['de'],
                           'eta_crit': e['eta_crit'], 'y_crit_m': e['eta_crit']*semi,
                           'estol_inicia_antes_do_aileron': 'sim' if dentro else 'nao',
                           'excesso_max_cl_norm': e['excesso_max']})
    grava_csv(os.path.join(saida, 'tab_estol.csv'), linhas)
    imprime('Item 4 -- secao critica (M 0,2)', linhas)

    fig, axs = plt.subplots(1, 2, figsize=(10.0, 4.2), sharey=True)
    for ax, cg in zip(axs, ('fwd', 'aft')):
        eta_ref = res[(cg, False)]['eta']
        y = eta_ref*semi
        ax.plot(y, clmax_local(eta_ref), color=TINTA2, lw=1.2, ls='--',
                label='cl_max do perfil (Lab 03)')
        for trim, ls, lw in ((True, '-', 2.4), (False, '--', 1.2)):
            e = res[(cg, trim)]
            cor = AZUL if not trim else LARANJA
            excl = np.isin(e['eta'], e['excluidas'])
            ye = e['eta']*semi
            rot = (f'{"com" if trim else "sem"} trimagem: α = {e["alfa"]:.2f}°, '
                   f'CLmax = {e["CL"]:.3f}')
            ax.plot(ye[~excl], e['cl_norm'][~excl], color=cor, lw=lw, ls=ls, label=rot)
            ax.scatter(ye[excl], e['cl_norm'][excl], color=CINZA, marker='x', s=14, zorder=3,
                       label='faixa excluída (raiz/junção)' if not trim else None)
            i = int(np.argmin(np.abs(e['eta'] - e['eta_crit'])))
            ax.plot(ye[i], e['cl_norm'][i], marker='*', ms=11, color=cor, mec=TINTA,
                    mew=0.6, ls='none', zorder=4)
        ax.axvline(ETA_AILERON[0]*semi, color=VERDE, lw=0.9, ls=':')
        ax.text(ETA_AILERON[0]*semi, 0.05, ' raiz do aileron\n (η = 0,56)', color=VERDE,
                fontsize=7, va='bottom')
        ax.set_xlabel('y [m]')
        ax.set_title(f'CG {NOME_CG[cg]} (it = {its[cg]:.2f}°), M 0,2 no estol\n'
                     '★ = início do estol', fontsize=9)
        ax.set_ylim(0.0, 2.1)
        ax.legend(fontsize=7, loc='lower left', bbox_to_anchor=(0.0, 0.12))
    axs[0].set_ylabel('cl_norm')
    fig.tight_layout()
    fig.savefig(os.path.join(saida, 'estol_cl_y.png'))
    plt.close(fig)

    notas = ['Item 4 -- caracteristicas de estol (metodo da secao critica, M 0,2):']
    for l in linhas:
        pos = 'para dentro' if l['estol_inicia_antes_do_aileron'] == 'sim' else 'sobre/para fora'
        notas.append(
            f'  CG {NOME_CG[l["cg"]]}, {l["trimagem"]} trimagem: alfa_max = {l["alfa_max_graus"]:.2f} graus, '
            f'CLmax = {l["CLmax"]:.3f}, de = {l["de_graus"]:.2f} graus; estol inicia em '
            f'eta = {l["eta_crit"]:.3f} (y = {l["y_crit_m"]:.2f} m), {pos} da raiz do aileron '
            f'(eta = {ETA_AILERON[0]:.2f}).')
    ok = all(l['estol_inicia_antes_do_aileron'] == 'sim' for l in linhas)
    notas.append('  Conclusao: ' + (
        'nos quatro casos o estol comeca na parte interna da asa, antes da raiz do aileron, '
        'como exigido. A margem, porem, e pequena: o inicio fica so ~2 m para dentro da raiz '
        'do aileron e toda a asa externa esta a poucos centesimos do clmax no estol '
        '(estol quase simultaneo na asa externa, efeito da restricao com desempate de 0,2 '
        'grau na otimizacao de torcao). Criterio atendido com margem pequena.' if ok else
        'em pelo menos um caso o estol comeca sobre o aileron ou para fora dele: perda de '
        'controle lateral no estol e tendencia a cair de asa -- caracteristica NAO adequada '
        '(corrigir com torcao/perfis de ponta).'))
    notas.append('  As faixas da raiz e da juncao com o winglet sao excluidas da busca '
                 '(singularidades de canto do VLM; ver estol.py).')
    return res, notas


# --------------------------------------------------------------------------
# itens 5, 6, 7 -- polares, CL x alfa, CL x de (M de cruzeiro)
# --------------------------------------------------------------------------

def varredura(arquivo, mach, it, trim, alfas):
    cmds = ['d2 pm 0' if trim else 'd2 d2 0']
    for a in alfas:
        cmds += [f'a a {a:.6f}', 'x']
    s = sessao(arquivo, mach, it, cmds)
    if 'Cannot trim' in s:
        raise RuntimeError(f'AVL nao compensou a arfagem em {arquivo}')
    rs = casos_ft(s)
    if len(rs) != len(alfas):
        raise RuntimeError(f'{len(rs)} casos lidos, {len(alfas)} pedidos ({arquivo})')
    for a, r in zip(alfas, rs):
        if abs(r['alfa'] - a) > 1e-4:
            raise RuntimeError(f'alfa lido {r["alfa"]} != pedido {a}')
        if trim and abs(r['Cm']) > 1e-4:
            raise RuntimeError(f'Cm = {r["Cm"]} com compensacao (alfa {a})')
        if not trim:
            r['de'] = 0.0
        r['CDT'] = r['CDvis'] + r['CDff']
    return rs


def itens567(arqs, av, its, estol_res, saida):
    M, CLp = av['M'], av['CL']
    ap = _analisa(av['Lc_h'])
    CD_dt, _, dd = aerodynamics(ap, Mach=M, altitude=av['h'], CL=CLp)

    pol, proj, linhas_pol, linhas_pp = {}, {}, [], []
    for chave, cg, trim, rot, cor, ls in CASOS:
        it = its[cg]
        CLmax = estol_res[(cg, trim)]['CL']
        r0 = caso(arqs[cg], M, cl=CL_MIN, it=it, trim=trim)
        r1 = caso(arqs[cg], M, cl=CLmax, it=it, trim=trim)
        a0, a1 = r0['alfa'], r1['alfa']
        meio = np.arange(np.ceil(a0/PASSO_ALFA)*PASSO_ALFA, a1, PASSO_ALFA)
        meio = meio[(meio - a0 > 0.05) & (a1 - meio > 0.05)]
        alfas = np.concatenate([[a0], meio, [a1]])
        rs = varredura(arqs[cg], M, it, trim, alfas)
        pol[chave] = {k: np.array([r[k] for r in rs]) for k in ('alfa', 'CL', 'CD', 'CDT',
                                                                'CDind', 'CDff', 'Cm', 'de')}
        for r in rs:
            linhas_pol.append({'caso': chave, 'cg': cg, 'trimagem': 'com' if trim else 'sem',
                               'M': M, 'it_graus': it, 'alfa_graus': r['alfa'], 'CL': r['CL'],
                               'CD_trefftz': r['CDT'], 'CD_campo_proximo': r['CD'],
                               'CDp': r['CDvis'], 'CDff': r['CDff'], 'CDind': r['CDind'],
                               'Cm': r['Cm'], 'de_graus': r['de']})
        rp = caso(arqs[cg], M, cl=CLp, it=it, trim=trim)
        rp['de'] = rp['de'] if trim else 0.0
        rp['CDT'] = rp['CDvis'] + rp['CDff']
        proj[chave] = rp
        CD_interp = float(np.interp(CLp, pol[chave]['CL'], pol[chave]['CDT']))
        linhas_pp.append({'caso': chave, 'descricao': rot, 'M': M, 'CL_projeto': CLp,
                          'alfa_graus': rp['alfa'], 'de_graus': rp['de'],
                          'CD_AVL': rp['CDT'], 'CD_AVL_interp_polar': CD_interp,
                          'CDp': rp['CDvis'], 'CDff': rp['CDff'],
                          'CD_AVL_campo_proximo': rp['CD'], 'CD_designTool': CD_dt,
                          'dif_counts': 1e4*(rp['CDT'] - CD_dt),
                          'dif_counts_campo_proximo': 1e4*(rp['CD'] - CD_dt),
                          'CLmax_fim_polar': CLmax, 'alfa_fim_polar_graus': a1})
    grava_csv(os.path.join(saida, 'polares.csv'), linhas_pol)
    grava_csv(os.path.join(saida, 'tab_polar_ponto_projeto.csv'), linhas_pp)
    imprime('Item 5 -- CD no CL de projeto (M 0,85)', linhas_pp)
    print(f'  designTool: CD = {CD_dt:.5f} (CD0 {dd["CD0"]:.5f} + CDind {dd["CDind"]:.5f} '
          f'+ CDwave {dd["CDwave"]:.5f})')

    # polar do designTool para referencia
    cl_dt = np.linspace(CL_MIN, max(e['CL'] for e in estol_res.values()), 60)
    cd_dt = np.array([aerodynamics(ap, Mach=M, altitude=av['h'], CL=c)[0] for c in cl_dt])

    # --- item 5: CD x CL
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(cl_dt, cd_dt, color=CINZA, lw=1.2, ls=':', label='designTool (M 0,85)')
    for chave, cg, trim, rot, cor, ls in CASOS:
        p = pol[chave]
        ax.plot(p['CL'], p['CDT'], color=cor, ls=ls, lw=1.4, label=rot)
        ax.plot(CLp, proj[chave]['CDT'], 'o', color=cor, mec=TINTA, mew=0.5, ms=5, zorder=4)
    ax.plot([], [], 'o', color='white', mec=TINTA, mew=0.6, label='ponto de projeto')
    ax.set_xlabel('CL')
    ax.set_ylabel('CD = CDp + CDi (Trefftz)')
    ax.set_title(f'Polares de arrasto, M {_v(M)}, h = {av["h"]:.0f} m '
                 f'(fim de cada curva: CLmax do item 4)', fontsize=9)
    ax.legend(fontsize=7.5, loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(saida, 'polar_cd_cl.png'))
    plt.close(fig)

    # --- item 6: CL x alfa
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for chave, cg, trim, rot, cor, ls in CASOS:
        p = pol[chave]
        ax.plot(p['alfa'], p['CL'], color=cor, ls=ls, lw=1.4, label=rot)
        ax.plot(proj[chave]['alfa'], CLp, 'o', color=cor, mec=TINTA, mew=0.5, ms=5, zorder=4)
    ax.plot([], [], 'o', color='white', mec=TINTA, mew=0.6, label='ponto de projeto')
    ax.axhline(CLp, color=CINZA, lw=0.7, ls=':')
    ax.set_xlabel('α [°]')
    ax.set_ylabel('CL')
    ax.set_title(f'Curvas de sustentação, M {_v(M)}', fontsize=9)
    ax.legend(fontsize=7.5, loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(saida, 'cl_alpha.png'))
    plt.close(fig)

    # --- item 7: CL x de
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for chave, cg, trim, rot, cor, ls in CASOS:
        p = pol[chave]
        if trim:
            ax.plot(p['de'], p['CL'], color=cor, ls=ls, lw=1.4, label=rot)
        else:
            dash = (0, (5, 3)) if cg == 'fwd' else (4, (5, 3))
            ax.plot(p['de'], p['CL'], color=cor, ls=dash, lw=1.4,
                    label=rot + ' (sobrepostas)')
        ax.plot(proj[chave]['de'], CLp, 'o', color=cor, mec=TINTA, mew=0.5, ms=5, zorder=4)
    ax.plot([], [], 'o', color='white', mec=TINTA, mew=0.6, label='ponto de projeto')
    ax.set_xlabel('δe [°]  (positivo = bordo de fuga para baixo)')
    ax.set_ylabel('CL')
    ax.set_title(f'Deflexão de profundor, M {_v(M)} (it do item 3 em cada CG)', fontsize=9)
    ax.legend(fontsize=7.5, loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(saida, 'cl_delta_e.png'))
    plt.close(fig)

    notas = ['Item 5 -- CD no CL de projeto:']
    for l in linhas_pp:
        notas.append(f'  ({l["caso"]}) CD_AVL = {l["CD_AVL"]:.5f}  designTool = {l["CD_designTool"]:.5f}'
                     f'  diferenca = {l["dif_counts"]:+.1f} counts')
    notas.append('  CD do AVL = CDp (CD0 do designTool, constante) + CDi do plano de Trefftz (CDff). '
                 'O CDi de campo proximo (CDind, que entra no CDtot do ft) nao e usado: com a '
                 'fuselagem e os aneis das naceles ele fica irrealisticamente baixo (e implicito > 1); '
                 'o CDtot fica em tab_polar_ponto_projeto.csv/polares.csv so para referencia.')
    notas.append('  O designTool '
                 f'soma ainda o arrasto de onda (CDwave = {dd["CDwave"]:.5f} no ponto de projeto) e '
                 f'usa o seu proprio fator de Oswald (e = {dd["e"]:.4f}).')
    notas.append('Item 7 -- deflexoes de profundor (casos compensados):')
    for chave, cg, trim, rot, cor, ls in CASOS:
        if not trim:
            continue
        p = pol[chave]
        dmax = float(np.max(np.abs(p['de'])))
        incl = np.polyfit(p['CL'], p['de'], 1)[0]
        notas.append(f'  ({chave}) CG {NOME_CG[cg]}: d(de)/dCL = {incl:+.2f} graus '
                     + ('(negativo: gradiente de profundor normal, estavel)' if incl < 0 else
                        '(positivo: gradiente invertido -- margem estatica negativa no Mach de '
                        'cruzeiro, ver item 8)') + '.')
        notas.append(f'  ({chave}) CG {NOME_CG[cg]}: de de {p["de"].min():+.2f} a {p["de"].max():+.2f} '
                     f'graus entre CL = {p["CL"].min():.2f} e {p["CL"].max():.3f}; '
                     f'no ponto de projeto de = {proj[chave]["de"]:+.3f} graus; '
                     f'|de|max = {dmax:.2f} graus = {100*dmax/DE_LIMITE:.0f}% de {DE_LIMITE:.0f} graus '
                     f'(amplitude {"adequada" if dmax < 0.8*DE_LIMITE else "com margem pequena"}'
                     + ('' if incl < 0 else
                        '; porem o gradiente invertido indica aeronave estaticamente instavel '
                        'no cruzeiro com este CG, segundo o AVL') + ').')
    notas.append(f'  Referencia: autoridade tipica de profundor +-{DE_LIMITE:.0f} graus; o uso '
                 'de ate ~80% dela no envelope de cruzeiro deixa margem para manobra e rajada.')
    return pol, proj, notas


# --------------------------------------------------------------------------
# item 8 -- ponto neutro e margem estatica
# --------------------------------------------------------------------------

def item8(arqs, av, its, saida):
    ap = _analisa(av['Lc_h'])
    xm_w, c = ap['geometry']['xm_w'], av['Cref']
    xcg = {'fwd': av['xcg_fwd'], 'aft': av['xcg_aft']}
    pc = lambda x: 100*(x - xm_w)/c  # noqa: E731
    linhas = []
    for M in (av['M'], MACH_BAIXO):
        for cg in ('fwd', 'aft'):
            r = caso(arqs[cg], M, cl=av['CL'], it=its[cg], trim=True, derivadas=True)
            xnp = derivadas(r['saida'], 'st')['Xnp']
            linhas.append({'fonte': 'AVL', 'M': M, 'cg': cg, 'arquivo': arqs[cg],
                           'xnp_m': xnp, 'xnp_pctCMA': pc(xnp), 'xcg_m': xcg[cg],
                           'xcg_pctCMA': pc(xcg[cg]), 'MS': (xnp - xcg[cg])/c,
                           'observacao': f'st no CL de projeto, it = {its[cg]:.3f}, compensado'})
        xnp_dt = av['xnp_dt'] if M == av['M'] else np_designtool(av['Lc_h'], M)['xnp']
        for cg in ('fwd', 'aft'):
            linhas.append({'fonte': 'designTool', 'M': M, 'cg': cg, 'arquivo': '',
                           'xnp_m': xnp_dt, 'xnp_pctCMA': pc(xnp_dt), 'xcg_m': xcg[cg],
                           'xcg_pctCMA': pc(xcg[cg]), 'MS': (xnp_dt - xcg[cg])/c,
                           'observacao': 'balance.py' if M == av['M'] else
                           'formulas de balance.py com Mach livre (ponto_neutro.np_designtool)'})
    grava_csv(os.path.join(saida, 'tab_ponto_neutro.csv'), linhas)
    imprime('Item 8 -- ponto neutro e margem estatica', linhas)
    return linhas


# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--fwd', default='fwd.avl')
    ap.add_argument('--aft', default='aft.avl')
    ap.add_argument('--saida', default=os.path.join('resultados', 'roteiro'))
    ap.add_argument('--lc', type=float, default=None)
    ap.add_argument('--gera-dev', action='store_true',
                    help='escreve resultados/_tmp/rot_fwd.avl e rot_aft.avl (sem torcao) e usa-os')
    a = ap.parse_args(argv)

    Lc_h = a.lc if a.lc is not None else lc_h_escolhido()
    av = aeronave(Lc_h)
    if a.gera_dev:
        a.fwd = escreve_avl(av, 'resultados/_tmp/rot_fwd.avl', cg='fwd')
        a.aft = escreve_avl(av, 'resultados/_tmp/rot_aft.avl', cg='aft')
    arqs = {'fwd': a.fwd.replace('\\', '/'), 'aft': a.aft.replace('\\', '/')}
    saida = os.path.join(AQUI, a.saida)
    os.makedirs(saida, exist_ok=True)
    print(f'Lc_h = {Lc_h:.4f}; arquivos: {arqs}; saida: {os.path.relpath(saida, AQUI)}')

    geo = item1(arqs, saida)
    item2(av, saida)
    its = item3(arqs, av, geo, saida)
    estol_res, notas4 = item4(arqs, av, its, saida)
    _, _, notas57 = itens567(arqs, av, its, estol_res, saida)
    l8 = item8(arqs, av, its, saida)
    notas8 = ['Item 8 -- ponto neutro e margem estatica (MS = (xnp - xcg)/cref):']
    for l in l8:
        notas8.append(f'  {l["fonte"]:10s} M {l["M"]:.2f}  CG {NOME_CG[l["cg"]]:9s}  xnp = {l["xnp_m"]:.3f} m '
                      f'({l["xnp_pctCMA"]:.1f}% CMA)  MS = {l["MS"]:+.3f}')

    notas = [
        f'Lab 04 -- roteiro (itens 1 a 8). Lc_h = {Lc_h:.4f}. Arquivos: {arqs["fwd"]}, {arqs["aft"]}.',
        '',
        'Escolhas:',
        f'  - Polares, CL x alfa e CL x de no Mach de cruzeiro (M {av["M"]:.2f}, h {av["h"]:.0f} m).',
        '  - CLmax pelo metodo da secao critica em M 0,2 (baixa velocidade), clmax dos perfis do '
        'Lab 03; cada polar vai de CL = -0,5 ate o CLmax do caso correspondente do item 4 '
        '(sem trimagem -> CLmax sem trimagem; compensada -> CLmax compensado). Em M 0,85 o alfa '
        'desse CL e menor que o alfa_max de M 0,2.',
        f'  - it de cada CG = o do item 3 (fwd {its["fwd"]:.3f} graus, aft {its["aft"]:.3f} graus), '
        'mantido nos itens 4 a 8.',
        f'  - Passo de alfa das varreduras: {PASSO_ALFA} grau.',
        '',
    ] + notas4 + [''] + notas57 + [''] + notas8
    with open(os.path.join(saida, 'notas_roteiro.txt'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(notas) + '\n')
    print('\n' + '\n'.join(notas))


if __name__ == '__main__':
    main(sys.argv[1:])
