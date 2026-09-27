'''
Sensibilidade da margem estatica traseira (MS_aft) a duas mudancas de
projeto, como recomendacao quantitativa para a proxima iteracao (item 8 do
Lab 04).

Motivacao: no projeto entregue (Lc_h = 4,8237, asa com a torcao otimizada)
o ponto neutro do AVL fica bem a frente do do designTool (36,3% CMA em
M 0,85 e 43,3% em M 0,2, contra 46,5% no designTool), e com o CG traseiro
em 39,2% CMA a margem traseira e -2,9% (AVL, M 0,85), +4,1% (AVL, M 0,2) e
+7,3% (designTool). ponto_neutro.py mostra que a EH concorda entre os
metodos; a diferenca vem da fuselagem (momento de Munk) e do anel da nacele.
Decisao: manter o projeto entregue e dizer quanto seria preciso mudar para
  (i)  MS_aft >= 5% pelo AVL em M 0,2   e
  (ii) MS_aft >= 0  pelo AVL em M 0,85,
e quanto isso custa em W0.

Duas estrategias, cada uma uma varredura 1-D:

A) EH maior (Cht ^), com a geometria escolhida fixa (Lc_h = 4,8237): o
   designTool recalcula S_h, W0 e a faixa de CG. Com Lc_h fixo o c/4 da CMA
   da EH fica no lugar e a corda cresce, entao o bordo de fuga da raiz da EH
   passa um pouco do limite da fuselagem (folga < 0,5 m); a folga e
   reportada na tabela.

B) Asa para tras (Delta xr_w), a passos de 0,25 m:
   - os motores vao junto com a asa (x_n += Delta, estao sob a asa);
   - o trem principal vai junto com a asa (x_mlg += Delta), para manter a
     posicao relativa ao caixao da asa (onde fica a fixacao do trem);
   - as empenagens ficam na mesma posicao absoluta: a EH esta no limite da
     fuselagem, entao para cada Delta recalcula-se Lc_h para manter o bordo
     de fuga da raiz da EH em L_f - 0,5 m (lc_h_maximo com os inputs
     modificados), e Lb_v e ajustado para manter o bordo de ataque da raiz
     da EV no mesmo x absoluto;
   - Cht e Cvt fixos: S_h e S_v crescem porque os bracos encurtam;
   - carga paga e tripulacao ficam onde estao (cabine); o CG do combustivel
     acompanha a asa automaticamente (tanque na asa, balance.py).

Para cada ponto: designTool (W0, S_h, S_v, CGs, xnp, SM) e AVL (modelo
completo com a torcao otimizada de resultados/torcao_otimizada.json, CG
traseiro; it que zera o profundor no cruzeiro; Xnp do `st` em M 0,2 e
M 0,85 no CL de projeto, compensado -- mesmo procedimento do item 8 de
lab04_roteiro.py). Os valores exigidos de Cht e Delta xr_w sao obtidos por
interpolacao linear entre os pontos da varredura.

Saidas em resultados/sensibilidade/: sens_cht.csv, sens_xrw.csv,
resumo.json, sensibilidade_ms.png, sensibilidade_w0.png. Temporarios do
AVL em resultados/_tmp/sens_*.avl.

Rodar de dentro de avl/:   python sensibilidade_estabilidade.py
'''

import csv
import json
import os

import numpy as np

from estilo import plt, AZUL, LARANJA, VERDE, TINTA, TINTA2, CINZA
from aeronave import analisa_mod, aeronave_mod, lc_h_maximo, FOLGA_FUSELAGEM
from analises import it_para_de_zero
from avl_run import caso
from designTool.constants import gravity
from designTool.geometry import geometry
from designTool.standard_airplane import standard_airplane
from gera_avl import escreve_avl, AQUI, ETAS

SAIDA = os.path.join(AQUI, 'resultados', 'sensibilidade')
ARQ_TORCAO = os.path.join(AQUI, 'resultados', 'torcao_otimizada.json')

MACHS = (0.2, 0.85)
MS_ALVO = {0.2: 0.05, 0.85: 0.0}      # criterios (i) e (ii)
W0_ENTREGUE_KGF = 288764.0             # W0 do projeto entregue (Lc_h = 4,8237)
# Valores entregues (item 8, lab04_roteiro.py / resultados/roteiro) para a
# conferencia do ponto de partida das varreduras
MS_ENTREGUE = {'AVL_M0.85': -0.029, 'AVL_M0.2': 0.041, 'designTool': 0.073}
TOL_CONFERENCIA = 0.002                # 0,2% CMA

CHT0, PASSO_CHT = 0.70, 0.05
PASSO_XRW = 0.25                       # [m]
MAX_PASSOS = 16


# ---------------------------------------------------------------------------
# Geometria das estrategias

def torcao_otimizada():
    with open(ARQ_TORCAO) as f:
        d = json.load(f)
    if list(d['etas']) != list(ETAS):
        raise ValueError('estacoes da torcao otimizada diferem de gera_avl.ETAS')
    return np.concatenate([[0.0], d['twist']]), d


def _geo(overrides):
    ap = standard_airplane('my_airplane')
    ap['inputs'].update(overrides)
    geometry(ap)
    return ap['inputs'], ap['geometry']


def lb_v_para_xr_v(xr_v_alvo, overrides, lo=0.2, hi=0.8, tol=1e-9):
    '''Lb_v que coloca o bordo de ataque da raiz da EV em xr_v_alvo (x
    absoluto), com os demais inputs em `overrides` (bissecao: xr_v cresce
    com Lb_v).'''
    def xr_v(lb):
        return _geo({**overrides, 'Lb_v': lb})[1]['xr_v']
    if not xr_v(lo) <= xr_v_alvo <= xr_v(hi):
        raise ValueError('xr_v alvo fora do intervalo de Lb_v')
    while hi - lo > tol:
        mid = 0.5*(lo + hi)
        if xr_v(mid) < xr_v_alvo:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo + hi)


def overrides_cht(cht):
    return {'Cht': float(cht)}


def overrides_xrw(dx, base_inp, base_geo):
    '''Inputs da estrategia B para um deslocamento dx [m] da asa.'''
    ov = {'xr_w': base_inp['xr_w'] + dx, 'x_n': base_inp['x_n'] + dx,
          'x_mlg': base_inp['x_mlg'] + dx}
    if dx != 0.0:
        ov['Lc_h'] = lc_h_maximo(lo=3.0, overrides=ov)
        ov['Lb_v'] = lb_v_para_xr_v(base_geo['xr_v'], ov)
    return ov


# ---------------------------------------------------------------------------
# Avaliacao de um ponto (designTool + AVL)

def avalia(overrides, tag, torcao, parametro):
    ap = analisa_mod(overrides)
    av = aeronave_mod(overrides)
    inp, geo = ap['inputs'], ap['geometry']
    bal, lg = ap['balance'], ap['landing_gear']
    c, xm = geo['cm_w'], geo['xm_w']

    def pc(x):
        return (x - xm)/c*100

    arq = escreve_avl(av, f'resultados/_tmp/sens_{tag}.avl', cg='aft', torcao=torcao)
    it, _ = it_para_de_zero(arq, av['M'], av['CL'])
    xnp_avl = {}
    for M in MACHS:
        r = caso(arq, M, cl=av['CL'], it=it, trim=True, derivadas=True)
        xnp_avl[M] = r['xnp']

    W0 = ap['thrust_matching']['W0']/gravity
    fim_eh = geo['xr_h'] + geo['cr_h']
    linha = {
        'parametro': parametro,
        'Cht': inp['Cht'], 'Lc_h': inp['Lc_h'], 'Lb_v': inp['Lb_v'],
        'xr_w_m': inp['xr_w'], 'x_n_m': inp['x_n'], 'x_mlg_m': inp['x_mlg'],
        'S_h_m2': geo['S_h'], 'S_v_m2': geo['S_v'],
        'xr_h_m': geo['xr_h'], 'fim_EH_m': fim_eh, 'folga_fus_m': inp['L_f'] - fim_eh,
        'xr_v_m': geo['xr_v'],
        'W0_kgf': W0, 'dW0_kgf': W0 - W0_ENTREGUE_KGF,
        'dW0_pct': (W0/W0_ENTREGUE_KGF - 1)*100,
        'xm_w_m': xm, 'cm_w_m': c,
        'xcg_fwd_m': bal['xcg_fwd'], 'xcg_aft_m': bal['xcg_aft'],
        'xcg_fwd_pctCMA': pc(bal['xcg_fwd']), 'xcg_aft_pctCMA': pc(bal['xcg_aft']),
        'xnp_dt_pctCMA': pc(bal['xnp']),
        'xnp_avl_M0.85_pctCMA': pc(xnp_avl[0.85]),
        'xnp_avl_M0.2_pctCMA': pc(xnp_avl[0.2]),
        'MS_aft_avl_M0.85': (xnp_avl[0.85] - bal['xcg_aft'])/c,
        'MS_aft_avl_M0.2': (xnp_avl[0.2] - bal['xcg_aft'])/c,
        'MS_aft_dt': bal['SM_aft'], 'MS_fwd_dt': bal['SM_fwd'],
        'MS_fwd_avl_M0.85': (xnp_avl[0.85] - bal['xcg_fwd'])/c,
        'MS_fwd_avl_M0.2': (xnp_avl[0.2] - bal['xcg_fwd'])/c,
        'it_graus': it,
        'tank_excess': bal['tank_excess'],
        'alpha_tipback_graus': np.degrees(lg['alpha_tipback']),
        'alpha_tailstrike_graus': np.degrees(lg['alpha_tailstrike']),
        'frac_nlg_fwd': lg['frac_nlg_fwd'], 'frac_nlg_aft': lg['frac_nlg_aft'],
        'phi_overturn_graus': np.degrees(lg['phi_overturn']),
    }
    return linha


def atende(l):
    return (l['MS_aft_avl_M0.2'] >= MS_ALVO[0.2]
            and l['MS_aft_avl_M0.85'] >= MS_ALVO[0.85])


def varre(valores, avalia_valor, nome):
    '''Avalia os valores em ordem ate os dois criterios serem atendidos, mais
    um passo alem.'''
    linhas, extra = [], None
    for k, v in enumerate(valores):
        l = avalia_valor(v)
        linhas.append(l)
        print(f"  {nome} = {v:6.3f}  W0 = {l['W0_kgf']:9.0f} kgf  "
              f"MS_aft AVL M0,85 = {l['MS_aft_avl_M0.85']*100:+5.1f}%  "
              f"AVL M0,2 = {l['MS_aft_avl_M0.2']*100:+5.1f}%  "
              f"dT = {l['MS_aft_dt']*100:+5.1f}%  (it = {l['it_graus']:+.2f})")
        if extra is not None:
            return linhas
        if atende(l):
            extra = k
    raise RuntimeError(f'{nome}: criterios nao atendidos ate {valores[-1]}')


# ---------------------------------------------------------------------------
# Interpolacao dos valores exigidos

def cruzamento(x, y, alvo):
    '''Primeiro x em que y (crescente com x) atinge `alvo`, por interpolacao
    linear. Devolve x[0] se o ponto inicial ja atende.'''
    x, y = np.asarray(x, float), np.asarray(y, float)
    if y[0] >= alvo:
        return float(x[0])
    for i in range(len(x) - 1):
        if y[i] < alvo <= y[i + 1]:
            return float(x[i] + (alvo - y[i])*(x[i + 1] - x[i])/(y[i + 1] - y[i]))
    return float('nan')


def exigidos(linhas, chave_x, nome):
    x = [l[chave_x] for l in linhas]
    out = {}
    for rot, M in (('AVL_M0.2_MS5', 0.2), ('AVL_M0.85_MS0', 0.85)):
        xe = cruzamento(x, [l[f'MS_aft_avl_M{M:g}'] for l in linhas], MS_ALVO[M])
        out[rot] = xe
    out['ambos'] = max(out.values())
    res = {}
    for rot, xe in out.items():
        W0 = float(np.interp(xe, x, [l['W0_kgf'] for l in linhas]))
        d = {nome: xe, 'W0_kgf': W0, 'dW0_kgf': W0 - W0_ENTREGUE_KGF,
             'dW0_pct': (W0/W0_ENTREGUE_KGF - 1)*100}
        for k in ('MS_aft_avl_M0.85', 'MS_aft_avl_M0.2', 'MS_aft_dt', 'S_h_m2', 'S_v_m2',
                  'Lc_h', 'alpha_tipback_graus', 'alpha_tailstrike_graus',
                  'frac_nlg_aft', 'folga_fus_m'):
            d[k] = float(np.interp(xe, x, [l[k] for l in linhas]))
        res[rot] = d
    return res


# ---------------------------------------------------------------------------
# Saidas

def grava_csv(caminho, linhas):
    with open(caminho, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        for l in linhas:
            w.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in l.items()})


def _virg(v, n=2):
    return f'{v:.{n}f}'.replace('.', ',')


SERIES = [('MS_aft_avl_M0.85', AZUL, 'AVL, M 0,85'),
          ('MS_aft_avl_M0.2', LARANJA, 'AVL, M 0,2'),
          ('MS_aft_dt', VERDE, 'designTool')]


def _painel_ms(ax, linhas, chave_x, chave_req, req, rotulo_x, titulo, fmt_x, unid_x):
    x = np.array([l[chave_x] for l in linhas])
    for chave, cor, rot in SERIES:
        y = np.array([l[chave] for l in linhas])*100
        ax.plot(x, y, '-o', color=cor, lw=1.4, ms=3.5, label=rot)
        ax.plot(x[0], y[0], 'o', ms=8, mfc='none', mec=TINTA, mew=1.0, zorder=5)
    ax.axhline(5.0, color=TINTA2, lw=0.9, ls='--', zorder=0)
    ax.axhline(0.0, color=TINTA2, lw=0.9, ls=':', zorder=0)
    xl = ax.get_xlim()
    ax.text(xl[1], 5.0, 'MS = 5%', ha='right', va='bottom', fontsize=7, color=TINTA2)
    ax.text(xl[1], 0.0, 'MS = 0', ha='right', va='bottom', fontsize=7, color=TINTA2)
    ax.set_xlim(xl)
    yl = ax.get_ylim()
    # ponto entregue (circulos pretos); rotulo acima do circulo do designTool
    y0 = linhas[0]['MS_aft_dt']*100
    ax.annotate('projeto entregue', (x[0], y0), xytext=(0, 10), textcoords='offset points',
                fontsize=7, color=TINTA2, ha='left', va='bottom')
    # cruzamentos: vertical ate o eixo x e rotulo abaixo-direita (as curvas
    # sobem para a direita, esse canto fica livre)
    for rot, cor, alvo in (('AVL_M0.85_MS0', AZUL, 0.0), ('AVL_M0.2_MS5', LARANJA, 5.0)):
        xe = req[rot][chave_req]
        ax.plot([xe, xe], [yl[0], alvo], color=cor, lw=0.8, ls='--', zorder=1)
        ax.plot(xe, alvo, 'D', color=cor, ms=5, mec='white', mew=0.8, zorder=6)
        ax.annotate(f'{fmt_x(xe)}{unid_x}\n$\\Delta W_0$ = {req[rot]["dW0_pct"]:+.2f}%'
                    .replace('.', ','),
                    (xe, alvo), xytext=(6, -5), textcoords='offset points',
                    fontsize=7, color=TINTA, ha='left', va='top')
    ax.set_ylim(yl)
    ax.set_xlabel(rotulo_x)
    ax.set_ylabel('margem estatica traseira MS_aft [% CMA]')
    ax.set_title(titulo, fontsize=9)
    ax.legend(fontsize=7, loc='upper left')


def figura_ms(lin_a, lin_b, req_a, req_b):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
    _painel_ms(a1, lin_a, 'Cht', 'Cht', req_a, 'coeficiente de volume da EH, Cht [-]',
               'A: EH maior (Lc_h = 4,8237 fixo)', lambda v: 'Cht = ' + _virg(v, 3), '')
    _painel_ms(a2, lin_b, 'parametro', 'dxr_w_m', req_b, 'deslocamento da asa para tras, Δxr_w [m]',
               'B: asa para tras (empenagens fixas, Cht e Cvt fixos)',
               lambda v: 'Δxr_w = ' + _virg(v, 2), ' m')
    fig.tight_layout()
    fig.savefig(os.path.join(SAIDA, 'sensibilidade_ms.png'))
    plt.close(fig)


def _painel_w0(ax, linhas, chave_x, chave_req, req, rotulo_x, titulo):
    x = np.array([l[chave_x] for l in linhas])
    W0 = np.array([l['W0_kgf'] for l in linhas])/1000
    ax.plot(x, W0, '-o', color=TINTA2, lw=1.4, ms=3.5, label='W0 (designTool)')
    ax.plot(x[0], W0[0], 'o', ms=8, mfc='none', mec=TINTA, mew=1.0, zorder=5,
            label='projeto entregue')
    for rot, cor, nome, dy in (('AVL_M0.85_MS0', AZUL, 'MS_aft = 0 (AVL, M 0,85)', -12),
                               ('AVL_M0.2_MS5', LARANJA, 'MS_aft = 5% (AVL, M 0,2)', -12)):
        r = req[rot]
        ax.plot(r[chave_req], r['W0_kgf']/1000, 'D', color=cor, ms=6, mec='white', mew=0.8,
                zorder=6, label=nome)
        ax.annotate(f'+{r["dW0_kgf"]:.0f} kgf ({r["dW0_pct"]:+.2f}%)'.replace('.', ','),
                    (r[chave_req], r['W0_kgf']/1000), xytext=(8, dy),
                    textcoords='offset points', fontsize=7, color=TINTA,
                    ha='left', va='center')
    ax.set_xlabel(rotulo_x)
    ax.set_ylabel('W0 [t]')
    ax.set_title(titulo, fontsize=9)
    ax.legend(fontsize=7, loc='upper left')


def figura_w0(lin_a, lin_b, req_a, req_b):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.0))
    _painel_w0(a1, lin_a, 'Cht', 'Cht', req_a, 'coeficiente de volume da EH, Cht [-]',
               'A: EH maior (Lc_h = 4,8237 fixo)')
    _painel_w0(a2, lin_b, 'parametro', 'dxr_w_m', req_b, 'deslocamento da asa para tras, Δxr_w [m]',
               'B: asa para tras (empenagens fixas, Cht e Cvt fixos)')
    fig.tight_layout()
    fig.savefig(os.path.join(SAIDA, 'sensibilidade_w0.png'))
    plt.close(fig)


# ---------------------------------------------------------------------------

def imprime_tabela(titulo, linhas, chave_x, nome_x):
    print(f'\n{titulo}')
    print(f"{nome_x:>7} {'Lc_h':>7} {'S_h':>6} {'S_v':>6} {'W0[kgf]':>9} {'dW0%':>6} "
          f"{'CGaft%':>7} {'NPdt%':>6} {'NP.85%':>7} {'NP.2%':>6} "
          f"{'MS.85':>6} {'MS.2':>6} {'MSdt':>6} {'folgaEH':>7} {'tipb':>5} {'tails':>5}")
    for l in linhas:
        print(f"{l[chave_x]:7.3f} {l['Lc_h']:7.4f} {l['S_h_m2']:6.2f} {l['S_v_m2']:6.2f} "
              f"{l['W0_kgf']:9.0f} {l['dW0_pct']:+6.2f} {l['xcg_aft_pctCMA']:7.2f} "
              f"{l['xnp_dt_pctCMA']:6.2f} {l['xnp_avl_M0.85_pctCMA']:7.2f} "
              f"{l['xnp_avl_M0.2_pctCMA']:6.2f} {l['MS_aft_avl_M0.85']*100:+6.2f} "
              f"{l['MS_aft_avl_M0.2']*100:+6.2f} {l['MS_aft_dt']*100:+6.2f} "
              f"{l['folga_fus_m']:7.3f} {l['alpha_tipback_graus']:5.1f} "
              f"{l['alpha_tailstrike_graus']:5.1f}")


def conferencia(l0):
    '''O primeiro ponto das duas varreduras e o projeto entregue.'''
    obtido = {'AVL_M0.85': l0['MS_aft_avl_M0.85'], 'AVL_M0.2': l0['MS_aft_avl_M0.2'],
              'designTool': l0['MS_aft_dt']}
    ok = all(abs(obtido[k] - MS_ENTREGUE[k]) <= TOL_CONFERENCIA for k in obtido)
    ok_w0 = abs(l0['W0_kgf'] - W0_ENTREGUE_KGF) < 1.0
    return {'MS_aft_obtido': obtido, 'MS_aft_entregue': MS_ENTREGUE,
            'tolerancia': TOL_CONFERENCIA, 'W0_kgf': l0['W0_kgf'],
            'ok': bool(ok and ok_w0)}


def main():
    os.makedirs(SAIDA, exist_ok=True)
    torcao, _ = torcao_otimizada()
    base_inp, base_geo = _geo({})
    print(f"Projeto entregue: Lc_h = {base_inp['Lc_h']:.4f}, Cht = {base_inp['Cht']:.2f}, "
          f"xr_w = {base_inp['xr_w']:.4f} m, bordo de fuga da raiz da EH em "
          f"{base_geo['xr_h'] + base_geo['cr_h']:.4f} m (L_f - {FOLGA_FUSELAGEM} m = "
          f"{base_inp['L_f'] - FOLGA_FUSELAGEM:.4f} m), xr_v = {base_geo['xr_v']:.4f} m")

    print('\nA) Varredura em Cht (Lc_h fixo):')
    chts = [round(CHT0 + k*PASSO_CHT, 4) for k in range(MAX_PASSOS)]
    lin_a = varre(chts, lambda v: avalia(overrides_cht(v), f'cht_{v:.3f}', torcao, v), 'Cht')

    print('\nB) Varredura em Delta xr_w (asa, motores e trem principal para tras):')
    dxs = [round(k*PASSO_XRW, 4) for k in range(MAX_PASSOS)]
    lin_b = varre(dxs, lambda v: avalia(overrides_xrw(v, base_inp, base_geo),
                                        f'xrw_{v:.2f}', torcao, v), 'dxr_w')

    grava_csv(os.path.join(SAIDA, 'sens_cht.csv'), lin_a)
    grava_csv(os.path.join(SAIDA, 'sens_xrw.csv'), lin_b)
    imprime_tabela('Tabela A (MS em % CMA, angulos em graus)', lin_a, 'Cht', 'Cht')
    imprime_tabela('Tabela B (MS em % CMA, angulos em graus)', lin_b, 'parametro', 'dxr_w')

    req_a = exigidos(lin_a, 'Cht', 'Cht')
    req_b = exigidos(lin_b, 'parametro', 'dxr_w_m')

    conf = {'A': conferencia(lin_a[0]), 'B': conferencia(lin_b[0])}
    print('\nConferencia do ponto de partida (projeto entregue):')
    for s, c in conf.items():
        o = c['MS_aft_obtido']
        print(f"  {s}: MS_aft AVL M0,85 = {o['AVL_M0.85']*100:+.2f}%  AVL M0,2 = "
              f"{o['AVL_M0.2']*100:+.2f}%  designTool = {o['designTool']*100:+.2f}%  "
              f"W0 = {c['W0_kgf']:.1f} kgf  -> {'OK' if c['ok'] else 'DIFERE'}")

    print('\nValores exigidos (interpolacao linear):')
    for s, req, chave in (('A', req_a, 'Cht'), ('B', req_b, 'dxr_w_m')):
        for rot, d in req.items():
            print(f"  {s} {rot:14s} {chave} = {d[chave]:.4f}  W0 = {d['W0_kgf']:.0f} kgf  "
                  f"dW0 = {d['dW0_kgf']:+.0f} kgf ({d['dW0_pct']:+.2f}%)  "
                  f"S_h = {d['S_h_m2']:.2f} m2  S_v = {d['S_v_m2']:.2f} m2  "
                  f"MS_dt = {d['MS_aft_dt']*100:+.1f}%")

    l0 = lin_b[0]
    print('\nTrem de pouso na estrategia B (designTool, landing_gear):')
    print(f"  entregue: tipback = {l0['alpha_tipback_graus']:.1f} graus, "
          f"tailstrike = {l0['alpha_tailstrike_graus']:.1f} graus, "
          f"NLG aft = {l0['frac_nlg_aft']*100:.1f}%, NLG fwd = {l0['frac_nlg_fwd']*100:.1f}%")
    for l in lin_b[1:]:
        print(f"  dxr_w = {l['parametro']:.2f} m: tipback = {l['alpha_tipback_graus']:.1f}, "
              f"tailstrike = {l['alpha_tailstrike_graus']:.1f}, "
              f"tipback - tailstrike = {l['alpha_tipback_graus'] - l['alpha_tailstrike_graus']:+.1f}, "
              f"NLG aft = {l['frac_nlg_aft']*100:.1f}%, NLG fwd = {l['frac_nlg_fwd']*100:.1f}%, "
              f"overturn = {l['phi_overturn_graus']:.1f}")
    trem_ok = all(l['alpha_tipback_graus'] > l['alpha_tailstrike_graus']
                  and l['phi_overturn_graus'] < 63.0 for l in lin_b)
    print(f"  tipback > tailstrike e overturn < 63 graus em todos os pontos: {trem_ok}")

    a0 = lin_a[-1]
    print(f"\nA: com Lc_h fixo a corda da EH cresce e a folga ate o fim da fuselagem cai "
          f"de {lin_a[0]['folga_fus_m']:.3f} m para {a0['folga_fus_m']:.3f} m "
          f"(Cht = {a0['Cht']:.2f}).")

    resumo = {
        'descricao': ('Cht e Delta xr_w exigidos para MS_aft >= 5% (AVL, M 0,2) e '
                      'MS_aft >= 0 (AVL, M 0,85), com a penalidade de W0 em relacao '
                      'ao projeto entregue. Interpolacao linear na varredura.'),
        'W0_entregue_kgf': W0_ENTREGUE_KGF,
        'criterios': {'AVL_M0.2_MS5': 'MS_aft >= 0,05 (AVL, M 0,2)',
                      'AVL_M0.85_MS0': 'MS_aft >= 0 (AVL, M 0,85)',
                      'ambos': 'os dois criterios'},
        'A_Cht': {'hipoteses': 'Lc_h = 4,8237 fixo; designTool recalcula S_h, W0 e CGs',
                  'exigido': req_a},
        'B_xr_w': {'hipoteses': ('asa, motores (x_n) e trem principal (x_mlg) deslocados '
                                 'juntos; bordo de fuga da raiz da EH em L_f - 0,5 m '
                                 '(Lc_h recalculado); bordo de ataque da raiz da EV no '
                                 'mesmo x (Lb_v recalculado); Cht e Cvt fixos; carga paga '
                                 'e tripulacao fixas'),
                   'exigido': req_b,
                   'trem_de_pouso_ok': bool(trem_ok)},
        'conferencia_ponto_entregue': conf,
    }
    with open(os.path.join(SAIDA, 'resumo.json'), 'w', encoding='utf-8') as f:
        json.dump(resumo, f, indent=2, ensure_ascii=False)

    figura_ms(lin_a, lin_b, req_a, req_b)
    figura_w0(lin_a, lin_b, req_a, req_b)


if __name__ == '__main__':
    main()
