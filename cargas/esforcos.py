'''
Item 3 da Atividade 06: diagramas de esforcos da semiasa direita.

Convencoes: y na envergadura projetada (0 no plano de simetria, b/2 na
ponta); x para tras, origem no nariz (designTool/AVL); forcas positivas para
cima; momento fletor positivo flexiona a ponta para cima; momento torcor
positivo cabra (bordo de ataque para cima), em torno do eixo elastico, uma
reta a 32,5% da corda (centro da caixa entre as longarinas de 12% e 53%).

Cargas sobre a semiasa:
  - sustentacao da asa: faixas do AVL, F = q cl Area, com a largura
    projetada de cada faixa = Area cos(diedro)/corda; F normalizada para a
    sustentacao da superficie dada pelo AVL (diferenca < 0,5%);
  - momento aerodinamico de cada faixa em torno do eixo elastico:
    q cm_c/4 c Area + F (x_EA - x_c/4);
  - sustentacao da nacele (anel do AVL), concentrada em y_n e transmitida
    pelo pilone, no ponto x_cp obtido de Cm/CL da superficie; sustentacao
    vertical do winglet, concentrada na ponta;
  - inercia (para baixo, n vezes o peso): estrutura da asa proporcional a
    corda de 0 a b/2, aplicada a 40% da corda; combustivel proporcional a
    corda no tanque (0 a 90% de b/2), a 32,5% da corda (centro do tanque);
    motor instalado e trem principal concentrados.

Integracao da ponta para a raiz:
  V(y) = int_y^{b/2} f dy' + soma das cargas concentradas fora de y
  M(y) = int_y^{b/2} f (y' - y) dy' + soma P_i (y_i - y)   (= int V dy')
  T(y) = int_y^{b/2} t dy' + soma dos torques concentrados fora de y
O peso do motor, como inercia, entra com sinal negativo (alivia V e M).

Momentos em eixos da corrente com a envergadura projetada; para o eixo
elastico enflechado (Lambda_EA), M_EA ~ M/cos(Lambda_EA) e T_EA ~ T cos(Lambda_EA).

Saidas em resultados/03_esforcos/.
'''

import os

import matplotlib.pyplot as plt
import numpy as np

import comum as cm
from estilo import AZUL, CINZA, LARANJA, TINTA, TINTA2, VERDE

N_CELULAS = 6000
CORES = {'inicio': AZUL, 'fim': LARANJA}


# ---------------------------------------------------------------------------
# integracao

def integra(y, f_cel, pontuais):
    '''
    Cortante e fletor nos nos `y` (crescentes, de 0 a b/2) para a carga
    distribuida f_cel [N/m] constante em cada celula e as cargas
    concentradas [(y_i, P_i)] [N].
    '''
    y = np.asarray(y, float)
    dy = np.diff(y)
    ym = 0.5*(y[:-1] + y[1:])
    F = np.asarray(f_cel, float)*dy
    S0 = np.append(np.cumsum(F[::-1])[::-1], 0.0)
    S1 = np.append(np.cumsum((F*ym)[::-1])[::-1], 0.0)
    V = S0.copy()
    M = S1 - y*S0
    for yi, P in pontuais:
        fora = yi > y
        V += P*fora
        M += P*(yi - y)*fora
    return V, M


def integra_torcao(y, t_cel, torques):
    y = np.asarray(y, float)
    Tq = np.asarray(t_cel, float)*np.diff(y)
    T = np.append(np.cumsum(Tq[::-1])[::-1], 0.0)
    for yi, Ti in torques:
        T += Ti*(yi > y)
    return T


def bordas_faixas(area, corda, diedro):
    '''Bordas das faixas do AVL na envergadura projetada.'''
    larg = np.asarray(area)/np.asarray(corda)*np.cos(diedro)
    return np.concatenate([[0.0], np.cumsum(larg)])


# ---------------------------------------------------------------------------
# cargas

def cargas(av, sim):
    '''Distribuicoes nas celulas e cargas concentradas de uma condicao.'''
    asa = av['asa']
    yt = asa['yt']
    e_ea = cm.e_eixo_elastico(av)
    n, q = sim['n'], sim['q_din']

    y = np.linspace(0.0, yt, N_CELULAS + 1)
    ym = 0.5*(y[:-1] + y[1:])
    dy = np.diff(y)
    c = cm.corda(ym, asa)
    x_ea = cm.x_corda(ym, e_ea, asa)

    # sustentacao e momento aerodinamico (faixas do AVL)
    fx = {k: np.asarray(v) for k, v in sim['faixas'].items()}
    bordas = bordas_faixas(fx['area'], fx['corda'], asa['diedro'])
    larg = np.diff(bordas)
    F = q*fx['cl']*fx['area']
    L_sup = sim['CL_asa_direita']*q*av['S']
    k_L = L_sup/F.sum()
    F = k_L*F
    x_c4 = cm.x_corda(fx['y'], 0.25, asa)
    x_ea_f = cm.x_corda(fx['y'], e_ea, asa)
    Mt = q*fx['cm_c4']*fx['corda']*fx['area'] + F*(x_ea_f - x_c4)
    banda = np.clip(np.searchsorted(bordas, ym) - 1, 0, len(larg) - 1)
    L_cel = (F/larg)[banda]
    ta_cel = (Mt/larg)[banda]
    # conserva exatamente os totais das faixas
    L_cel *= F.sum()/np.sum(L_cel*dy)
    ta_cel *= Mt.sum()/np.sum(ta_cel*dy)

    # inercia distribuida (proporcional a corda)
    W_est = av['W_asa']/2
    ws = c*W_est/np.sum(c*dy)
    t = av['tanque']
    no_tanque = (ym >= t['eta0']*yt) & (ym <= t['eta1']*yt)
    W_cmb = sim['W_comb']/2
    wf = np.where(no_tanque, c, 0.0)
    wf = wf*W_cmb/np.sum(wf*dy)
    x_est = cm.x_corda(ym, cm.X_ESTRUTURA_C, asa)
    x_cmb = cm.x_corda(ym, t['x0_c'] + t['dc']/2, asa)
    is_cel = -n*ws
    if_cel = -n*wf
    ts_cel = is_cel*(x_ea - x_est)
    tf_cel = if_cel*(x_ea - x_cmb)

    # cargas concentradas: (nome, y, x, P)
    def x_ea_em(yy):
        return float(cm.x_corda(yy, e_ea, asa))

    pontuais = [
        ('motor (inércia)', av['y_motor'], av['x_motor'], -n*av['W_motor']),
        ('trem principal (inércia)', av['y_tdp'], av['x_tdp'], -n*av['W_tdp']),
        ('nacele (sustentação)', av['y_motor'], sim['x_nacele_cp'], sim['L_nacele_direita']),
        ('winglet (sustentação)', yt, None, sim['L_winglet_direito']),
    ]
    conc = []
    for nome, yy, xx, P in pontuais:
        braco = 0.0 if xx is None else x_ea_em(yy) - xx
        conc.append({'nome': nome, 'y': yy, 'x': xx, 'P': P, 'braco': braco, 'T': P*braco})

    return {'y': y, 'ym': ym, 'dy': dy, 'k_L': k_L, 'e_ea': e_ea,
            'y_faixa': 0.5*(bordas[:-1] + bordas[1:]), 'L_faixa': F/larg, 't_faixa': Mt/larg,
            'L': L_cel, 'inercia_estrutura': is_cel, 'inercia_combustivel': if_cel,
            't_aero': ta_cel, 't_estrutura': ts_cel, 't_combustivel': tf_cel,
            'concentradas': conc,
            'lambda_ea': np.arctan((cm.x_corda(yt, e_ea, asa) - cm.x_corda(0.0, e_ea, asa))/yt)}


def esforcos(cg):
    '''V, M e T totais e por componente nos nos.'''
    y = cg['y']
    comp_f = {'sustentação da asa': cg['L'],
              'inércia da estrutura': cg['inercia_estrutura'],
              'inércia do combustível': cg['inercia_combustivel']}
    comp_t = {'sustentação da asa': cg['t_aero'],
              'inércia da estrutura': cg['t_estrutura'],
              'inércia do combustível': cg['t_combustivel']}
    zero = np.zeros_like(cg['ym'])
    out = {}
    for nome, f in comp_f.items():
        V, M = integra(y, f, [])
        out[nome] = {'V': V, 'M': M, 'T': integra_torcao(y, comp_t[nome], [])}
    for p in cg['concentradas']:
        V, M = integra(y, zero, [(p['y'], p['P'])])
        out[p['nome']] = {'V': V, 'M': M, 'T': integra_torcao(y, zero, [(p['y'], p['T'])])}
    out['total'] = {k: sum(c[k] for c in out.values()) for k in ('V', 'M', 'T')}
    return out


# ---------------------------------------------------------------------------
# tabelas e figuras

def _valor_em(y, arr, yy):
    return float(np.interp(yy, y, arr))


def tab_concentradas(av, sims, cgs):
    linhas = []
    for k, cg in cgs.items():
        n = sims[k]['n']
        for nome, W, frac in (('estrutura da asa (semiasa)', av['W_asa']/2, cm.X_ESTRUTURA_C),
                              ('combustível (semiasa)', sims[k]['W_comb']/2, cg['e_ea'])):
            linhas.append({'condicao': k, 'carga': nome, 'peso_N': W, 'peso_kgf': W/cm.GRAVIDADE,
                           'n': n, 'forca_N': -n*W, 'y_m': 'distribuído (∝ c)',
                           'x': f'{frac*100:.1f}% da corda', 'braco_EA_m': '', 'torque_Nm': ''})
        for p in cg['concentradas']:
            peso = -p['P']/n if 'inércia' in p['nome'] else ''
            linhas.append({'condicao': k, 'carga': p['nome'],
                           'peso_N': peso, 'peso_kgf': (peso/cm.GRAVIDADE if peso != '' else ''),
                           'n': n if peso != '' else '', 'forca_N': p['P'], 'y_m': p['y'],
                           'x': '' if p['x'] is None else p['x'],
                           'braco_EA_m': p['braco'], 'torque_Nm': p['T']})
    return linhas


def tab_raiz(av, cgs, efs):
    linhas = []
    Rf = av['asa']['R_fus']
    for k, ef in efs.items():
        y = cgs[k]['y']
        for estacao, yy in (('plano de simetria (y = 0)', 0.0),
                            (f'lateral da fuselagem (y = {Rf:.2f} m)', Rf)):
            for comp, d in ef.items():
                linhas.append({'condicao': k, 'estacao': estacao, 'componente': comp,
                               'V_kN': _valor_em(y, d['V'], yy)/1e3,
                               'M_kNm': _valor_em(y, d['M'], yy)/1e3,
                               'T_kNm': _valor_em(y, d['T'], yy)/1e3})
    return linhas


def _marcas(ax, av, rotulos=True):
    for yy, nome in ((av['asa']['R_fus'], 'fuselagem'), (av['y_tdp'], 'trem'),
                     (av['y_motor'], 'motor')):
        ax.axvline(yy, color=CINZA, lw=0.7, ls=':')
        if rotulos:
            ax.text(yy, 0.02, nome, transform=ax.get_xaxis_transform(), rotation=90,
                    ha='right', va='bottom', fontsize=7, color=TINTA2)


def _seta(ax, y, P, h, texto, cor):
    '''Carga concentrada desenhada como seta de comprimento h (simbolico).'''
    ax.annotate('', xy=(y, h), xytext=(y, 0), arrowprops=dict(arrowstyle='-|>', color=cor, lw=1.6))
    ax.text(y + 0.4, h, texto, fontsize=7.5, color=TINTA, va='center')


def figura_cargas(av, sim, cg, caminho):
    ym = cg['ym']
    kN = 1e-3
    conc = {p['nome']: p for p in cg['concentradas']}
    fig, axs = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    a, b, c_, d = axs.ravel()

    yf = cg['y_faixa']
    a.plot(yf, cg['L_faixa']*kN, color=AZUL, label='sustentação da asa (faixas do AVL)')
    p = conc['nacele (sustentação)']
    _seta(a, p['y'], p['P'], 0.7*np.max(cg['L_faixa']*kN),
          f'nacele\n+{cm.virgula(p["P"]*kN, 0)} kN\n(concentrada)', LARANJA)
    a.set_ylabel('força [kN/m]')
    a.set_title('Sustentação', fontsize=9)
    a.legend(fontsize=7.5, loc='upper right')

    b.plot(ym, cg['inercia_estrutura']*kN, color=AZUL, label='estrutura da asa')
    esc = 3.0*np.max(np.abs(cg['inercia_estrutura']*kN))
    Pmax = max(abs(conc[k]['P']) for k in ('motor (inércia)', 'trem principal (inércia)'))
    for nome, rot in (('motor (inércia)', 'motor'), ('trem principal (inércia)', 'trem')):
        p = conc[nome]
        _seta(b, p['y'], p['P'], -esc*abs(p['P'])/Pmax,
              f'{rot}\n{cm.virgula(p["P"]*kN, 0)} kN', LARANJA)
    b.set_ylim(-1.25*esc, 0.25*esc)
    b.set_title(f'Inércia: estrutura, motor e trem ($n$ = {cm.virgula(sim["n"], 1)})', fontsize=9)
    b.legend(fontsize=7.5, loc='lower right')

    c_.plot(ym, cg['inercia_combustivel']*kN, color=AZUL)
    c_.set_ylabel('força [kN/m]')
    c_.set_title(f'Inércia do combustível ({cm.virgula(sim["W_comb"]/cm.GRAVIDADE/1e3, 1)} t na aeronave)',
                 fontsize=9)
    lo = min(c_.get_ylim()[0], -1.0)
    c_.set_ylim(lo, 0.15*abs(lo))

    d.plot(yf, cg['t_faixa']*kN, color=AZUL, label='aerodinâmico (faixas do AVL)')
    d.plot(ym, cg['t_estrutura']*kN, color=VERDE, ls='--', label='inércia da estrutura')
    d.plot(ym, cg['t_combustivel']*kN, color=TINTA2, ls=':',
           label='inércia do combustível (no eixo: nulo)')
    linhas = []
    for nome, rot in (('motor (inércia)', 'motor'), ('trem principal (inércia)', 'trem'),
                      ('nacele (sustentação)', 'nacele')):
        p = conc[nome]
        d.plot([p['y']], [0], marker='o', color=LARANJA, ms=5, ls='none')
        sinal = '+' if p['T'] > 0 else '−'
        linhas.append(f'{rot} (y = {cm.virgula(p["y"], 1)} m): {sinal}{cm.virgula(abs(p["T"])*kN, 0)} kN·m')
    d.text(0.98, 0.97, 'concentrados (●):\n' + '\n'.join(linhas), transform=d.transAxes,
           ha='right', va='top', fontsize=7.5, color=TINTA,
           bbox=dict(boxstyle='round,pad=0.4', fc='white', ec=CINZA, lw=0.6))
    d.set_title('Momento torçor distribuído (em torno do eixo elástico, + cabra)', fontsize=9)
    d.set_ylabel('momento [kN·m/m]')
    d.legend(fontsize=7, loc='lower right')

    for ax in axs.ravel():
        ax.axhline(0, color=CINZA, lw=0.8)
        _marcas(ax, av, rotulos=(ax is a))
        ax.set_xlim(0, av['asa']['yt'])
    for ax in (c_, d):
        ax.set_xlabel('y [m]')
    fig.suptitle(f'{sim["nome"]}: $V_A$ = {cm.virgula(sim["V_EAS"], 1)} m/s, $n$ = '
                 f'{cm.virgula(sim["n"], 1)}, nível do mar', fontsize=10)
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def figura_esforcos(av, sims, cgs, efs, caminho):
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.9))
    escalas = (('V', 1e-3, 'força cortante [kN]'), ('M', 1e-6, 'momento fletor [MN·m]'),
               ('T', 1e-6, 'momento torçor [MN·m]'))
    for ax, (k, s, nome) in zip(axs, escalas):
        for cond, ef in efs.items():
            y = cgs[cond]['y']
            v0 = ef['total'][k][0]*s
            casas = 0 if k == 'V' else 2
            ax.plot(y, ef['total'][k]*s, color=CORES[cond],
                    label=f'{sims[cond]["nome"]} (raiz: {cm.virgula(v0, casas)})')
        ax.axhline(0, color=CINZA, lw=0.8)
        _marcas(ax, av, rotulos=(ax is axs[0]))
        ax.set_xlim(0, av['asa']['yt'])
        ax.set_xlabel('y [m]')
        ax.set_ylabel(nome)
        ax.legend(fontsize=7, loc='upper right' if k != 'T' else 'lower right')
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def main():
    av = cm.aeronave()
    sims = cm.le_json(os.path.join(cm.DIR_AVL, 'simulacoes.json'))
    cgs = {k: cargas(av, s) for k, s in sims.items()}
    efs = {k: esforcos(cg) for k, cg in cgs.items()}

    for k, cg in cgs.items():
        figura_cargas(av, sims[k], cg, os.path.join(cm.DIR_ESFORCOS, f'cargas_{k}.png'))
        dist = [{'y_m': cg['ym'][i], 'sustentacao_N_m': cg['L'][i],
                 'inercia_estrutura_N_m': cg['inercia_estrutura'][i],
                 'inercia_combustivel_N_m': cg['inercia_combustivel'][i],
                 't_aero_Nm_m': cg['t_aero'][i], 't_estrutura_Nm_m': cg['t_estrutura'][i],
                 't_combustivel_Nm_m': cg['t_combustivel'][i]}
                for i in range(0, len(cg['ym']), 20)]
        cm.salva_csv(dist, os.path.join(cm.DIR_ESFORCOS, f'distribuicoes_{k}.csv'))
        ef = efs[k]['total']
        cm.salva_csv([{'y_m': cg['y'][i], 'V_N': ef['V'][i], 'M_Nm': ef['M'][i], 'T_Nm': ef['T'][i]}
                      for i in range(0, len(cg['y']), 20)],
                     os.path.join(cm.DIR_ESFORCOS, f'esforcos_{k}.csv'))
    figura_esforcos(av, sims, cgs, efs, os.path.join(cm.DIR_ESFORCOS, 'esforcos_internos.png'))
    cm.salva_csv(tab_concentradas(av, sims, cgs), os.path.join(cm.DIR_ESFORCOS, 'tab_cargas.csv'))
    cm.salva_csv(tab_raiz(av, cgs, efs), os.path.join(cm.DIR_ESFORCOS, 'tab_raiz.csv'))

    for k, cg in cgs.items():
        ef = efs[k]
        print(f'== {sims[k]["nome"]}  (k_L = {cg["k_L"]:.4f}, eixo elastico a '
              f'{cg["e_ea"]*100:.1f}% da corda, Lambda_EA = {np.degrees(cg["lambda_ea"]):.1f} graus)')
        for comp, d in ef.items():
            print(f'   {comp:28s} V0 = {d["V"][0]/1e3:9.1f} kN   M0 = {d["M"][0]/1e6:7.3f} MN.m   '
                  f'T0 = {d["T"][0]/1e6:7.3f} MN.m')


if __name__ == '__main__':
    main()
