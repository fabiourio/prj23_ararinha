'''
Item 2 da Atividade 06: AVL em V_A com o fator de carga limite.

Para cada peso (inicio e fim do cruzeiro), CG dianteiro (fwd.avl), nivel do
mar:
  - V_A = V_S1*sqrt(n_max) [EAS = TAS ao nivel do mar];
  - CL total = n W/(q S) (= CLmax, por definicao de V_A);
  - compensado pelo profundor (d2 -> Cm = 0), it do Lab 04 para o CG dianteiro;
  - taxa de arfagem q = g (n - 1)/V_TAS, passada ao AVL como qc/2V.

Saidas em resultados/02_avl/: simulacoes.json, tab_simulacoes.csv,
faixas_<condicao>.csv, saida_avl_<condicao>.txt e cl_cm_envergadura.png.
'''

import os
import re

import matplotlib.pyplot as plt
import numpy as np

import comum as cm
import vn
from avl_run import le_faixas, le_resultados, roda
from estilo import AZUL, CINZA, LARANJA, TINTA2

CORES = {'inicio': AZUL, 'fim': LARANJA}
_NUM = r'-?\d+\.\d+'


def superficies(s):
    '''
    CL e Cm (referidos a Sref, Cref, Xref) de cada superficie, lidos do
    cabecalho de cada superficie no "fs" (5 casas; o "fn" so tem 4), e da
    fuselagem ("fb").
    '''
    i = s.rfind('Surface and Strip Forces by surface')
    CL, Cm = {}, {}
    for m in re.finditer(r'Surface #\s*(\d+)\s+(\S[^\n]*?)\s*\n.*?CLsurf\s*=\s*(\S+).*?Cmsurf\s*=\s*(\S+)',
                         s[i:], flags=re.S):
        CL[m.group(2)] = float(m.group(3))
        Cm[m.group(2)] = float(m.group(4))
    k = s.rfind('Body Forces')
    m = re.search(rf'^\s*1\s+{_NUM}\s+{_NUM}\s+{_NUM}\s+({_NUM})\s+{_NUM}\s+({_NUM}).*Fuselage',
                  s[k:], flags=re.M)
    CL['Fuselage'], Cm['Fuselage'] = float(m.group(1)), float(m.group(2))
    return CL, Cm


def roda_va(M, CL, qhat):
    cmds = ['load ' + cm.ARQ_AVL, 'oper', 'm', f'mn {M:.5f}', '',
            'de', f'1 {cm.IT_FWD}', '', 'd2 pm 0', f'a c {CL:.6f}', f'p p {qhat:.6f}',
            'x', 'fn', '', 'fb', '', 'fs', '', '', 'quit']
    s = roda(cmds)
    if 'Cannot trim' in s or '*******' in s:
        raise RuntimeError('AVL sem solucao valida (Cannot trim ou overflow):\n' + s[-2500:])
    r = le_resultados(s)
    if abs(r['CL'] - CL) > 1e-4 or abs(r['Cm']) > 1e-4:
        raise RuntimeError(f'AVL nao convergiu: CL = {r["CL"]}, Cm = {r["Cm"]}')
    qc2v = float(re.findall(r'qc/2V\s*=\s*([-\d.Ee+]+)', s)[-1])
    if abs(qc2v - qhat) > 1e-5:
        raise RuntimeError(f'qc/2V aplicado ({qc2v}) difere do pedido ({qhat})')
    r['superficies'], r['Cm_superficies'] = superficies(s)
    r['faixas'] = le_faixas(s, superficie=1)
    r['saida'] = s
    return r


def simula(k, cond, av):
    W = cond['W']
    n = vn.n_max_manobra(W)
    V = vn.v_estol(W, av['S'], cm.CLMAX)*np.sqrt(n)     # V_A [EAS = TAS no nivel do mar]
    qdin = 0.5*cm.RHO0*V**2
    CL = n*W/(qdin*av['S'])
    M = V/cm.A0
    q = cm.GRAVIDADE*(n - 1)/V
    qhat = q*av['c_ref']/(2*V)
    r = roda_va(M, CL, qhat)

    sup = r['superficies']
    L = {nome: CLs*qdin*av['S'] for nome, CLs in sup.items()}
    L_asa = L['Wing'] + L['Wing (YDUP)']
    L_winglet = L['Winglet'] + L['Winglet (YDUP)']
    L_emp = L['Horizontal tail'] + L['Horizontal tail (YDUP)'] + L['Vertical tail']
    L_nac = L['Nacelle'] + L['Nacelle (YDUP)']
    L_fus = L['Fuselage']
    # ponto de aplicacao da sustentacao da nacele direita: Cm*q*S*c = L*(Xref - x)
    CLn, Cmn = sup['Nacelle'], r['Cm_superficies']['Nacelle']
    x_nac = av['xcg'] - Cmn*av['c_ref']/CLn
    with open(os.path.join(cm.DIR_AVL, f'saida_avl_{k}.txt'), 'w', encoding='utf-8') as f:
        f.write(r['saida'])
    fx = r['faixas']
    cm.salva_csv([{c: float(fx[c][i]) for c in ('y', 'corda', 'area', 'cl', 'cl_norm', 'cm_c4')}
                  for i in range(len(fx['y']))], os.path.join(cm.DIR_AVL, f'faixas_{k}.csv'))
    return {
        'nome': cond['nome'], 'W': W, 'W_comb': cond['W_comb'], 'n': n,
        'h': 0.0, 'rho': cm.RHO0, 'V_EAS': V, 'V_TAS': V, 'M': M, 'q_din': qdin,
        'CL_alvo': CL, 'CL': r['CL'], 'alfa': r['alfa'], 'de': r['de'], 'Cm': r['Cm'],
        'q_rad_s': q, 'q_graus_s': np.degrees(q), 'qc_2V': qhat,
        'CL_superficies': sup, 'Cm_superficies': r['Cm_superficies'],
        'L_nacele_direita': L['Nacelle'], 'x_nacele_cp': x_nac,
        'L_winglet_direito': L['Winglet'],
        'L_asa': L_asa, 'L_winglet': L_winglet, 'L_empenagem': L_emp,
        'L_EH': L['Horizontal tail'] + L['Horizontal tail (YDUP)'],
        'L_naceles': L_nac, 'L_fuselagem': L_fus, 'L_total': r['CL']*qdin*av['S'],
        'nW': n*W,
        'CL_asa_direita': sup['Wing'],
        'faixas': {c: fx[c] for c in ('y', 'corda', 'area', 'cl', 'cl_norm', 'cm_c4')},
    }


def tabela(av, res):
    ini, fim = res['inicio'], res['fim']
    g = cm.GRAVIDADE
    linhas = []

    def add(par, un, a, b, obs=''):
        linhas.append({'parametro': par, 'unidade': un, 'inicio_cruzeiro': a,
                       'fim_cruzeiro': b, 'observacao': obs})

    add('peso da aeronave', 'N', ini['W'], fim['W'])
    add('peso da aeronave', 'kgf', ini['W']/g, fim['W']/g)
    add('peso da asa (estrutura, duas semiasas)', 'N', av['W_asa'], av['W_asa'], 'designTool (Raymer 15.25)')
    add('peso de combustivel', 'N', ini['W_comb'], fim['W_comb'], 'combustivel restante da missao')
    add('peso de combustivel', 'kgf', ini['W_comb']/g, fim['W_comb']/g)
    add('peso de cada trem de pouso principal', 'N', av['W_tdp'], av['W_tdp'],
        'W_mlg do designTool / 2')
    add('peso de cada motor (instalado)', 'N', av['W_motor'], av['W_motor'],
        f'motor seco {av["W_motor_seco"]:.0f} N x 1,3 (designTool)')
    add('altitude', 'm', ini['h'], fim['h'], 'nivel do mar (caso critico)')
    add('fator de carga n', '-', ini['n'], fim['n'], '25.337(b)')
    add('velocidade equivalente V_A', 'm/s', ini['V_EAS'], fim['V_EAS'])
    add('velocidade verdadeira V_A', 'm/s', ini['V_TAS'], fim['V_TAS'], 'igual a EAS no nivel do mar')
    add('Mach', '-', ini['M'], fim['M'])
    add('pressao dinamica', 'Pa', ini['q_din'], fim['q_din'])
    add('CL total da aeronave', '-', ini['CL'], fim['CL'], 'n W/(q S) = CLmax')
    add('angulo de ataque', 'graus', ini['alfa'], fim['alfa'])
    add('profundor de compensacao', 'graus', ini['de'], fim['de'], f'it = {cm.IT_FWD:.3f} graus')
    add('taxa de arfagem q', 'rad/s', ini['q_rad_s'], fim['q_rad_s'], 'g (n - 1)/V_TAS')
    add('taxa de arfagem q', 'graus/s', ini['q_graus_s'], fim['q_graus_s'])
    add('taxa de arfagem adimensional qc/2V', '-', ini['qc_2V'], fim['qc_2V'], 'entrada do AVL (p p)')
    add('sustentacao da asa', 'N', ini['L_asa'], fim['L_asa'], 'superficies Wing + Wing (YDUP)')
    add('sustentacao da empenagem', 'N', ini['L_empenagem'], fim['L_empenagem'], 'EH + EV')
    add('sustentacao das naceles', 'N', ini['L_naceles'], fim['L_naceles'])
    add('sustentacao dos winglets', 'N', ini['L_winglet'], fim['L_winglet'])
    add('sustentacao da fuselagem', 'N', ini['L_fuselagem'], fim['L_fuselagem'])
    add('sustentacao total', 'N', ini['L_total'], fim['L_total'], 'confere com n W')
    add('n W', 'N', ini['nW'], fim['nW'])
    return linhas


def figura(av, res, caminho):
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    yt, Rf = av['asa']['yt'], av['asa']['R_fus']
    for k, r in res.items():
        fx = r['faixas']
        rot = f'{r["nome"]} ($V_A$ = {cm.virgula(r["V_EAS"], 1)} m/s)'
        axs[0].plot(fx['y'], fx['cl'], color=CORES[k], label=rot)
        axs[1].plot(fx['y'], fx['cm_c4'], color=CORES[k], label=rot)
    for ax, nome in zip(axs, ('$c_\\ell$ local', '$c_m$ local (c/4)')):
        ax.axvspan(0, Rf, color=CINZA, alpha=0.15, lw=0)
        ax.set_xlim(0, yt)
        ax.set_xlabel('y [m]')
        ax.set_ylabel(nome)
    axs[0].text(Rf/2, 0.03, 'fuselagem', rotation=90, ha='center', va='bottom',
                transform=axs[0].get_xaxis_transform(), fontsize=7.5, color=TINTA2)
    axs[0].legend(fontsize=7.5, loc='lower center')
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def main():
    av = cm.aeronave()
    res = {k: simula(k, c, av) for k, c in av['condicoes'].items()}
    cm.salva_json(res, os.path.join(cm.DIR_AVL, 'simulacoes.json'))
    cm.salva_csv(tabela(av, res), os.path.join(cm.DIR_AVL, 'tab_simulacoes.csv'))
    figura(av, res, os.path.join(cm.DIR_AVL, 'cl_cm_envergadura.png'))
    for k, r in res.items():
        print(f'{k:6s} V_A = {r["V_EAS"]:.1f} m/s  M = {r["M"]:.3f}  CL = {r["CL"]:.4f}  '
              f'alfa = {r["alfa"]:.2f}  de = {r["de"]:.2f}  qc/2V = {r["qc_2V"]:.5f}  '
              f'L_asa = {r["L_asa"]/1e3:.0f} kN  L_emp = {r["L_empenagem"]/1e3:.1f} kN  '
              f'L_nac = {r["L_naceles"]/1e3:.0f} kN  L_fus = {r["L_fuselagem"]/1e3:.1f} kN  '
              f'soma/nW = {r["L_total"]/r["nW"]:.5f}')


if __name__ == '__main__':
    main()
