'''
Lab 04 -- secao 3: derivadas de estabilidade e controle (CG traseiro, ponto
de projeto, SEM restricao de trimagem).

Saidas em resultados/derivadas/:
  tab4_dados_gerais.csv   Tabela 4 (designTool / ponto de projeto)
  tab6_alpha0.csv         Tabela 6 (alfa = 0, it = 0, de = 0, M 0,85; ft)
  tab7_ajuste_cd.csv      Tabela 7 (ajuste CD = CD0 + CDa*alfa + CDa2*alfa^2
                          da polar 5(b), alfa em rad) + ajuste_cd_alpha.png
  tab9_derivadas.csv      Tabela 9 completa (valor bruto do AVL, valor final,
                          conversao aplicada e fonte)
  verificacao_q.csv       diferencas finitas em qc/2V (CLq e Cmq conferidos
                          contra o st; CDq, que o AVL nao imprime)
  saida_st_sb.txt         saida bruta do AVL (ft + st + sb) do ponto de
                          projeto e do caso alfa = 0, para rastreabilidade

Convencoes (enunciado):
  - derivadas de controle do AVL (de, da, dr, it) estao em 1/grau: x 180/pi;
  - sinal invertido (convencao do programa de MVO): CYb, CYp, CYr, Cl_da,
    Cl_dr, Cn_da, Cn_dr;
  - Clp, Clr, Cnp, Cnr, Cl_da, Cl_dr, Cn_da, Cn_dr do "sb" (eixos do corpo);
    Clb e Cnb do "st"; o resto do "st".
  - CL_it, CM_it, CD_it: bloco da variavel de projeto (g1 = it); CD_it e
    CD_de sao os de Trefftz (CDffg1, CDffd2), unicos que o AVL imprime.

Escolhas:
  - Ponto de projeto: "a c CL_projeto", "d2 d2 0" e it = incidencia do CG
    traseiro do item 3 (tab_incidencia.csv do roteiro), M 0,85.
  - CDq: o AVL nao imprime; diferenca finita central em q^ = qc/2V = +-0,005
    com alfa FIXO no alfa do ponto de projeto ("p p <valor>" no OPER; o
    comando "q q" nao existe no AVL 3.37). Usa o CD de Trefftz (CDp + CDff),
    a mesma definicao de CD0/CDalfa/CDalfa2 (ajuste da polar) e de CDit/CDde
    (CDffg1/CDffd2): a tabela inteira fica numa so definicao de arrasto. O
    valor de campo proximo (CDvis + CDind) vai para verificacao_q.csv como
    sensibilidade: a razao campo proximo/Trefftz (~2,7x) e a mesma de CDit e
    CDde, ou seja, e diferenca de definicao, nao efeito da rotacao. A mesma
    diferenca finita reproduz o CLq e o Cmq do st (verificacao_q.csv).
  - Tabela 4: momentos de inercia do designTool (moment_of_inertia) com a
    fracao de combustivel do ponto de projeto e 100% de carga paga, em torno
    do CG dessa condicao (x) e da linha z = 0 do designTool (eixos do
    designTool: x para tras, z para cima; Ixz = soma m (x-xcg) z, que tem o
    mesmo sinal nos eixos do corpo x para frente, z para baixo, pois os
    dois eixos trocam de sinal). ip = 0 (o designTool nao tem incidencia de
    motor). xp, zp: centro da nacele (x_n + L_n/2, z_n), coordenadas REAIS
    do designTool com origem no nariz da fuselagem; zp com sinal invertido
    (z para baixo). Tmax = T0 (empuxo estatico ao nivel do mar dos dois
    motores, designTool) x kT(M, h) do modelo de motor do designTool
    (engineTSFC, Howe/Scholz) no Mach e altitude de projeto.

Uso (de dentro de avl/):
  python derivadas_estabilidade.py [--aft aft.avl] [--roteiro resultados/roteiro]
                                   [--saida resultados/derivadas]
'''

import argparse
import csv
import os
import sys

import numpy as np

from estilo import plt, AZUL, LARANJA
from aeronave import aeronave, _analisa
from avl_run import le_resultados
from avl_saida import derivadas, sessao, casos_ft, converte
from gera_avl import AQUI
from lab04_roteiro import lc_h_escolhido, grava_csv, imprime
from designTool.moment_of_inertia import moment_of_inertia
from designTool.propulsion import engineTSFC
from designTool.constants import gravity

DQ = 0.005   # passo em qc/2V da diferenca finita


def tabela4(av, ap):
    inp = ap['inputs']
    moment_of_inertia(ap, fuel_frac=av['fuel_frac'], payload_frac=1.0)
    moi = ap['moment_of_inertia']
    _, kT = engineTSFC(av['M'], av['h'], ap)
    T0 = ap['thrust_matching']['T0']
    xp = inp['x_n'] + inp['L_n']/2
    zp = -inp['z_n']
    m = av['W']/gravity
    obs_moi = (f'moment_of_inertia(fuel_frac = {av["fuel_frac"]:.4f}, payload_frac = 1,0); '
               'origem no CG da condicao, z = 0 do designTool')
    linhas = [
        ('Sref', 'area de referencia [m2]', av['Sref'], 'designTool', ''),
        ('cref', 'corda de referencia [m]', av['Cref'], 'designTool', 'CMA da asa'),
        ('bref', 'envergadura de referencia [m]', av['Bref'], 'designTool', ''),
        ('m', 'massa da aeronave [kg]', m, 'designTool', f'W/g no ponto de projeto, W = {av["W"]:.1f} N'),
        ('Ixx', 'momento de inercia [kg.m2]', moi['Ixx'], 'designTool', obs_moi),
        ('Iyy', 'momento de inercia [kg.m2]', moi['Iyy'], 'designTool', obs_moi),
        ('Izz', 'momento de inercia [kg.m2]', moi['Izz'], 'designTool', obs_moi),
        ('Ixz', 'momento de inercia [kg.m2]', moi['Ixz'], 'designTool',
         obs_moi + '; Ixz = soma m (x-xcg) z, mesmo sinal nos eixos do corpo (x frente, z baixo)'),
        ('ip', 'incidencia do motor [graus]', 0.0, 'designTool',
         'o designTool nao define incidencia de motor: eixo do motor paralelo ao eixo x'),
        ('xp', 'posicao longitudinal do motor [m]', xp, 'designTool',
         f'centro da nacele x_n + L_n/2 = {inp["x_n"]} + {inp["L_n"]}/2; origem no nariz, x para tras; '
         f'xp - xcg_aft = {xp - av["xcg_aft"]:+.4f} m'),
        ('zp', 'posicao vertical do motor [m] (sinal invertido)', zp, 'designTool',
         f'-z_n, z_n = {inp["z_n"]} m (designTool, z para cima, origem no nariz); '
         'positivo = abaixo da linha z = 0'),
        ('Tmax', 'tracao maxima [N]', T0*kT, 'designTool',
         f'2 motores: T0 = {T0:.0f} N (nivel do mar, estatico) x kT = {kT:.5f} '
         f'(engineTSFC, M {av["M"]}, h {av["h"]:.0f} m)'),
        ('V', 'velocidade de voo [m/s]', av['V'], 'ponto de projeto', ''),
        ('h', 'altitude de voo [m]', av['h'], 'ponto de projeto', ''),
    ]
    return [{'parametro': p, 'explicacao': e, 'valor': v, 'fonte': f, 'observacao': o}
            for p, e, v, f, o in linhas]


def tabela6(arquivo, av):
    s = sessao(arquivo, av['M'], 0.0, ['a a 0', 'd2 d2 0', 'x', 'ft', ''])
    r = le_resultados(s)
    linhas = [{'parametro': 'CL0', 'explicacao': 'CL para alfa = 0', 'valor': r['CL'],
               'fonte': 'AVL ft'},
              {'parametro': 'CM0', 'explicacao': 'CM para alfa = 0', 'valor': r['Cm'],
               'fonte': 'AVL ft'}]
    for l in linhas:
        l['condicao'] = (f'alfa = 0, it = 0, de = 0, beta = 0, p = r = 0, da = dr = 0, '
                         f'M {av["M"]}, Xref = {av["xcg_aft"]:.4f} m (CG traseiro)')
    return linhas, s


def tabela7(roteiro, saida):
    with open(os.path.join(roteiro, 'polares.csv'), encoding='utf-8') as f:
        rows = [r for r in csv.DictReader(f) if r['caso'] == 'b']
    alfa = np.radians([float(r['alfa_graus']) for r in rows])
    cd = np.array([float(r['CD_trefftz']) for r in rows])
    c2, c1, c0 = np.polyfit(alfa, cd, 2)
    res = cd - (c0 + c1*alfa + c2*alfa**2)
    rms = float(np.sqrt(np.mean(res**2)))

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    a_f = np.linspace(alfa.min(), alfa.max(), 200)
    ax.plot(np.degrees(alfa), cd, 'o', color=LARANJA, ms=4, label='AVL, polar 5(b)')
    ax.plot(np.degrees(a_f), c0 + c1*a_f + c2*a_f**2, color=AZUL, lw=1.4,
            label=f'ajuste: CD = {c0:.5f} {c1:+.5f}·α {c2:+.4f}·α²  (α em rad)')
    ax.set_xlabel('α [°]')
    ax.set_ylabel('CD = CDp + CDi (Trefftz)')
    ax.set_title(f'Ajuste quadrático da polar não trimada, CG traseiro, M 0,85 '
                 f'(RMS {1e4*rms:.2f} counts)', fontsize=9)
    ax.legend(fontsize=7.5, loc='upper center')
    fig.tight_layout()
    fig.savefig(os.path.join(saida, 'ajuste_cd_alpha.png'))
    plt.close(fig)

    obs = (f'{len(rows)} pontos da polar 5(b) (CG traseiro, de = 0, M 0,85), '
           f'alfa de {np.degrees(alfa.min()):.2f} a {np.degrees(alfa.max()):.2f} graus; '
           f'RMS do residuo {1e4*rms:.3f} counts')
    return [{'parametro': 'CD0', 'explicacao': 'termo constante da polar', 'valor': c0,
             'fonte': 'exercicio 5.b', 'observacao': obs},
            {'parametro': 'CDa', 'explicacao': 'termo linear da polar [1/rad]', 'valor': c1,
             'fonte': 'exercicio 5.b', 'observacao': obs},
            {'parametro': 'CDa2', 'explicacao': 'termo quadratico da polar [1/rad2]', 'valor': c2,
             'fonte': 'exercicio 5.b', 'observacao': obs}]


def ponto_projeto(arquivo, av, it):
    s = sessao(arquivo, av['M'], it, [f'a c {av["CL"]}', 'd2 d2 0', 'x', 'ft', '', 'st', '',
                                       'sb', ''])
    r = le_resultados(s)
    if abs(r['CL'] - av['CL']) > 1e-4:
        raise RuntimeError(f'CL do ponto de projeto nao convergiu: {r["CL"]}')
    return r, derivadas(s, 'st'), derivadas(s, 'sb'), s


def dif_finita_q(arquivo, av, it, alfa, dq=DQ):
    '''Diferenca finita central em q^ = qc/2V com alfa fixo.'''
    s = sessao(arquivo, av['M'], it, ['d2 d2 0', f'a a {alfa:.6f}', f'p p {dq}', 'x',
                                       f'p p {-dq}', 'x'])
    rp, rm = casos_ft(s)
    if abs(rp['alfa'] - alfa) > 1e-4 or abs(rm['alfa'] - alfa) > 1e-4:
        raise RuntimeError('alfa nao ficou fixo na diferenca finita em q')
    d = lambda k: (rp[k] - rm[k])/(2*dq)  # noqa: E731
    return {'CLq': d('CL'), 'Cmq': d('Cm'), 'CDq': d('CDind'), 'CDq_trefftz': d('CDff'),
            'CLq_trefftz': d('CLff')}, s


def tabela9(t4, t6, t7, st, sb, fd):
    linhas = []

    def add(par, expl, bruto, final, conv, fonte):
        # valor_4as: 4 algarismos significativos (o AVL imprime alguns com so 3,
        # p. ex. CDffd2 = 0.000117), e o que vai para o relatorio
        try:
            v4 = float(f'{float(final):.4g}')
        except (TypeError, ValueError):
            v4 = final
        linhas.append({'parametro': par, 'explicacao': expl, 'valor_avl_bruto': bruto,
                       'valor_final': final, 'valor_4as': v4, 'conversao': conv,
                       'fonte': fonte})

    for l in t4:
        add(l['parametro'], l['explicacao'], '', l['valor'],
            'sinal invertido' if l['parametro'] == 'zp' else '', 'Tabela 4 (' + l['fonte'] + ')')
    v6 = {l['parametro']: l['valor'] for l in t6}
    v7 = {l['parametro']: l['valor'] for l in t7}

    def av_st(par, expl, nome, inverte=False, por_grau=False, bloco='st'):
        conv = ', '.join(c for c, on in (('sinal invertido', inverte), ('x180/pi', por_grau)) if on)
        bruto = (st if bloco == 'st' else sb)[nome]
        add(par, expl, bruto, converte(bruto, inverte, por_grau), conv, f'AVL {bloco} ({nome})')

    add('CL0', 'CL para alfa = 0', v6['CL0'], v6['CL0'], '', 'Tabela 6 (AVL ft)')
    av_st('CLa', 'dCL/dalfa [1/rad]', 'CLa')
    av_st('CLq', 'dCL/dq [1/rad]', 'CLq')
    av_st('CLit', 'dCL/dit [1/rad]', 'CLg1', por_grau=True)
    av_st('CLde', 'dCL/dde [1/rad]', 'CLd2', por_grau=True)
    add('CD0', 'CD para alfa = 0 (ajuste)', '', v7['CD0'], '', 'Tabela 7 (exercicio 5.b)')
    add('CDa', 'termo linear da polar [1/rad]', '', v7['CDa'], '', 'Tabela 7 (exercicio 5.b)')
    add('CDa2', 'termo quadratico da polar [1/rad2]', '', v7['CDa2'], '', 'Tabela 7 (exercicio 5.b)')
    add('CDq', 'dCD/dq [1/rad]', '', fd['CDq_trefftz'],
        f'diferenca finita central em qc/2V = +-{DQ} (alfa fixo), CDp + CDff (Trefftz), '
        'consistente com CD0/CDalfa/CDit/CDde',
        'AVL ft ("p p <valor>")')
    av_st('CDit', 'dCD/dit [1/rad]', 'CDffg1', por_grau=True)
    av_st('CDde', 'dCD/dde [1/rad]', 'CDffd2', por_grau=True)
    add('CM0', 'CM para alfa = 0', v6['CM0'], v6['CM0'], '', 'Tabela 6 (AVL ft)')
    av_st('CMa', 'dCM/dalfa [1/rad]', 'Cma')
    av_st('CMq', 'dCM/dq [1/rad]', 'Cmq')
    av_st('CMit', 'dCM/dit [1/rad]', 'Cmg1', por_grau=True)
    av_st('CMde', 'dCM/dde [1/rad]', 'Cmd2', por_grau=True)
    av_st('CYb', 'dCY/dbeta [1/rad]', 'CYb', inverte=True)
    av_st('CYp', 'dCY/dp [1/rad]', 'CYp', inverte=True)
    av_st('CYr', 'dCY/dr [1/rad]', 'CYr', inverte=True)
    av_st('CYdr', 'dCY/ddr [1/rad]', 'CYd3', por_grau=True)
    av_st('Clb', 'dCl/dbeta [1/rad]', 'Clb')
    av_st('Clp', 'dCl/dp [1/rad]', 'Clp', bloco='sb')
    av_st('Clr', 'dCl/dr [1/rad]', 'Clr', bloco='sb')
    av_st('Clda', 'dCl/dda [1/rad]', 'Cld1', inverte=True, por_grau=True, bloco='sb')
    av_st('Cldr', 'dCl/ddr [1/rad]', 'Cld3', inverte=True, por_grau=True, bloco='sb')
    av_st('Cnb', 'dCn/dbeta [1/rad]', 'Cnb')
    av_st('Cnp', 'dCn/dp [1/rad]', 'Cnp', bloco='sb')
    av_st('Cnr', 'dCn/dr [1/rad]', 'Cnr', bloco='sb')
    av_st('Cnda', 'dCn/dda [1/rad]', 'Cnd1', inverte=True, por_grau=True, bloco='sb')
    av_st('Cndr', 'dCn/ddr [1/rad]', 'Cnd3', inverte=True, por_grau=True, bloco='sb')
    return linhas


def it_do_roteiro(roteiro):
    with open(os.path.join(roteiro, 'tab_incidencia.csv'), encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['cg'] == 'aft':
                return float(r['it_graus']), r['arquivo']
    raise RuntimeError('tab_incidencia.csv sem a linha do CG traseiro')


def main(argv=None):
    p = argparse.ArgumentParser(description='Lab 04 -- derivadas de estabilidade (secao 3)')
    p.add_argument('--aft', default=None, help='default: o arquivo aft do roteiro')
    p.add_argument('--roteiro', default=os.path.join('resultados', 'roteiro'))
    p.add_argument('--saida', default=os.path.join('resultados', 'derivadas'))
    p.add_argument('--lc', type=float, default=None)
    a = p.parse_args(argv)

    roteiro = os.path.join(AQUI, a.roteiro)
    saida = os.path.join(AQUI, a.saida)
    os.makedirs(saida, exist_ok=True)
    Lc_h = a.lc if a.lc is not None else lc_h_escolhido()
    av, ap = aeronave(Lc_h), _analisa(Lc_h)
    it, arq_rot = it_do_roteiro(roteiro)
    arquivo = (a.aft or arq_rot).replace('\\', '/')
    print(f'Lc_h = {Lc_h:.4f}; arquivo {arquivo}; it (CG traseiro, item 3) = {it:.4f} graus')

    t4 = tabela4(av, ap)
    grava_csv(os.path.join(saida, 'tab4_dados_gerais.csv'), t4)
    imprime('Tabela 4', t4)

    t6, s6 = tabela6(arquivo, av)
    grava_csv(os.path.join(saida, 'tab6_alpha0.csv'), t6)
    imprime('Tabela 6', t6)

    t7 = tabela7(roteiro, saida)
    grava_csv(os.path.join(saida, 'tab7_ajuste_cd.csv'), t7)
    imprime('Tabela 7', t7)

    r, st, sb, s9 = ponto_projeto(arquivo, av, it)
    fd, sq = dif_finita_q(arquivo, av, it, r['alfa'])
    ver = [{'grandeza': k, 'diferenca_finita': fd[k], 'AVL_st': st.get(k, float('nan')),
            'erro_rel': (fd[k]/st[k] - 1) if k in st else float('nan')}
           for k in ('CLq', 'Cmq', 'CDq', 'CDq_trefftz', 'CLq_trefftz')]
    for v in ver:
        v['observacao'] = {'CDq': 'sensibilidade (campo proximo: CDvis + CDind; ~2,7x o de '
                                  'Trefftz, mesma razao de CDit e CDde)',
                           'CDq_trefftz': 'adotado na Tabela 9 (CDp + CDff, mesma definicao '
                                          'das demais derivadas de arrasto)',
                           'CLq_trefftz': 'Trefftz nao reproduz o CLq do st com rotacao'}.get(
            v['grandeza'], 'conferencia do metodo contra o st')
    grava_csv(os.path.join(saida, 'verificacao_q.csv'), ver)
    imprime(f'Diferenca finita em qc/2V = +-{DQ} (alfa = {r["alfa"]:.4f} graus fixo)', ver)
    if abs(fd['CLq']/st['CLq'] - 1) > 0.02:
        raise RuntimeError('diferenca finita em q nao reproduz o CLq do st (>2%)')

    t9 = tabela9(t4, t6, t7, st, sb, fd)
    grava_csv(os.path.join(saida, 'tab9_derivadas.csv'), t9)
    imprime(f'Tabela 9 (CG traseiro, CL = {av["CL"]:.4f}, alfa = {r["alfa"]:.4f} graus, '
            f'it = {it:.4f} graus, de = 0, M {av["M"]})', t9)

    with open(os.path.join(saida, 'saida_st_sb.txt'), 'w', encoding='utf-8') as f:
        f.write(f'# Lab 04 -- derivadas de estabilidade. Arquivo {arquivo}, it = {it:.6f} graus\n'
                f'# 1) ponto de projeto: a c {av["CL"]}, d2 d2 0, M {av["M"]} (ft, st, sb)\n')
        f.write(s9)
        f.write('\n\n# 2) alfa = 0, it = 0, de = 0 (Tabela 6)\n')
        f.write(s6)
        f.write(f'\n\n# 3) diferenca finita em qc/2V = +-{DQ}, alfa = {r["alfa"]:.6f} fixo\n')
        f.write(sq)


if __name__ == '__main__':
    main(sys.argv[1:])
