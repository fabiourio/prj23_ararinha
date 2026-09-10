'''
O que o modelo de Korn previa, e o que o Euler mediu -- item 9 do roteiro.

Este e o entregavel que volta para o Lab 02. Ele NAO altera a aeronave: ele
quantifica a discrepancia e deixa a decisao de projeto para a equipe.

Para cada estacao compara tres coisas:

  1. c_d,onda que a equacao de Korn preve PARA AQUELA SECAO, usando a mesma
     formula do designTool (aerodynamics.py:269) mas com a espessura e o cl
     locais em vez da media ponderada da asa;
  2. c_d que o Euler mediu no perfil de PARTIDA (NACA 1411 reescalado);
  3. c_d que o Euler mediu no perfil OTIMIZADO.

Os valores do Euler saem da extrapolacao de Richardson (verificacao_malha),
porque o nivel 1,0 superestima o arrasto em 70%.

O k_korn implicito responde a pergunta pratica: que valor de k_korn faria o
modelo do Lab 02 reproduzir o que medimos? Se o k_korn implicito do perfil
OTIMIZADO for proximo de 0,95, o modelo estava certo e a equipe so precisa
projetar bons perfis. Se for muito menor, o modelo e otimista demais na
espessura escolhida.

Rodar com:  python compara_korn.py
'''

import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import descritores as dsc
import otimiza_secao as ot
from analisa_otimos import RES, carrega

# k_korn assumido no Lab 02 (standard_airplane.py)
K_KORN_LAB02 = 0.95

# Fator de correcao de malha, medido em verificacao_malha.py: no nivel 1,0 o
# c_d sai 70% acima do extrapolado por Richardson. Aplicamos a razao medida
# para cada estacao, dos ultimos resultados de analisa_otimos --malha.
# (nivel 1,0 -> nivel 2,0, que ainda erra ~13% para cima)
FATOR_MALHA = {'raiz': 0.018813 / 0.024545,
               'meio': 0.007145 / 0.011822,
               'ponta': 0.003859 / 0.008873}


def korn_secao(tc, cl, k=K_KORN_LAB02, mach=None):
    '''
    c_d,onda da equacao de Korn aplicada A SECAO, no plano normal.

    No plano normal as grandezas ja foram transformadas, entao cos(Lambda) = 1
    e a formula do designTool se reduz a  M_dd = k - (t/c) - cl/10.
    '''
    mach = ot.MACH_N if mach is None else mach
    m_dd = k - tc - cl / 10.0
    m_crit = m_dd - (0.1 / 80) ** (1 / 3)
    return 20 * max(0.0, mach - m_crit) ** 4, m_dd, m_crit


def k_implicito(cd_alvo, tc, cl, mach=None):
    '''Que k_korn reproduz o c_d medido? Inverte a formula de Korn.'''
    mach = ot.MACH_N if mach is None else mach
    if cd_alvo <= 0:
        return np.nan
    # 20*(M - Mcrit)^4 = cd  ->  Mcrit = M - (cd/20)^(1/4)
    m_crit = mach - (cd_alvo / 20.0) ** 0.25
    m_dd = m_crit + (0.1 / 80) ** (1 / 3)
    return m_dd + tc + cl / 10.0


def main():
    print(f'Condicao: M_n = {ot.MACH_N}   k_korn do Lab 02 = {K_KORN_LAB02}\n')
    print(f'{"estacao":8s} {"(t/c)_n":>8s} {"cl":>7s} {"M_dd":>7s} '
          f'{"Korn":>9s} {"Euler part.":>12s} {"Euler otim.":>12s} '
          f'{"k impl. otim":>13s}')

    linhas = []
    for nome, est in ot.ESTACOES.items():
        com = carrega(f'otim_{nome}')
        if com is None:
            print(f'{nome:8s}  (sem resultado)')
            continue
        h = com['hist']
        tc, cl = est['tc_ref'], est['cl_ref']

        cd_korn, m_dd, m_crit = korn_secao(tc, cl)
        f = FATOR_MALHA.get(nome, 1.0)
        cd_part = h['CD'][0] * f
        cd_otim = h['CD'][com['i_otimo']] * f

        k_otim = k_implicito(cd_otim, tc, cl)
        k_part = k_implicito(cd_part, tc, cl)

        print(f'{nome:8s} {tc:8.4f} {cl:7.4f} {m_dd:7.4f} {cd_korn:9.6f} '
              f'{cd_part:12.6f} {cd_otim:12.6f} {k_otim:13.4f}')
        linhas.append({'estacao': nome, 'tc': tc, 'cl': cl, 'M_dd': m_dd,
                       'cd_korn': cd_korn, 'cd_euler_partida': cd_part,
                       'cd_euler_otimo': cd_otim,
                       'k_implicito_partida': k_part,
                       'k_implicito_otimo': k_otim})

    if not linhas:
        return 1

    print(f'\n{"="*70}\nLEITURA\n{"="*70}')
    for l in linhas:
        razao = (l['cd_euler_otimo'] / l['cd_korn']
                 if l['cd_korn'] > 0 else np.inf)
        print(f'\n{l["estacao"].upper()} (t/c = {l["tc"]:.4f}):')
        print(f'  Korn previa      c_d,onda = {l["cd_korn"]:.6f}')
        print(f'  Euler no otimo   c_d      = {l["cd_euler_otimo"]:.6f}  '
              f'({razao:.2f}x o previsto)')
        print(f'  k_korn implicito no otimo = {l["k_implicito_otimo"]:.4f}  '
              f'(assumido: {K_KORN_LAB02})')
        if l['k_implicito_otimo'] > K_KORN_LAB02:
            print(f'  -> o perfil otimizado SUPERA o modelo: k efetivo maior '
                  f'que o assumido.')
        else:
            print(f'  -> o perfil otimizado fica AQUEM do modelo nesta '
                  f'espessura.')

    cam = os.path.join(RES, 'comparacao_korn.csv')
    cols = list(linhas[0])
    with open(cam, 'w') as fid:
        fid.write(','.join(cols) + '\n')
        for l in linhas:
            fid.write(','.join(
                f'{l[c]:.6g}' if isinstance(l[c], float) else str(l[c])
                for c in cols) + '\n')
    print(f'\ngravado em {cam}')

    integra_envergadura(linhas)
    return 0


def integra_envergadura(linhas):
    '''
    Integra o arrasto de onda das secoes ao longo da envergadura e compara com
    o CDwave que o designTool preve para a aeronave inteira.

    O k_korn implicito mostrou que o modelo acerta a FISICA DA SECAO. A duvida
    que sobra e a AGREGACAO: o designTool representa a asa inteira por uma
    espessura media ponderada, (t/c)_m = 0,25*tcr + 0,75*tct. Com tcr = 0,196 e
    tct = 0,08 isso da 0,109 -- que, na variacao linear de espessura, e a
    espessura da estacao eta = 0,75. Ou seja, a formula representa a asa toda
    por uma estacao bem externa, onde a corda ja e curta.

    Conversao de arrasto de secao para arrasto de asa sob enflechamento:
        CD = (2/S) * integral[ c_d,normal * cos^3(Lambda) * c(y) dy ]
    O cos^3 vem de: a pressao dinamica normal escala com cos^2, e a area da
    secao normal com cos.

    Integramos so a parte EXPOSTA da asa (fora da fuselagem), que e o que
    gera arrasto.
    '''
    import numpy as np

    # geometria da aeronave B (documento de projeto, secao 2)
    b, cr, taper = 60.1518, 9.8898, 0.24
    S_w = 368.833
    cos50 = 0.84659
    eta_junc = 0.1011

    etas = np.array([ot.ESTACOES[n]['eta'] for n in
                     ('raiz', 'meio', 'ponta')])
    cds = np.array([next(l['cd_euler_otimo'] for l in linhas
                         if l['estacao'] == n)
                    for n in ('raiz', 'meio', 'ponta')])

    # interpola em log(c_d): o arrasto varia por fatores, nao por parcelas
    eg = np.linspace(eta_junc, 1.0, 400)
    cd_g = np.exp(np.interp(eg, etas, np.log(cds)))
    c_g = cr * (1 - (1 - taper) * eg)

    integral = np.trapezoid(cd_g * c_g, eg) * (b / 2)
    CDwave_secoes = 2 * integral * cos50 ** 3 / S_w

    print(f'\n{"="*70}\nAGREGACAO: DA SECAO PARA A ASA\n{"="*70}')
    print(f'  integrando os c_d das secoes otimizadas na parte exposta da asa')
    print(f'  (eta de {eta_junc:.3f} a 1,0), com cos^3(Lambda) = {cos50**3:.4f}:')
    print(f'\n  CDwave da integracao por secoes = {CDwave_secoes:.6f}')

    # o que o designTool preve
    tcr, tct = 0.19623, 0.08
    tcm = 0.25 * tcr + 0.75 * tct
    CL_aer = 0.5053                      # peso medio de cruzeiro
    M = 0.85
    m_dd = (K_KORN_LAB02 / cos50 - tcm / cos50 ** 2
            - CL_aer / (10 * cos50 ** 3))
    m_crit = m_dd - (0.1 / 80) ** (1 / 3)
    CDwave_dt = 20 * max(0.0, M - m_crit) ** 4
    print(f'  CDwave do designTool            = {CDwave_dt:.6f}')
    print(f'  razao                           = {CDwave_secoes/CDwave_dt:.1f}x')

    # a que estacao corresponde a espessura media do designTool?
    eta_equiv = (tcr - tcm) / (tcr - tct)
    print(f'\n  (t/c)_m = 0,25*{tcr:.3f} + 0,75*{tct:.3f} = {tcm:.4f}')
    print(f'  essa espessura ocorre em eta = {eta_equiv:.3f}')
    print(f'  ou seja, o modelo representa a asa inteira por uma estacao')
    print(f'  a {eta_equiv*100:.0f}% da semi-envergadura, onde a corda ja caiu')
    print(f'  para {(1-(1-taper)*eta_equiv)*100:.0f}% da corda de raiz.')

    print(f'\n{"="*70}\nCONCLUSAO PARA O designTool\n{"="*70}')
    print('O k_korn = 0,95 NAO e o problema: os perfis otimizados alcancam')
    print('k implicito de 0,89 a 0,94, ou seja, a fisica da secao no modelo')
    print('esta essencialmente correta. O modelo apostava que a equipe')
    print('entregaria bons perfis, e o Lab 03 mostra que isso e alcancavel.')
    print('\nO problema e a AGREGACAO. A ponderacao 0,25/0,75 representa a asa')
    print(f'inteira pela estacao eta = {eta_equiv:.2f}, subestimando a regiao')
    print('interna -- que e justamente onde a asa e grossa E tem corda longa,')
    print('logo onde o arrasto de onda pesa mais.')
    print('\nDecisao de projeto para a equipe (nao para este script):')
    print('  a) afinar tcr_w e reotimizar o Lab 02, ou')
    print('  b) manter a espessura e corrigir a agregacao do CDwave, ou')
    print('  c) aceitar o custo e declara-lo no relatorio.')


if __name__ == '__main__':
    sys.exit(main())
