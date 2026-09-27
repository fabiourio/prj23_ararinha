'''
Etapa 1 -- varredura da posicao longitudinal da EH (Lc_h), asa sem torcao.

Para cada Lc_h: designTool (W0, S_h, CL de projeto, CGs, dε/dα) -> .avl
fwd/aft -> it que zera o profundor no cruzeiro, CD compensado, margem
estatica -> estol (secao critica, compensado) -> sombreamento e leme.
Escolhe o MAIOR Lc_h viavel (fuselagem e profundor). O sombreamento e
informativo aqui: e exigido na etapa 2, com o alfa de estol da asa ja
torcida (a asa sem torcao estola na ponta cedo demais).

Margem estatica (MS_aft): NAO entra na escolha de Lc_h. Cht = eta_h*S_h e
mantido fixo pelo designTool (dimensionamento do pouso), entao S_h cai com
Lc_h de forma que o produto (braco x area) e a contribuicao da EH ao NP
praticamente nao mudam com Lc_h -- MS_aft nao discrimina entre pontos da
grade (fica em -2,1% a -2,9% no AVL em M 0,85, para qualquer Lc_h testado).
A causa da MS negativa (instavel) foi investigada e e reportada como
achado do item 8 do relatorio (ver ponto_neutro.py): o NP do AVL (36,1%
CMA em M 0,85; 43,2% em M 0,2) fica bem a re do NP do designTool (46,3%
CMA) por causa (a) do momento de Munk da fuselagem tipo corpo esbelto no
AVL, 2-3x o valor empirico de Raymer/Gilruth usado no designTool, e
amplificado pelo fator de Prandtl-Glauert em M 0,85, e (b) do anel da
nacele (-6/-7% CMA), nao modelado no designTool; a contribuicao da EH
concorda entre os dois (+37/+31% AVL vs +38/+34% designTool) e o a.c.
isolado da asa no AVL fica a re de c/4 (34-39% CMA). Por isso MS_aft e
mantida na tabela/figura so como informacao (coluna `ok_ms`, que NAO
entra em `viavel`), junto com MS_aft_M02 (AVL em M 0,2) e MS_aft_dt
(designTool), para comparacao das tres formas de estimar a margem.

Regra de monotonicidade do CD: o AVL tem ruido de malha de ~0,2 count
(ver convergencia_malha.py). CD_aft so e considerado nao-monotono se
algum passo AUMENTAR mais que TOL_CD_COUNT (0,3 count) em relacao ao
ponto anterior; aumentos menores sao ruido de malha e nao derrubam a
escolha pelo maior Lc_h. W0 precisa decrescer estritamente (sem ruido:
vem direto do designTool, nao do AVL).

Rodar de dentro de avl/:   python varredura_eh.py
'''

import csv
import json
import os

import numpy as np

from estilo import plt, AZUL, LARANJA, VERDE, TINTA2
from aeronave import aeronave, lc_h_maximo, FOLGA_FUSELAGEM
from analises import it_para_de_zero
from avl_run import caso
from estol import estol
from gera_avl import escreve_avl, AQUI
from sombra import alfa_saida, fracao_leme_encoberto

LC_MIN, PASSO = 4.0, 0.1
FOLGA_SOMBRA = 2.0      # [graus] exigido na etapa 2; aqui so referencia
DE_MAX = 20.0           # [graus] profundor na compensacao no CLmax, CG dianteiro
MS_MIN = 0.05           # margem estatica minima, CG traseiro
TOL_CD_COUNT = 0.3      # [count] ruido de malha do AVL (~0,2 count); aumentos
                        # menores que isso entre pontos vizinhos nao contam
                        # como nao-monotonicidade de CD_aft
SAIDA = os.path.join(AQUI, 'resultados')


def ponto(lc):
    av = aeronave(lc)
    r = {'Lc_h': lc, 'S_h': av['EH']['S'], 'W0_kgf': av['W0_kgf'], 'CL': av['CL'],
         'folga_fus': av['fuselagem']['L'] - (av['EH']['xr'] + av['EH']['cr'])}
    arq_aft = None
    for cg in ('fwd', 'aft'):
        arq = escreve_avl(av, f'resultados/varredura/{cg}_{lc:.3f}.avl', cg=cg)
        it, c = it_para_de_zero(arq, av['M'], av['CL'])
        e = estol(arq, it=it, trim=True, av=av)
        xcg = av[f'xcg_{cg}']
        # CD = CDp + CDff (Trefftz), a mesma definicao do resto do Lab 04; o
        # CDtot de campo proximo fica so como referencia (CDnf_*)
        r.update({f'it_{cg}': it, f'CD_{cg}': 1e4*(av['CD0'] + c['CDff']),
                  f'CDnf_{cg}': 1e4*c['CD'], f'CDff_{cg}': 1e4*c['CDff'],
                  f'e_{cg}': c['e'], f'MS_{cg}': (c['xnp'] - xcg)/av['Cref'],
                  f'alfa_estol_{cg}': e['alfa'], f'CLmax_{cg}': e['CL'],
                  f'de_estol_{cg}': e['de'], f'eta_estol_{cg}': e['eta_crit']})
        if cg == 'aft':
            arq_aft = arq
    r['alfa_estol'] = min(r['alfa_estol_fwd'], r['alfa_estol_aft'])
    r['alfa_saida'] = alfa_saida(av)
    r['folga_sombra'] = r['alfa_estol'] - r['alfa_saida']
    r['leme_encoberto'] = fracao_leme_encoberto(av)
    r['ok_fus'] = r['folga_fus'] >= FOLGA_FUSELAGEM - 1e-6
    r['ok_sombra'] = r['folga_sombra'] >= FOLGA_SOMBRA
    r['ok_profundor'] = abs(r['de_estol_fwd']) <= DE_MAX
    # MS_aft_M02 (AVL, baixa velocidade) e MS_aft_dt (designTool): so
    # informacao, para comparar com MS_aft (AVL, cruzeiro M 0,85) -- ver
    # docstring do modulo. de=0 (nao compensado): so as derivadas importam.
    m02 = caso(arq_aft, 0.2, alfa=2.0, it=r['it_aft'], de=0.0, derivadas=True)
    r['MS_aft_M02'] = (m02['xnp'] - av['xcg_aft'])/av['Cref']
    r['MS_aft_dt'] = (av['xnp_dt'] - av['xcg_aft'])/av['Cref']
    r['ok_ms'] = r['MS_aft'] >= MS_MIN
    # a sombra so e exigida depois da torcao (etapa 2), e a margem estatica
    # nao discrimina Lc_h (ver docstring): nenhuma das duas entra em viavel.
    r['viavel'] = r['ok_fus'] and r['ok_profundor']
    return r


def figura(linhas, escolhido):
    lc = np.array([l['Lc_h'] for l in linhas])
    paineis = [
        ('W0 [kgf]', [('W0_kgf', AZUL, None)], None),
        ('CD compensado, CDp + CDff (Trefftz) [count]', [('CD_fwd', LARANJA, 'CG dianteiro'),
                                   ('CD_aft', AZUL, 'CG traseiro')], None),
        ('δe no CLmax, CG dianteiro [°]', [('de_estol_fwd', AZUL, None)], -DE_MAX),
        ('margem estática, CG traseiro\n(informativo, não entra em viável)',
         [('MS_aft', AZUL, 'AVL M 0,85'), ('MS_aft_M02', LARANJA, 'AVL M 0,2'),
          ('MS_aft_dt', VERDE, 'designTool')], MS_MIN),
        ('α_estol − α_saída da sombra [°] (asa sem torção)', [('folga_sombra', AZUL, None)], FOLGA_SOMBRA),
        ('fração do leme encoberta', [('leme_encoberto', AZUL, None)], None),
    ]
    fig, eixos = plt.subplots(2, 3, figsize=(10.5, 5.9), sharex=True)
    for ax, (titulo, series, limite) in zip(eixos.flat, paineis):
        for chave, cor, rot in series:
            ax.plot(lc, [l[chave] for l in linhas], '-o', color=cor, label=rot)
        if limite is not None:
            ax.axhline(limite, color=TINTA2, lw=0.8, ls='--')
        for l in linhas:
            if not l['viavel']:
                ax.axvspan(l['Lc_h'] - PASSO/2, l['Lc_h'] + PASSO/2, color='#f1f0ea', lw=0, zorder=0)
        if escolhido is not None:
            ax.axvline(escolhido, color=TINTA2, lw=1.0)
        ax.set_title(titulo, fontsize=9)
        if any(rot for _, _, rot in series):
            ax.legend(fontsize=7)
    for ax in eixos[1]:
        ax.set_xlabel('Lc_h (braço da EH / CMA)')
    fig.suptitle('Varredura da posição da EH (fundo claro: inviável; linha vertical: escolhida; tracejado: limite)')
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.savefig(os.path.join(SAIDA, 'varredura_eh.png'))


def main():
    lc_max = lc_h_maximo()
    grade = list(np.round(np.arange(LC_MIN, lc_max - 1e-9, PASSO), 3)) + [round(lc_max, 4)]
    print(f'Lc_h de {LC_MIN} a {lc_max:.4f} (limite da fuselagem, folga {FOLGA_FUSELAGEM} m)')
    linhas = []
    for lc in grade:
        l = ponto(float(lc))
        linhas.append(l)
        print(f'  Lc_h {lc:6.3f}  W0 {l["W0_kgf"]:9.0f}  CD_aft {l["CD_aft"]:7.2f}'
              f'  de_estol {l["de_estol_fwd"]:6.2f}  MS_aft {l["MS_aft"]:6.3f}'
              f'  folga_sombra {l["folga_sombra"]:5.2f}  leme {l["leme_encoberto"]:.2f}'
              f'  {"VIAVEL" if l["viavel"] else "inviavel"}')

    viaveis = [l for l in linhas if l['viavel']]
    escolhido = max(l['Lc_h'] for l in viaveis) if viaveis else None
    w0 = [l['W0_kgf'] for l in linhas]
    cd = [l['CD_aft'] for l in linhas]
    # W0 precisa decrescer estritamente. CD_aft: so conta como nao-monotono
    # se algum passo AUMENTAR mais que TOL_CD_COUNT (ruido de malha do AVL,
    # ver docstring do modulo).
    w0_monotono = bool(np.all(np.diff(w0) < 0))
    cd_monotono = bool(np.all(np.diff(cd) < TOL_CD_COUNT))
    monotono = w0_monotono and cd_monotono
    if viaveis and not monotono:
        escolhido = min(viaveis, key=lambda l: l['CD_aft'])['Lc_h']
        print(f'W0 ou CD nao caem monotonamente (tolerancia de {TOL_CD_COUNT} count'
              ' de ruido de malha em CD_aft): escolha pelo menor CD viavel')
    print(f'\nLc_h escolhido: {escolhido}   (W0 e CD monotonos: {monotono})')

    os.makedirs(SAIDA, exist_ok=True)
    with open(os.path.join(SAIDA, 'varredura_eh.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)
    with open(os.path.join(SAIDA, 'lc_h_escolhido.json'), 'w') as f:
        json.dump({'Lc_h': escolhido, 'monotono': monotono, 'lc_max': lc_max,
                   'restricoes': {'folga_fus_m': FOLGA_FUSELAGEM,
                                  'folga_sombra_graus_etapa2': FOLGA_SOMBRA,
                                  'de_max_graus': DE_MAX, 'ms_min': MS_MIN,
                                  'tol_cd_count': TOL_CD_COUNT}}, f, indent=2)
    figura(linhas, escolhido)


if __name__ == '__main__':
    main()
