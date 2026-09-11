'''
Monta a pasta ENTREGA/ -- so o que serve para o relatorio.

Organizada pelos itens do roteiro, com nomes que dizem o que a figura mostra.
Fica fora: log de execucao, pickle de historico, CSV intermediario de DOE,
rodadas arquivadas, subproduto de solver.

Rodar com:  python monta_pasta_entrega.py
'''

import os
import shutil
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(AQUI, 'resultados')
LAB = os.path.normpath(os.path.join(AQUI, '..'))
DEST = os.path.join(LAB, 'ENTREGA')

EST = ('raiz', 'meio', 'ponta')

# (pasta destino, arquivo de origem, nome final)
PLANO = []

PLANO += [('01_ponto_de_projeto', 'tabelas_relatorio.md',
           'tabelas_1_e_2.md')]

for e in EST:
    PLANO += [
        ('04_convergencia', f'otim_{e}_convergencia.png',
         f'convergencia_{e}.png'),
        ('05_geometria', f'otim_{e}_geometria.png', f'geometria_{e}.png'),
        ('06_cp_e_mach', f'otim_{e}_cp_mach.png', f'cp_mach_{e}.png'),
        ('07_polar_transonica', f'polar_transonica_{e}.png',
         f'polar_transonica_{e}.png'),
        ('07_polar_transonica', f'polar_transonica_{e}.csv',
         f'polar_transonica_{e}.csv'),
        ('08_subsonico_xfoil', f'subsonico_{e}.png', f'subsonico_{e}.png'),
        ('09_batentes', f'item9_{e}.png', f'batente_vs_cusp_{e}.png'),
    ]

PLANO += [
    ('03_restricao_clmax', 'doe_clmax_corte.png',
     'clmax_vs_afiamento_do_nariz.png'),
    ('03_restricao_clmax', 'doe_clmax_dispersao.png',
     'qual_descritor_preve_clmax.png'),
    ('03_restricao_clmax', 'doe_clmax_limiar.png',
     'limiar_da_restricao.png'),
    ('03_restricao_clmax', 'limiares_bluntez.csv',
     'limiares_por_estacao.csv'),
    ('10_verificacao', 'verificacao_malha.csv', 'convergencia_de_malha.csv'),
    ('10_verificacao', 'verificacao_adjunto.csv',
     'adjunto_vs_diferencas_finitas.csv'),
    ('10_verificacao', 'gradiente_clmax.csv',
     'gradiente_do_clmax_e_ruido.csv'),
    ('11_realimentacao_lab02', 'comparacao_korn.csv',
     'korn_previsto_vs_euler_medido.csv'),
]

DOCS = [('RESULTADOS.md', 'RESULTADOS.md'),
        ('README.md', 'COMO_RODAR.md')]


def main():
    shutil.rmtree(DEST, ignore_errors=True)
    copiados, faltando = 0, []

    for pasta, origem, destino in PLANO:
        src = os.path.join(RES, origem)
        if not os.path.isfile(src):
            faltando.append(origem)
            continue
        d = os.path.join(DEST, pasta)
        os.makedirs(d, exist_ok=True)
        shutil.copyfile(src, os.path.join(d, destino))
        copiados += 1

    for origem, destino in DOCS:
        src = os.path.join(AQUI, origem)
        if os.path.isfile(src):
            shutil.copyfile(src, os.path.join(DEST, destino))
            copiados += 1
        else:
            faltando.append(origem)

    doc = os.path.join(os.path.normpath(os.path.join(LAB, '..')),
                       'docs', 'superpowers', 'specs',
                       '2026-09-10-otimizacao-aerofolio-design.md')
    if os.path.isfile(doc):
        shutil.copyfile(doc, os.path.join(DEST, 'DEFINICAO_DO_PROBLEMA.md'))
        copiados += 1

    # coordenadas dos perfis otimizados, no formato do XFoil
    exporta_perfis(os.path.join(DEST, '02_perfis_otimizados'))

    print(f'ENTREGA/ montada: {copiados} arquivos')
    if faltando:
        print(f'\nFALTANDO ({len(faltando)}):')
        for f in faltando:
            print(f'  {f}')
    print()
    for base, _, arqs in sorted(os.walk(DEST)):
        rel = os.path.relpath(base, DEST)
        if arqs:
            print(f'  {rel if rel != "." else "(raiz)"}/')
            for a in sorted(arqs):
                kb = os.path.getsize(os.path.join(base, a)) / 1024
                print(f'      {a:44s} {kb:7.0f} KB')


def exporta_perfis(pasta):
    '''Coordenadas dos perfis otimizados, prontas para carregar no XFoil.'''
    sys.path.insert(0, AQUI)
    import otimiza_secao as ot
    import xfoil_runner as xr
    from analisa_otimos import carrega

    os.makedirs(pasta, exist_ok=True)
    for e in EST:
        for suf, rot in [('', 'roteiro'), ('_cusp', 'cusp_liberado')]:
            d = carrega(f'otim_{e}{suf}')
            if d is None:
                continue
            Al, Au, alpha = ot.desmonta(d['xx_otimo'])
            foil = xr.cst_coords(Au, Al)
            xr._export_quiet(foil, os.path.join(pasta, f'{e}_{rot}.dat'))
            # formato pronto para colar num script do professor
            fmt = lambda v: '[' + ', '.join(f'{x:.8f}' for x in v) + ']'
            i = d['i_otimo']
            with open(os.path.join(pasta, f'{e}_{rot}_cst.txt'), 'w') as fid:
                fid.write(f'# Estacao {e} -- batentes: {rot}\n')
                fid.write(f'# Aeronave B, M_n = {ot.MACH_N}, '
                          f'c_l de projeto = {ot.ESTACOES[e]["cl_ref"]}\n#\n')
                fid.write(f'Al = {fmt(Al)}\n')
                fid.write(f'Au = {fmt(Au)}\n')
                fid.write(f'alpha = {np.degrees(alpha):.6f}  # graus\n#\n')
                fid.write(f'# c_l obtido  = {d["hist"]["CL"][i]:.6f}\n')
                fid.write(f'# c_d obtido  = {d["hist"]["CD"][i]:.8f}\n')
                fid.write(f'# t/c maximo  = {d["hist"]["maxt"][i]:.6f} '
                          f'(minimo exigido {ot.ESTACOES[e]["tc_ref"]})\n')


import numpy as np  # noqa: E402  (usado em exporta_perfis)

if __name__ == '__main__':
    main()
