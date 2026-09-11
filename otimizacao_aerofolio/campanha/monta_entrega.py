'''
Monta o zip da entrega -- Ararinha_PRJ23_Lab03.zip

O professor re-executa o codigo, entao o zip precisa ser auto-contido e
pronto-para-rodar: material dele, nosso codigo, resultados e o passo a passo.

Ficam de fora os subprodutos das execucoes (solution.vtk, grid.xyz, wall.dat,
derivatives.dat, pickles de historico) e as pastas de rodadas arquivadas,
que existem so para a discussao de convergencia e nao para o professor
re-rodar.

Rodar com:  python monta_entrega.py
'''

import os
import sys
import zipfile

AQUI = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.normpath(os.path.join(AQUI, '..'))
RAIZ = os.path.normpath(os.path.join(LAB, '..'))

NOME = 'Ararinha_PRJ23_Lab03.zip'

# subprodutos e material que nao deve ir
IGNORAR_NOMES = {'solution.vtk', 'grid.xyz', 'grid.png', 'wall.dat',
                 'derivatives.dat', 'historico.pickle', '.gitignore'}
IGNORAR_EXT = {'.pyc', '.zip'}
IGNORAR_PASTAS = {'__pycache__', 'build', 'ftol_1e-6', 'ftol_1e-5_singular',
                  '.egg-info', 'paraview_layouts'}


def deve_ir(caminho):
    partes = caminho.replace('\\', '/').split('/')
    if any(p in IGNORAR_PASTAS or p.endswith('.egg-info') for p in partes):
        return False
    nome = partes[-1]
    if nome in IGNORAR_NOMES:
        return False
    if os.path.splitext(nome)[1].lower() in IGNORAR_EXT:
        return False
    return True


def main():
    destino = os.path.join(RAIZ, NOME)
    n, total = 0, 0
    with zipfile.ZipFile(destino, 'w', zipfile.ZIP_DEFLATED) as z:
        for base, _, arquivos in os.walk(LAB):
            for a in arquivos:
                cheio = os.path.join(base, a)
                rel = os.path.relpath(cheio, RAIZ)
                if not deve_ir(rel):
                    continue
                z.write(cheio, rel)
                n += 1
                total += os.path.getsize(cheio)
        # o documento de definicao do problema vai junto
        doc = os.path.join(RAIZ, 'docs', 'superpowers', 'specs',
                           '2026-09-10-otimizacao-aerofolio-design.md')
        if os.path.isfile(doc):
            z.write(doc, 'otimizacao_aerofolio/DEFINICAO_DO_PROBLEMA.md')
            n += 1

    tam = os.path.getsize(destino) / 1e6
    print(f'{NOME}: {n} arquivos, {total/1e6:.1f} MB crus -> {tam:.1f} MB')
    print(f'gravado em {destino}')

    # confere que o essencial entrou
    with zipfile.ZipFile(destino) as z:
        nomes = z.namelist()
    essenciais = ['campanha/otimiza_secao.py', 'campanha/README.md',
                  'campanha/xfoil_runner.py', 'eulerblock/package',
                  'Lab03_PRJ23_2026.pdf', 'DEFINICAO_DO_PROBLEMA.md']
    print('\nconferencia:')
    for e in essenciais:
        ok = any(e in n for n in nomes)
        print(f'  {"ok " if ok else "FALTA"} {e}')
    figs = [n for n in nomes if n.endswith('.png')]
    print(f'  {len(figs)} figuras')
    return 0


if __name__ == '__main__':
    sys.exit(main())
