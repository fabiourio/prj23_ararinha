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
IGNORAR_EXT = {'.pyc', '.zip',
               # subprodutos de compilacao do LaTeX
               '.aux', '.log', '.out', '.toc', '.synctex.gz', '.fls',
               '.fdb_latexmk'}
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
    # o relatorio fica na raiz do repositorio, fora de LAB, mas as instrucoes
    # do roteiro pedem TUDO num zip so -- entao ele entra junto
    RELATORIO = os.path.join(RAIZ, 'relatorio')
    with zipfile.ZipFile(destino, 'w', zipfile.ZIP_DEFLATED) as z:
        for raiz_busca in (LAB, RELATORIO):
            if not os.path.isdir(raiz_busca):
                continue
            for base, _, arquivos in os.walk(raiz_busca):
                for a in arquivos:
                    cheio = os.path.join(base, a)
                    rel = os.path.relpath(cheio, RAIZ)
                    if not deve_ir(rel):
                        continue
                    z.write(cheio, rel)
                    n += 1
                    total += os.path.getsize(cheio)
        # a definicao do problema ja vive em ENTREGA/ e entra pelo os.walk
        # acima; a copia que vinha de docs/ foi removida do repositorio

    tam = os.path.getsize(destino) / 1e6
    print(f'{NOME}: {n} arquivos, {total/1e6:.1f} MB crus -> {tam:.1f} MB')
    print(f'gravado em {destino}')

    # confere que o essencial entrou
    with zipfile.ZipFile(destino) as z:
        nomes = z.namelist()
    essenciais = ['campanha/otimiza_secao.py',
                  'campanha/xfoil_runner.py', 'eulerblock/package',
                  'ENTREGA/DEFINICAO_DO_PROBLEMA.md',
                  'ENTREGA/RESULTADOS.md',
                  # o pipeline que gera o que esta no relatorio
                  'campanha/roda_tudo.py',
                  'campanha/ativ2_alpha_naca.py',
                  'campanha/curvas_relatorio.py',
                  'campanha/figuras_relatorio.py',
                  'campanha/tabelas_do_relatorio.py',
                  # figuras e dados numericos do relatorio (o texto fica no
                  # Overleaf; aqui ficam so os entregaveis que o codigo gera)
                  'relatorio/images/05_transonico/mach.png',
                  'relatorio/tables/tab3_aerofolios.csv']
    print('\nconferencia:')
    for e in essenciais:
        ok = any(e in n for n in nomes)
        print(f'  {"ok " if ok else "FALTA"} {e}')
    figs = [n for n in nomes if n.endswith('.png')]
    print(f'  {len(figs)} figuras')
    return 0


if __name__ == '__main__':
    sys.exit(main())
