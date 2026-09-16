import os
import subprocess
import sys
import time

AQUI = os.path.dirname(os.path.abspath(__file__))
RELATORIO = os.path.normpath(os.path.join(AQUI, '..', '..', 'relatorio'))

# (script, argumentos, descricao, custo estimado, e_caro)
ETAPAS = [
    ('otimiza_secao.py', ['meio'],
     'otimizacao SLSQP da MAC', '~18 min', True),
    ('ativ2_alpha_naca.py', [],
     'item 2: alpha do NACA 1411 no c_l de projeto', '~6 min', True),
    ('curvas_relatorio.py', [],
     'curvas do Euler e do XFoil -> CSV', '~8 min', True),
    ('figuras_relatorio.py', [],
     'as 9 figuras, a partir dos CSVs', '~5 s', False),
    ('tabelas_do_relatorio.py', [],
     'as 6 tabelas, a partir dos CSVs', '~5 s', False),
]


def roda(script, args):
    print(f'\n{"=" * 68}\n>>> {script} {" ".join(args)}\n{"=" * 68}', flush=True)
    t0 = time.time()
    r = subprocess.run([sys.executable, '-u', os.path.join(AQUI, script)]
                       + args, cwd=AQUI)
    dt = time.time() - t0
    if r.returncode != 0:
        print(f'\n!!! {script} falhou (codigo {r.returncode}) '
              f'depois de {dt / 60:.1f} min', flush=True)
        return False
    print(f'--- {script} ok ({dt / 60:.1f} min)', flush=True)
    return True


def main():
    flags = [a for a in sys.argv[1:] if a.startswith('--')]

    if '--listar' in flags:
        print('etapas do pipeline:\n')
        for i, (s, a, d, c, caro) in enumerate(ETAPAS, 1):
            marca = '$' if caro else ' '
            print(f'  {i}. {marca} {s:26s} {c:>8s}  {d}')
        print('\n  ($ = chama solver)')
        return 0

    etapas = ETAPAS
    if '--so-desenho' in flags:
        etapas = [e for e in ETAPAS if not e[4]]
    elif '--pular-otim' in flags:
        etapas = [e for e in ETAPAS if e[0] != 'otimiza_secao.py']

    print(f'pipeline do relatorio -- {len(etapas)} etapa(s)')
    print(f'saida em {RELATORIO}')
    t0 = time.time()
    for script, args, desc, custo, _ in etapas:
        if not roda(script, args):
            print('\npipeline interrompido.', flush=True)
            return 1

    print(f'\n{"=" * 68}')
    print(f'pipeline completo em {(time.time() - t0) / 60:.1f} min')
    img = os.path.join(RELATORIO, 'images')
    tab = os.path.join(RELATORIO, 'tables')
    n_img = sum(len([f for f in fs if f.endswith('.png')])
                for _, _, fs in os.walk(img)) if os.path.isdir(img) else 0
    n_tab = len([f for f in os.listdir(tab)
                 if f.endswith('.csv')]) if os.path.isdir(tab) else 0
    print(f'  {n_img} figuras em relatorio/images/')
    print(f'  {n_tab} tabelas em relatorio/tables/')
    print('=' * 68)
    return 0


if __name__ == '__main__':
    sys.exit(main())
