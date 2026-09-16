'''
O poco de arrasto do perfil otimizado e estreito de verdade, ou e artefato?

A polar do item 7 (polares_finais.py) varre alpha de 2 em 2 graus com passo de
0,5 grau e mostra uma coisa esquisita na estacao do meio:

    alpha   1,889   2,389   2,889
    c_d    0,0186  0,0116  0,0205

Uma queda de 40% num unico ponto, justamente o de projeto, ladeada por dois
valores altos. E a curva c_l x alpha tem ali um salto de inclinacao de 22/rad
contra ~10/rad no resto da varredura -- o dobro.

Duas leituras possiveis, com consequencias bem diferentes:

  (a) e FISICA. Perfil supercritico mono-ponto: no ponto de projeto o
      escoamento e quase sem choque, e meio grau fora o choque se forma de
      supetao. O poco e estreito mesmo, e o perfil so serve se o c_l de
      projeto estiver certo -- o que depende da hipotese de distribuicao
      eliptica, que so o Lab 04 (torcao) confirma.

  (b) e NUMERICO. Os vizinhos nao convergiram e o solver devolveu um valor
      finito mas errado. Nesse caso o poco e uma miragem e a polar da entrega
      esta errada.

O jeito de separar as duas e olhar o RESIDUO, que o run_cst nao devolve. Aqui
capturamos a saida do executavel, registramos se cada ponto atingiu a
tolerancia ou estourou o limite de iteracoes, e varremos alpha com passo de
0,1 grau -- cinco vezes mais fino que a polar original.

Se for (a), a curva fina e suave e todo ponto converge.
Se for (b), os pontos ruins aparecem como nao-convergidos.

Rodar com:  python poco_de_arrasto.py [estacao]
'''

import os
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
PACOTE = os.path.normpath(os.path.join(AQUI, '..', 'eulerblock', 'package'))
sys.path.insert(0, PACOTE)
sys.path.insert(0, AQUI)

import otimiza_secao as ot
from analisa_otimos import RES, carrega

DALPHA = np.round(np.arange(-1.0, 1.001, 0.1), 3)   # graus, em torno do otimo
SAIDA_EXE = '_saida_exe.txt'


def _captura(cmd, *a, **k):
    '''Substitui subprocess.run para guardar o stdout do eulerblock.exe.'''
    p = subprocess.__dict__['_run_original'](cmd, capture_output=True,
                                             text=True)
    with open(SAIDA_EXE, 'w') as fid:
        fid.write(p.stdout or '')
    return p


def _le_convergencia():
    '''Le o stdout do solver: convergiu? em quantas iteracoes? que residuo?'''
    try:
        txt = open(SAIDA_EXE).read()
    except OSError:
        return {'convergiu': None, 'n_iter': np.nan, 'residuo': np.nan}
    convergiu = 'Reached desired residual tolerance' in txt
    n_iter, residuo = np.nan, np.nan
    for linha in reversed(txt.splitlines()):
        if 'Final iteration' in linha:
            partes = linha.replace(':', ' ').split()
            try:
                n_iter = float(partes[2])
                residuo = float(partes[-1].replace('E', 'e'))
            except (IndexError, ValueError):
                pass
            break
    return {'convergiu': convergiu, 'n_iter': n_iter, 'residuo': residuo}


def _um_ponto(args):
    from eulerblock import euler_mod as eb

    # instala a captura dentro do processo filho
    if '_run_original' not in subprocess.__dict__:
        subprocess.__dict__['_run_original'] = subprocess.run
    eb.subprocess.run = _captura

    Al, Au, alpha, cfl = args
    tmp = tempfile.mkdtemp(prefix='poco_')
    cwd = os.getcwd()
    try:
        os.chdir(tmp)
        r = eb.run_cst(Al, Au, ot.NCHORD, alpha, ot.MACH_N, gamma=1.4,
                       order=2, iter=ot.ITER, dt=ot.DT, CFL=cfl,
                       use_local_dt=1, res_NK=ot.RES_NK, res_tol=ot.RES_TOL,
                       reinitialize=0, adj_funcs=[], plot=False,
                       NJ=ot.NJ, s0=ot.S0)
        conv = _le_convergencia()
        mach = np.asarray(r['distrib']['Mach'])
        return dict(alpha=alpha, CL=float(r['CL']), CD=float(r['CD']),
                    CM=float(r['CM']), mach_max=float(np.nanmax(mach)),
                    **conv)
    except Exception as exc:                                   # noqa: BLE001
        return dict(alpha=alpha, CL=np.nan, CD=np.nan, CM=np.nan,
                    mach_max=np.nan, convergiu=False, n_iter=np.nan,
                    residuo=np.nan, erro=str(exc)[:60])
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    estacao = sys.argv[1] if len(sys.argv) > 1 else 'meio'
    d = carrega(f'otim_{estacao}')
    if d is None:
        print(f'{estacao}: sem resultado')
        return 1

    cfl = ot.ESTACOES[estacao].get('cfl', 0.20)
    Al, Au, a0 = ot.desmonta(d['xx_otimo'])
    cd0 = d['hist']['CD'][d['i_otimo']]
    alphas = a0 + np.radians(DALPHA)

    print(f'Poco de arrasto da estacao {estacao}')
    print(f'  otimo: alpha = {np.degrees(a0):.4f} deg, c_d = {cd0:.6f}')
    print(f'  {len(alphas)} pontos, passo de 0,1 deg '
          f'(a polar do item 7 usa 0,5)\n', flush=True)

    with Pool(4) as pool:
        saidas = pool.map(_um_ponto, [(Al, Au, a, cfl) for a in alphas])
    saidas.sort(key=lambda s: s['alpha'])

    print(f'{"d_alpha":>8s} {"c_l":>9s} {"c_d":>10s} {"vs otimo":>9s} '
          f'{"M_max":>7s} {"iter":>7s} {"residuo":>10s}  convergiu')
    for s, da in zip(saidas, DALPHA):
        if not np.isfinite(s['CD']):
            print(f'{da:+8.1f} {"falhou":>9s}')
            continue
        marca = 'sim' if s['convergiu'] else '*** NAO ***'
        print(f'{da:+8.1f} {s["CL"]:9.5f} {s["CD"]:10.6f} '
              f'{(s["CD"]/cd0-1)*100:+8.1f}% {s["mach_max"]:7.4f} '
              f'{s["n_iter"]:7.0f} {s["residuo"]:10.2e}  {marca}')

    cam = os.path.join(RES, f'poco_de_arrasto_{estacao}.csv')
    cols = ['alpha', 'CL', 'CD', 'CM', 'mach_max', 'convergiu', 'n_iter',
            'residuo']
    with open(cam, 'w') as fid:
        fid.write('d_alpha_deg,' + ','.join(cols) + '\n')
        for s, da in zip(saidas, DALPHA):
            fid.write(f'{da:.1f},' + ','.join(
                f'{s[c]:.8g}' if isinstance(s[c], float) else str(s[c])
                for c in cols) + '\n')

    ruins = [s for s in saidas if not s['convergiu']]
    finitos = [s for s in saidas if np.isfinite(s['CD'])]
    print(f'\n{"="*74}')
    print('LEITURA')
    print('=' * 74)
    print(f'  pontos que NAO convergiram: {len(ruins)} de {len(saidas)}')
    if ruins:
        print('  -> o poco estreito da polar do item 7 e, ao menos em parte,')
        print('     artefato numerico. A polar da entrega precisa ser refeita.')
    else:
        cds = np.array([s['CD'] for s in finitos])
        i = int(np.argmin(cds))
        print('  todo ponto convergiu -> o poco e FISICO.')
        print(f'  minimo da varredura fina: d_alpha = {DALPHA[i]:+.1f} deg, '
              f'c_d = {cds[i]:.6f}')
        # largura do poco: faixa em que c_d fica ate 10% acima do minimo
        dentro = DALPHA[cds <= 1.10 * cds[i]]
        if len(dentro):
            print(f'  largura do poco (c_d ate +10% do minimo): '
                  f'{dentro[0]:+.1f} a {dentro[-1]:+.1f} deg')
    print(f'\n  gravado em {cam}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
