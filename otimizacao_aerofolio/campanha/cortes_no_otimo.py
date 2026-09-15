'''
Cortes 1-a-1 no otimo: QUAL restricao barra cada direcao?

O Lab 02 fez essa figura para a aeronave e ela virou a narrativa do relatorio
("SM_aft barra AR baixo, tank_excess barra raiz fina..."). Aqui e o
equivalente para o perfil.

Saber quais restricoes estao ativas diz POUCO sozinho. O que explica o otimo
e o que acontece ao tentar andar em cada variavel: o arrasto piora, ou alguma
restricao e violada? Se para toda direcao houver uma das duas coisas, o ponto
e de fato um otimo local bem caracterizado.

Para cada coeficiente CST, varremos um intervalo em torno do otimo mantendo
os demais fixos, re-trimando alpha para segurar o c_l no alvo (um passo de
Newton com a derivada do adjunto), e registramos c_d e as restricoes.

Rodar com:  python cortes_no_otimo.py [estacao]
'''

import os
import sys
from multiprocessing import Pool

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
PACOTE = os.path.normpath(os.path.join(AQUI, '..', 'eulerblock', 'package'))
sys.path.insert(0, PACOTE)
sys.path.insert(0, AQUI)

import descritores as dsc
import otimiza_secao as ot
from analisa_otimos import RES, carrega

PASSOS = np.array([-0.06, -0.03, -0.015, 0.0, 0.015, 0.03, 0.06])


def avalia(args):
    '''Avalia um ponto, re-trimando alpha para manter o c_l no alvo.'''
    import shutil
    import tempfile

    from eulerblock import euler_mod as eb

    Al, Au, alpha0, cl_alvo, cfl, rot = args
    tmp = tempfile.mkdtemp(prefix='corte_')
    cwd = os.getcwd()
    alpha = alpha0
    try:
        os.chdir(tmp)
        r = None
        for _ in range(3):
            r = eb.run_cst(Al, Au, ot.NCHORD, alpha, ot.MACH_N,
                           gamma=1.4, order=2, iter=ot.ITER, dt=ot.DT,
                           CFL=cfl, use_local_dt=1, res_NK=ot.RES_NK,
                           res_tol=ot.RES_TOL, reinitialize=0,
                           adj_funcs=['cl_jlow', 'cd_jlow'], plot=False,
                           NJ=ot.NJ, s0=ot.S0)
            if not np.isfinite(r['CL']):
                break
            erro = r['CL'] - cl_alvo
            if abs(erro) < 2e-3:
                break
            dcl = r['grads']['cl_jlow']['alpha']
            if not np.isfinite(dcl) or abs(dcl) < 1e-6:
                break
            alpha = float(np.clip(alpha - erro / dcl,
                                  ot.ALPHA_MIN, ot.ALPHA_MAX))
        g = dsc.descritores(Au, Al)
        return {'rot': rot, 'CL': float(r['CL']), 'CD': float(r['CD']),
                'alpha': alpha, 'tmax': g['t_max'], 'tmin': g['t_min'],
                't01': ot.t01_de(Al, Au)}
    except Exception as exc:                                   # noqa: BLE001
        return {'rot': rot, 'CL': np.nan, 'CD': np.nan, 'alpha': alpha,
                'tmax': np.nan, 'tmin': np.nan, 't01': np.nan,
                'erro': str(exc)[:60]}
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    estacao = sys.argv[1] if len(sys.argv) > 1 else 'meio'
    d = carrega(f'otim_{estacao}')
    if d is None:
        print(f'{estacao}: sem resultado')
        return 1

    est = ot.ESTACOES[estacao]
    cfl = est.get('cfl', 0.20)
    lim = ot.limiar_bluntez(estacao)
    Al0, Au0, a0 = ot.desmonta(d['xx_otimo'])
    cd0 = d['hist']['CD'][d['i_otimo']]

    tarefas = []
    for i in range(ot.NVAR):
        for p in PASSOS:
            if p == 0.0:
                continue
            Al, Au = Al0.copy(), Au0.copy()
            Al[i] = np.clip(Al[i] + p, ot.AL_LOWER[i], ot.AL_UPPER[i])
            if abs(Al[i] - Al0[i]) < 1e-9:
                continue
            tarefas.append((Al, Au, a0, est['cl_ref'], cfl, (f'Al{i+1}', p)))
    for i in range(ot.NVAR):
        for p in PASSOS:
            if p == 0.0:
                continue
            Al, Au = Al0.copy(), Au0.copy()
            Au[i] = np.clip(Au[i] + p, ot.AU_LOWER[i], ot.AU_UPPER[i])
            if abs(Au[i] - Au0[i]) < 1e-9:
                continue
            tarefas.append((Al, Au, a0, est['cl_ref'], cfl, (f'Au{i+1}', p)))

    print(f'Cortes 1-a-1 no otimo da estacao {estacao}')
    print(f'  c_d do otimo = {cd0:.6f}   limiar de bluntez = {lim:.4f}')
    print(f'  {len(tarefas)} avaliacoes do Euler\n', flush=True)

    with Pool(4) as pool:
        saidas = pool.map(avalia, tarefas)

    por_var = {}
    for s in saidas:
        var, p = s['rot']
        por_var.setdefault(var, []).append((p, s))

    linhas = []
    print(f'{"var":5s} {"passo":>7s} {"c_d":>10s} {"vs otimo":>9s} '
          f'{"t/c":>8s} {"bluntez":>8s}  barrado por')
    for var in sorted(por_var):
        for p, s in sorted(por_var[var]):
            if not np.isfinite(s['CD']):
                print(f'{var:5s} {p:+7.3f} {"divergiu":>10s}')
                continue
            bl = ot.bluntez(s['t01'], s['tmax'])
            viola = []
            if s['tmax'] < est['tc_ref'] - 1e-4:
                viola.append('espessura')
            if s['tmin'] > ot.MINT_MAX + 1e-4:
                viola.append('BF grosso')
            if bl < lim - 1e-4:
                viola.append('bluntez')
            if abs(s['CL'] - est['cl_ref']) > 5e-3:
                viola.append('c_l fora do alvo')
            marca = ', '.join(viola) if viola else (
                'nada -- c_d piorou' if s['CD'] > cd0 else 'NADA: c_d MELHOROU')
            print(f'{var:5s} {p:+7.3f} {s["CD"]:10.6f} '
                  f'{(s["CD"]/cd0-1)*100:+8.2f}% {s["tmax"]:8.5f} '
                  f'{bl:8.5f}  {marca}')
            linhas.append({'variavel': var, 'passo': p, 'CD': s['CD'],
                           'CL': s['CL'], 'tmax': s['tmax'], 'bluntez': bl,
                           'barrado_por': marca})
        print()

    cam = os.path.join(RES, f'cortes_{estacao}.csv')
    with open(cam, 'w') as fid:
        cols = list(linhas[0])
        fid.write(','.join(cols) + '\n')
        for l in linhas:
            fid.write(','.join(
                f'{l[c]:.8g}' if isinstance(l[c], float) else str(l[c])
                for c in cols) + '\n')

    melhores = [l for l in linhas if 'MELHOROU' in l['barrado_por']]
    print('=' * 72)
    if melhores:
        print(f'ATENCAO: {len(melhores)} direcao(oes) melhoram o c_d sem')
        print('violar restricao. O ponto NAO e um otimo local bem resolvido:')
        for l in melhores[:5]:
            print(f'  {l["variavel"]} {l["passo"]:+.3f}: '
                  f'c_d {l["CD"]:.6f} contra {cd0:.6f}')
    else:
        print('Toda direcao ou piora o c_d ou viola alguma restricao.')
        print('O ponto e um otimo local bem caracterizado.')
    print(f'\ngravado em {cam}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
