'''
Interface em batch com o XFoil (Lab 03 - PRJ-23, equipe Ararinha).

O XFoil e dirigido por stdin, com os graficos desligados (PLOP/G F) para
rodar sem abrir janela. Cada analise acontece num diretorio temporario
proprio, porque o XFoil escreve o arquivo de polar no diretorio corrente --
isso permite rodar varias analises em paralelo sem que uma sobrescreva a
outra.

Uso principal:
    res = clmax_cst(Au, Al, Re=4.2e7, Mach=0.0)
    res['clmax'], res['alpha_clmax'], res['r_LE_upper'], ...

Escrito para a campanha de projeto de perfil da aeronave ararinha_23.
'''

import os
import shutil
import subprocess
import tempfile

import numpy as np

import airfoil_mod as am

# ---------------------------------------------------------------------------

XFOIL_EXE = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 '..', 'xfoil', 'xfoil.exe'))

# numero de paineis e concentracao no bordo de fuga sugeridos no roteiro
N_PANEL = 200
TE_BUNCH = 0.45

# malha de pontos para gerar as coordenadas do perfil a partir do CST
N_COORD = 161


def cst_coords(Au, Al, n=N_COORD):
    '''Coordenadas do perfil CST com espacamento cosseno (adensado nos bordos).'''
    x = (1 - np.cos(np.linspace(0, 1, n) * np.pi)) / 2
    return am.cstfoil(np.asarray(Au), np.asarray(Al), x)


def le_radius(A0):
    '''
    Raio de bordo de ataque adimensional da parametrizacao CST.

    Com a funcao de classe N1=0.5, N2=1.0, a superficie tende a
    y -> A[0]*sqrt(x) quando x -> 0. A parabola y^2 = A^2*x tem raio de
    curvatura no vertice igual a A^2/2, logo r_LE/c = A[0]^2/2.
    '''
    return A0 ** 2 / 2.0


def _write_commands(fid, cmds):
    fid.write('\n'.join(cmds) + '\n')


def run_xfoil(airfoil_file, Re, Mach=0.0, alpha_seq=(-2.0, 20.0, 0.25),
              n_crit=9.0, iter_max=200, workdir=None, timeout=300):
    '''
    Roda uma varredura de angulo de ataque no XFoil e devolve a polar.

    airfoil_file: caminho para o .dat no formato do XFoil
    Re: numero de Reynolds baseado na corda
    Mach: numero de Mach (0 desliga a correcao de compressibilidade)
    alpha_seq: (alpha_inicial, alpha_final, passo) em graus
    n_crit: fator de amplificacao do metodo e^N (9 = tunel limpo)

    Devolve um dicionario com arrays alpha, CL, CD, CDp, CM e metadados.
    '''

    if not os.path.isfile(XFOIL_EXE):
        raise FileNotFoundError(f'xfoil.exe nao encontrado em {XFOIL_EXE}')

    tmp = workdir or tempfile.mkdtemp(prefix='xfoil_')
    created_tmp = workdir is None

    try:
        local_foil = os.path.join(tmp, 'foil.dat')
        shutil.copyfile(airfoil_file, local_foil)
        polar_file = os.path.join(tmp, 'polar.txt')

        a0, a1, da = alpha_seq

        cmds = [
            'PLOP',       # opcoes de plotagem
            'G F',        # desliga os graficos
            '',           # volta ao menu principal
            'LOAD foil.dat',
            'PPAR',       # ajuste dos paineis, como recomenda o roteiro
            'T',
            f'{TE_BUNCH}',
            'N',
            f'{N_PANEL}',
            '',
            '',
            'OPER',
            f'ITER {iter_max}',
            'VPAR',       # parametros viscosos
            f'N {n_crit}',
            '',
        ]

        if Mach > 0:
            cmds.append(f'MACH {Mach}')

        cmds += [
            f'VISC {Re:.6g}',
            'PACC',           # acumula a polar num arquivo
            'polar.txt',
            '',               # sem arquivo de dump
            f'ASEQ {a0} {a1} {da}',
            'PACC',           # fecha a acumulacao
            '',               # sai do OPER
            'QUIT',
        ]

        script = '\n'.join(cmds) + '\n'

        proc = subprocess.run(
            [XFOIL_EXE], input=script, cwd=tmp, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=timeout)

        polar = _read_polar(polar_file)
        polar['returncode'] = proc.returncode
        polar['stdout'] = proc.stdout
        return polar

    except subprocess.TimeoutExpired:
        return {'alpha': np.array([]), 'CL': np.array([]), 'CD': np.array([]),
                'CDp': np.array([]), 'CM': np.array([]),
                'returncode': None, 'stdout': 'TIMEOUT', 'timeout': True}

    finally:
        if created_tmp:
            shutil.rmtree(tmp, ignore_errors=True)


def _read_polar(path):
    '''Le o arquivo de polar do XFoil (cabecalho de 12 linhas).'''
    empty = {'alpha': np.array([]), 'CL': np.array([]), 'CD': np.array([]),
             'CDp': np.array([]), 'CM': np.array([])}
    if not os.path.isfile(path):
        return empty

    rows = []
    with open(path) as fid:
        lines = fid.readlines()

    # os dados comecam depois da linha de tracos que segue o cabecalho
    start = None
    for ii, line in enumerate(lines):
        if line.strip().startswith('---'):
            start = ii + 1
            break
    if start is None:
        return empty

    for line in lines[start:]:
        parts = line.split()
        if len(parts) >= 5:
            try:
                rows.append([float(v) for v in parts[:5]])
            except ValueError:
                continue

    if not rows:
        return empty

    arr = np.array(rows)
    order = np.argsort(arr[:, 0])
    arr = arr[order]
    return {'alpha': arr[:, 0], 'CL': arr[:, 1], 'CD': arr[:, 2],
            'CDp': arr[:, 3], 'CM': arr[:, 4]}


def clmax_from_polar(polar, alpha_seq):
    '''
    Extrai cl_max da polar e classifica a confiabilidade do resultado.

    O XFoil, ao nao convergir num angulo, pula para o proximo. Esses buracos
    aparecem tambem na SUBIDA da curva (tropecos de convergencia) e nao
    significam estol -- por isso nao cortamos a polar neles. O criterio de
    estol e o comportamento fisico: o cl atinge um pico e depois cai.

    'situacao' diz o quanto confiar no numero:
      queda         -> depois do pico o cl cai pelo menos 2%: estol capturado
      topo_achatado -> a curva achatou no fim (inclinacao caiu para menos de
                       1/4 da linear): estamos no joelho, cl_max confiavel
      falha_conv    -> a varredura parou antes do fim com o cl AINDA SUBINDO
                       na inclinacao linear: o XFoil desistiu por problema
                       numerico, nao por estol -- NAO confiavel
      nao_atingido  -> chegou ao angulo final ainda subindo: cl_max e so um
                       LIMITE INFERIOR, aumente alpha_seq
      vazia         -> nenhum ponto convergiu

    A distincao entre 'queda'/'topo_achatado' e 'falha_conv' e o ponto
    delicado: nos dois casos a varredura termina cedo. O que separa e a
    INCLINACAO no fim. Sem isso, uma falha numerica do XFoil no meio da faixa
    linear entra no conjunto de dados como se fosse um cl_max baixo de
    verdade, e contamina qualquer regressao feita em cima.
    '''
    alpha, CL = polar['alpha'], polar['CL']
    a_end = alpha_seq[1]
    passo = abs(alpha_seq[2])

    vazio = {'clmax': np.nan, 'alpha_clmax': np.nan, 'n_pontos': 0,
             'alpha_max_convergido': np.nan, 'situacao': 'vazia',
             'confiavel': False, 'n_buracos': 0,
             'incl_linear': np.nan, 'incl_fim': np.nan, 'razao_incl': np.nan}
    if len(CL) < 3:
        return vazio

    n_buracos = int(np.sum(np.diff(alpha) > 1.5 * passo))

    imax = int(np.argmax(CL))
    cl_pico = float(CL[imax])
    caiu = bool(np.any(CL[imax + 1:] < 0.98 * cl_pico))

    # inclinacao da regiao linear (primeiros 6 graus convergidos)
    lin = alpha <= alpha[0] + 6.0
    incl_linear = (float(np.polyfit(alpha[lin], CL[lin], 1)[0])
                   if lin.sum() >= 3 else np.nan)
    # inclinacao nos ultimos pontos convergidos
    incl_fim = (float(np.polyfit(alpha[-4:], CL[-4:], 1)[0])
                if len(alpha) >= 4 else np.nan)

    razao = (incl_fim / incl_linear
             if np.isfinite(incl_linear) and incl_linear > 1e-6
             and np.isfinite(incl_fim) else np.nan)

    chegou_ao_fim = alpha[-1] >= a_end - 1e-6

    if caiu:
        situacao = 'queda'
    elif np.isfinite(razao) and razao < 0.25:
        situacao = 'topo_achatado'
    elif chegou_ao_fim:
        situacao = 'nao_atingido'
    else:
        situacao = 'falha_conv'

    return {'clmax': cl_pico,
            'alpha_clmax': float(alpha[imax]),
            'n_pontos': int(len(CL)),
            'alpha_max_convergido': float(alpha[-1]),
            'situacao': situacao,
            'confiavel': situacao in ('queda', 'topo_achatado'),
            'n_buracos': n_buracos,
            'incl_linear': incl_linear,
            'incl_fim': incl_fim,
            'razao_incl': float(razao) if np.isfinite(razao) else np.nan}


def clmax_cst(Au, Al, Re, Mach=0.0, alpha_seq=(-2.0, 20.0, 0.25),
              n_crit=9.0, iter_max=200, timeout=300, keep_dir=None,
              tentativas=3):
    '''
    Gera o perfil CST, roda o XFoil e devolve cl_max junto com as
    propriedades geometricas relevantes.

    Se a varredura nao capturar o pico (situacao 'falha_conv'), repete com
    passo de angulo pela metade e mais iteracoes viscosas: marchar mais
    devagar costuma atravessar o ponto em que o XFoil tropeca. Sem esse
    reforco, falhas numericas entram no conjunto de dados como cl_max baixos
    e contaminam a analise -- foi exatamente o que aconteceu no primeiro
    corte controlado deste DOE.
    '''
    airfoil = cst_coords(Au, Al)

    tmp = keep_dir or tempfile.mkdtemp(prefix='xfoil_')
    created = keep_dir is None
    try:
        foil_path = os.path.join(tmp, 'airfoil.dat')
        _export_quiet(airfoil, foil_path)

        seq, it = alpha_seq, iter_max
        for n in range(tentativas):
            polar = run_xfoil(foil_path, Re, Mach, seq, n_crit,
                              it, workdir=tmp, timeout=timeout)
            res = clmax_from_polar(polar, seq)
            res['tentativas'] = n + 1
            res['passo_final'] = seq[2]
            if res['situacao'] != 'falha_conv':
                break
            seq = (seq[0], seq[1], seq[2] / 2.0)
            it = int(it * 1.5)

        res.update({
            'r_LE_upper': le_radius(Au[0]),
            'r_LE_lower': le_radius(Al[0]),
            'Au0': float(Au[0]),
            'Al0': float(Al[0]),
            'max_thickness': float(np.real(airfoil['max_thickness'])),
            'x_max_thickness': float(np.real(airfoil['x_max_thickness'])),
            'min_thickness': float(np.real(airfoil['min_thickness'])),
            'max_camber': float(np.real(airfoil['max_camber'])),
            'x_max_camber': float(np.real(airfoil['x_max_camber'])),
            'polar': polar,
        })
        return res
    finally:
        if created:
            shutil.rmtree(tmp, ignore_errors=True)


def _export_quiet(airfoil, path, tol=1e-12):
    '''
    Escreve as coordenadas no formato do XFoil.

    Nao usamos o export_airfoil do airfoil_mod porque ele imprime em stdout e,
    com o padrao close_te=True, repete o ponto do bordo de fuga. O cstfoil ja
    devolve o contorno fechado (primeiro ponto = ultimo, que e a convencao do
    XFoil para bordo de fuga afiado), entao a copia extra vira um terceiro
    ponto coincidente -- um painel de comprimento zero. Diante disso o XFoil
    encerra em silencio, com codigo de retorno 0 e sem mensagem de erro, logo
    depois de imprimir "Clockwise ordering". Ver test_xfoil_runner.py.
    '''
    xf = np.real(np.asarray(airfoil['x_coord'], dtype=float))
    yf = np.real(np.asarray(airfoil['y_coord'], dtype=float))

    # descarta pontos consecutivos coincidentes, preservando o fechamento
    keep = [0]
    for ii in range(1, len(xf)):
        dx, dy = xf[ii] - xf[keep[-1]], yf[ii] - yf[keep[-1]]
        if dx * dx + dy * dy > tol ** 2:
            keep.append(ii)
    xf, yf = xf[keep], yf[keep]

    # fecha o contorno apenas se ele ainda nao estiver fechado
    dx, dy = xf[0] - xf[-1], yf[0] - yf[-1]
    fechar = dx * dx + dy * dy > tol ** 2

    with open(path, 'w') as fid:
        fid.write('AIRFOIL\n')
        for xi, yi in zip(xf, yf):
            fid.write('%.13f %.13f\n' % (xi, yi))
        if fechar:
            fid.write('%.13f %.13f\n' % (xf[0], yf[0]))
