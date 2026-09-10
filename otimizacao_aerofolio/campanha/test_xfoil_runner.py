'''
Testes de regressao da interface com o XFoil.

O teste central e o do ponto de bordo de fuga duplicado: o XFoil encerra
silenciosamente (codigo de retorno 0, sem mensagem de erro) quando o arquivo
de coordenadas repete o ponto do bordo de fuga, porque isso gera um painel de
comprimento zero. Como a falha e muda, sem este teste ela reaparece calada.

Rodar com:  python -m pytest test_xfoil_runner.py -v
        ou:  python test_xfoil_runner.py
'''

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xfoil_runner as xr
import descritores as dsc

# NACA 1411 -- ponto de partida da otimizacao do roteiro
AL_1411 = [-0.1489439, -0.10330027, -0.10305128, -0.10514982]
AU_1411 = [0.16146332, 0.18349204, 0.14126241, 0.18194397]


def test_descritores_reproduzem_o_naca_1411():
    '''
    O NACA 1411 tem, por definicao: arqueamento 1% em x/c = 0.4, espessura
    11% em x/c = 0.30. Serve de aferimento dos descritores -- e trava o bug
    de indexacao do x_max_thickness do cstfoil, que devolve 0.117.
    '''
    d = dsc.descritores(AU_1411, AL_1411)
    assert abs(d['t_max'] - 0.11) < 0.01, f"t_max = {d['t_max']:.4f}"
    assert abs(d['x_tmax'] - 0.30) < 0.03, \
        f"x_tmax = {d['x_tmax']:.4f}, esperado ~0.30 (o cstfoil devolve 0.117)"
    assert abs(d['c_max'] - 0.01) < 0.004, f"c_max = {d['c_max']:.4f}"
    assert abs(d['x_cmax'] - 0.40) < 0.05, f"x_cmax = {d['x_cmax']:.4f}"


def test_delta_y_do_naca_simetrico():
    '''
    Perfis NACA de 4 digitos com mesma espessura tem o mesmo nariz: o Delta y
    do 0012 e do 1411 devem ser proximos, e da ordem de 3% da corda.
    '''
    au0012 = [0.16941, 0.15145, 0.13904, 0.13988]
    dy = dsc.delta_y(au0012)
    assert 2.0 < dy < 5.0, f'Delta y fora da faixa de handbook: {dy:.3f}'


def test_export_nao_repete_bordo_de_fuga(tmp_path=None):
    '''O arquivo exportado nao pode ter o ponto do BF tres vezes.'''
    import tempfile
    d = tmp_path or tempfile.mkdtemp()
    path = os.path.join(str(d), 'foil.dat')

    foil = xr.cst_coords(AU_1411, AL_1411)
    xr._export_quiet(foil, path)

    coords = []
    with open(path) as fid:
        for line in fid.readlines()[1:]:      # pula o titulo
            parts = line.split()
            if len(parts) == 2:
                coords.append((float(parts[0]), float(parts[1])))

    # o contorno fecha: primeiro ponto igual ao ultimo (convencao do XFoil)
    assert np.allclose(coords[0], coords[-1], atol=1e-12), \
        'o contorno deve fechar no bordo de fuga'

    # mas nao pode haver pontos consecutivos coincidentes (painel de comprimento zero)
    arr = np.array(coords)
    d2 = np.sum(np.diff(arr, axis=0) ** 2, axis=1)
    assert np.all(d2 > 1e-24), \
        f'ha {np.sum(d2 <= 1e-24)} painel(eis) de comprimento zero no arquivo'


def test_xfoil_converge_e_gera_polar():
    '''O XFoil precisa rodar ate o fim e devolver uma polar nao vazia.'''
    res = xr.clmax_cst(AU_1411, AL_1411, Re=9.0e6, Mach=0.0,
                       alpha_seq=(0.0, 6.0, 2.0))
    polar = res['polar']

    assert len(polar['CL']) > 0, \
        f"polar vazia -- o XFoil nao chegou ao fim. stdout final:\n{polar.get('stdout','')[-800:]}"
    assert np.isfinite(res['clmax']), 'clmax deveria ser finito'
    # NACA 1411 a Re 9e6: cl entre 0 e 1 nessa faixa de alpha
    assert 0.0 < res['clmax'] < 1.5, f"cl fora do esperado: {res['clmax']}"


def test_captura_o_pico_de_clmax():
    '''
    Com faixa de alpha suficiente, o estol precisa ser CAPTURADO (o cl sobe,
    atinge o pico e cai). Se a situacao for 'nao_atingido', cl_max e so um
    limite inferior e o DOE inteiro fica viesado para baixo.
    '''
    res = xr.clmax_cst(AU_1411, AL_1411, Re=4.2e7, Mach=0.0,
                       alpha_seq=(0.0, 30.0, 0.5))
    assert res['situacao'] == 'queda', \
        f"esperava capturar o pico, veio '{res['situacao']}' " \
        f"(clmax={res['clmax']:.4f} em alpha={res['alpha_clmax']:.1f})"
    assert res['alpha_clmax'] < 28.0, \
        'o pico ficou colado no fim da varredura; aumente alpha_seq'


def _polar_sintetica(alpha, CL):
    a = np.asarray(alpha, dtype=float)
    c = np.asarray(CL, dtype=float)
    z = np.zeros_like(a)
    return {'alpha': a, 'CL': c, 'CD': z, 'CDp': z, 'CM': z}


def test_classifica_falha_de_convergencia_como_nao_confiavel():
    '''
    Polar truncada no meio da faixa linear (cl ainda subindo na inclinacao
    cheia) e falha numerica do XFoil, NAO estol. Se isso passar como cl_max
    valido, entra ruido no DOE -- foi o bug que estragou o primeiro corte.
    '''
    a = np.arange(0.0, 11.0, 0.5)
    cl = 0.12 + 0.11 * a                      # reta, sem sinal de estol
    res = xr.clmax_from_polar(_polar_sintetica(a, cl), (0.0, 30.0, 0.5))
    assert res['situacao'] == 'falha_conv', res['situacao']
    assert not res['confiavel']


def test_classifica_estol_real_como_confiavel():
    '''Polar que sobe, atinge o pico e cai: estol capturado.'''
    # pico em a = 16.75 (0.11 - 0.04*(a-14) = 0), bem dentro da faixa varrida
    a = np.arange(0.0, 22.5, 0.5)
    cl = 0.12 + 0.11 * a - 0.02 * np.maximum(0.0, a - 14.0) ** 2
    res = xr.clmax_from_polar(_polar_sintetica(a, cl), (0.0, 30.0, 0.5))
    assert res['situacao'] == 'queda', res['situacao']
    assert res['confiavel']
    assert 15.0 < res['alpha_clmax'] < 19.0, res['alpha_clmax']


def test_repeticao_recupera_ponto_que_falhava():
    '''
    Perfil que na varredura de passo 0.5 dava cl_max ~1.34 por falha de
    convergencia, quando os vizinhos davam ~2.0. Com passo refinado tem de
    fechar acima de 1.9.
    '''
    k = 1.4210526315789473                    # idx 9 do corte controlado
    au = list(np.array(AU_1411) * np.array([k, 1, 1, 1]))
    al = list(np.array(AL_1411) * np.array([k, 1, 1, 1]))
    res = xr.clmax_cst(au, al, Re=4.2e7, Mach=0.266,
                       alpha_seq=(0.0, 30.0, 0.5), tentativas=3)
    assert res['confiavel'], res['situacao']
    assert res['clmax'] > 1.9, \
        f"clmax = {res['clmax']:.4f} (esperado ~2.0); a repeticao nao atuou"


def test_raio_bordo_ataque_bate_com_naca():
    '''r_LE/c = A0^2/2 deve reproduzir o raio teorico do NACA de 4 digitos.'''
    r_cst = xr.le_radius(AU_1411[0])
    r_naca = 1.1019 * 0.11 ** 2          # NACA 4 digitos, t/c = 0.11
    assert abs(r_cst - r_naca) / r_naca < 0.05, \
        f'r_LE do CST ({r_cst:.5f}) longe do teorico ({r_naca:.5f})'


if __name__ == '__main__':
    falhas = 0
    for nome, fn in sorted(globals().items()):
        if nome.startswith('test_') and callable(fn):
            try:
                fn()
                print(f'PASS  {nome}')
            except AssertionError as exc:
                falhas += 1
                print(f'FALHA {nome}\n      {exc}')
            except Exception as exc:                      # noqa: BLE001
                falhas += 1
                print(f'ERRO  {nome}\n      {type(exc).__name__}: {exc}')
    print(f'\n{falhas} falha(s)')
    sys.exit(1 if falhas else 0)
