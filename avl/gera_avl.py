'''
Escreve os arquivos do AVL a partir do dicionario de aeronave.py.

Deslocamentos em Z (exigencia do professor, para manter distancia entre os
paineis): asa -1,20 m, EH +1,85 m, EV +0,85 m, nacele acompanhando a asa.
Nao sao alteracao de projeto e existem so aqui.

Winglet: superficie propria no COMPONENT 1 (mesmo componente da asa), com
toe zero; geometria do designTool (Raymer, Fig. 7.34): vertical, altura e
corda de raiz iguais a corda da ponta, afilamento 0,21, bordo de fuga reto.

Nacele: anel sustentador em COMPONENT proprio (default 2), separado da asa.
Na mesma COMPONENT que a asa (como no 737.avl da disciplina) o anel faz
CDff oscilar +-0,2 count com a malha da envergadura da asa (paineis da asa
e do anel compartilhando o mesmo componente de sustentacao no AVL); em
componente proprio essa oscilacao some.

Fuselagem: BFILE reamostrado (fuselage_reamostrada.dat, gerado por
reamostra_fuselagem) em vez do fuselage_nondim.dat original. O arquivo
original tem so 41 pontos e fecha o bordo de fuga de forma abrupta (a
superficie inferior avanca ~0,5 m em x nos ultimos 0,5% do comprimento e
"anda para tras"); isso faz a spline de comprimento de arco do AVL dar um
laco perto da cauda, e CDff/it variam com Nbody em vez de convergir. O
reamostrado interpola linearmente cada lado (superior/inferior, separados
no x minimo) em ~81 pontos com espacamento cosseno em x, preservando
x em [0, 1] e os pontos extremos; com ele CDff/it convergem a partir de
nbody ~= 80.
'''

import json
import os

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))

DZ = {'asa': -1.20, 'EH': +1.85, 'EV': +0.85}
TAPER_WINGLET = 0.21

# Estacoes da asa (fracao da semienvergadura) e perfis do Lab 03
ETAS = [0.0, 0.1011, 0.22, 0.32, 0.398, 0.48, 0.56, 0.70, 0.82, 0.90, 1.0]
PERFIS = ['raiz.dat', 'raiz.dat', 'gerado/eta_0220.dat', 'gerado/eta_0320.dat',
          'gerado/eta_0398.dat', 'gerado/eta_0480.dat', 'gerado/eta_0560.dat',
          'gerado/eta_0700.dat', 'gerado/eta_0820.dat', 'ponta.dat',
          'ponta.dat']
ETA_AILERON = (0.56, 0.90)

MALHA_PADRAO = {'asa_nc': 8, 'asa_ns': 20, 'winglet_ns': 8, 'eh_nc': 16,
                'eh_ns': 10, 'ev_nc': 8, 'ev_ns': 6, 'nbody': 20}
ARQ_MALHA = os.path.join(AQUI, 'resultados', 'malha_adotada.json')

_ANEL = [(0.0, 1.0), (0.5, 0.866), (0.866, 0.5), (1.0, 0.0), (0.866, -0.5),
         (0.5, -0.866), (0.0, -1.0), (-0.5, -0.866), (-0.866, -0.5),
         (-1.0, 0.0), (-0.866, 0.5), (-0.5, 0.866), (0.0, 1.0)]


def reamostra_fuselagem(origem='fuselage_nondim.dat', destino='fuselage_reamostrada.dat',
                        n=81):
    '''
    Reamostra o contorno da fuselagem (origem, formato BFILE do AVL: linha
    de cabecalho + pontos x,z fechando TE -> superficie de cima -> LE ->
    superficie de baixo -> TE) em `n` pontos por lado, com espacamento
    cosseno em x e interpolacao linear de cada lado sobre o contorno
    original. Evita a spline de comprimento de arco do AVL "laçar" perto
    da cauda quando o arquivo original fecha o bordo de fuga com poucos
    pontos (ver docstring do modulo). Preserva x em [0, 1] e os pontos
    extremos (TE em x=1). Devolve `destino` (caminho relativo a avl/).
    '''
    origem_abs = os.path.join(AQUI, origem)
    destino_abs = os.path.join(AQUI, destino)
    if (os.path.exists(destino_abs)
            and os.path.getmtime(destino_abs) >= os.path.getmtime(origem_abs)):
        return destino

    with open(origem_abs) as f:
        cabecalho = f.readline().strip()
    d = np.loadtxt(origem_abs, skiprows=1)
    i0 = int(np.argmin(d[:, 0]))
    cima = d[:i0 + 1][::-1]      # LE -> TE, superficie de cima
    baixo = d[i0:]               # LE -> TE, superficie de baixo

    xs = 0.5*(1 - np.cos(np.linspace(0.0, np.pi, n)))
    xs[0], xs[-1] = 0.0, 1.0
    z_cima = np.interp(xs, cima[:, 0], cima[:, 1])
    z_baixo = np.interp(xs, baixo[:, 0], baixo[:, 1])

    with open(destino_abs, 'w', encoding='utf-8') as f:
        f.write(cabecalho + '\n')
        for x, z in zip(xs[::-1], z_cima[::-1]):
            f.write(f'{x:.6f} {z:.6f}\n')
        for x, z in zip(xs[1:], z_baixo[1:]):
            f.write(f'{x:.6f} {z:.6f}\n')
    return destino


def malha_adotada():
    '''Malha da convergencia (resultados/malha_adotada.json) ou a padrao.'''
    m = dict(MALHA_PADRAO)
    if os.path.exists(ARQ_MALHA):
        with open(ARQ_MALHA) as f:
            m.update(json.load(f))
    return m


def perfil(eta):
    i = int(np.argmin([abs(eta - e) for e in ETAS]))
    return 'airfoils/' + PERFIS[i]


def secoes_asa(av, torcao=None, etas=ETAS, diedro=True):
    w = av['asa']
    etas = list(etas)
    torcao = np.zeros(len(etas)) if torcao is None else np.asarray(torcao, float)
    if len(torcao) != len(etas):
        raise ValueError('torcao precisa ter um valor por estacao')
    out = []
    for eta, tw in zip(etas, torcao):
        z_real = w['zr'] + eta*(w['zt'] - w['zr']) if diedro else w['zr']
        out.append({'eta': eta, 'y': eta*w['yt'],
                    'x': w['xr'] + eta*(w['xt'] - w['xr']),
                    'z': z_real + DZ['asa'],
                    'c': w['cr'] + eta*(w['ct'] - w['cr']),
                    'ainc': float(tw), 'perfil': perfil(eta),
                    'aileron': ETA_AILERON[0] - 1e-9 <= eta <= ETA_AILERON[1] + 1e-9})
    return out


def _secao(x, y, z, c, ainc=0.0):
    return ('\nSECTION\n#Xle      Yle      Zle      Chord    Ainc\n'
            f'{x:.4f}  {y:.4f}  {z:.4f}  {c:.4f}  {ainc:.4f}\n')


def _superficie(nome, nc, ns, extra=''):
    return ('#' + '-'*64 + f'\nSURFACE\n{nome}\n#Nchordwise Cspace Nspanwise Sspace\n'
            f'{nc} 1.0 {ns} 1.0\n{extra}')


def _asa(secoes, malha, controles):
    txt = _superficie('Wing', malha['asa_nc'], malha['asa_ns'],
                      'YDUPLICATE\n0.0\nANGLE\n0.0\n')
    for s in secoes:
        txt += _secao(s['x'], s['y'], s['z'], s['c'], s['ainc'])
        txt += f'AFILE\n{s["perfil"]}\nCLAF\n1.0000\n'
        if controles and s['aileron']:
            txt += 'CONTROL\naileron  1.0   0.73    0. 0. 0.    -1.0\n'
    return txt


def _winglet(secoes, malha):
    p = secoes[-1]
    ct = p['c']
    txt = _superficie('Winglet', malha['asa_nc'], malha['winglet_ns'],
                      'COMPONENT\n1\nYDUPLICATE\n0.0\nANGLE\n0.0\n')
    txt += _secao(p['x'], p['y'], p['z'], ct) + 'AFILE\nairfoils/ponta.dat\nCLAF\n1.0000\n'
    txt += (_secao(p['x'] + ct - TAPER_WINGLET*ct, p['y'], p['z'] + ct,
                   TAPER_WINGLET*ct) + 'AFILE\nairfoils/ponta.dat\nCLAF\n1.0000\n')
    return txt


def _eh(av, malha):
    h = av['EH']
    txt = _superficie('Horizontal tail', malha['eh_nc'], malha['eh_ns'],
                      'YDUPLICATE\n0.0\nANGLE\n0.0\n')
    for x, y, z, c in ((h['xr'], 0.0, h['zr'], h['cr']),
                       (h['xt'], h['yt'], h['zt'], h['ct'])):
        txt += _secao(x, y, z + DZ['EH'], c)
        txt += ('NACA\n0010\nCLAF\n1.0000\nDESIGN\nit      1.0\n'
                'CONTROL\nelevator 1.0   0.70    0. 0. 0.     1.0\n')
    return txt


def _ev(av, malha):
    v = av['EV']
    txt = _superficie('Vertical tail', malha['ev_nc'], malha['ev_ns'], 'ANGLE\n0.0\n')
    for x, z, c in ((v['xr'], v['zr'], v['cr']), (v['xt'], v['zt'], v['ct'])):
        txt += _secao(x, 0.0, z + DZ['EV'], c)
        txt += ('NACA\n0010\nCLAF\n1.0000\n'
                'CONTROL\nrudder   1.0   0.70    0. 0. 0.     0.0\n')
    return txt


def _corpo(av, malha, fuselagem='reamostrada'):
    f = av['fuselagem']
    if fuselagem == 'reamostrada':
        bfile = reamostra_fuselagem()
    elif fuselagem == 'original':
        bfile = 'fuselage_nondim.dat'
    else:
        raise ValueError("fuselagem precisa ser 'reamostrada' ou 'original'")
    return ('#' + '-'*64 + f'\nBODY\nFuselage\n# Nbody Bspace\n{malha["nbody"]} 1.0\n'
            f'SCALE\n{f["L"]:.4f} {f["D"]:.4f} {f["D"]:.4f}\n'
            f'BFILE\n{bfile}\n')


def _nacele(av, componente=2):
    n = av['nacele']
    txt = ('#' + '-'*64 + '\nSURFACE\nNacelle\n#Nchordwise  Cspace   Nspanwise  Sspace\n'
           f'6            1.0      12          0.0\nCOMPONENT\n{componente}\nYDUPLICATE\n0.0\n'
           f'SCALE\n{n["L"]:.4f}  {n["D"]/2:.4f}  {n["D"]/2:.4f}\n'
           f'TRANSLATE\n{n["x"]:.4f}  {n["y"]:.4f}  {n["z"] + DZ["asa"]:.4f}\n')
    for y, z in _ANEL:
        txt += ('\nSECTION\n#Xle   Yle    Zle      Chord   Ainc  Nspanwise  Sspace\n'
                f' 0.00  {y:.3f}  {z:.3f}  1.0  0.  1  0.\n')
    return txt


def escreve_avl(av, caminho, cg='aft', torcao=None, malha=None, etas=ETAS,
                so_asa=False, diedro=True, cdp=None, nacele_componente=2,
                fuselagem='reamostrada'):
    '''
    Escreve o .avl em avl/<caminho> e devolve `caminho` (relativo a avl/).
    so_asa=True: so a asa, sem winglet, controles, corpos ou CDp
    (verificacao eliptica).
    nacele_componente: COMPONENT do anel da nacele (default 2, separado da
    asa; ver docstring do modulo). fuselagem: 'reamostrada' (default) ou
    'original' (fuselage_nondim.dat, so para reproduzir o modelo congelado
    de referencia).
    '''
    malha = malha_adotada() if malha is None else malha
    secoes = secoes_asa(av, torcao, etas, diedro)
    xref = av['xcg_fwd'] if cg == 'fwd' else av['xcg_aft']
    cdp = (0.0 if so_asa else av['CD0']) if cdp is None else cdp

    txt = (f'Ararinha Lc_h={av["Lc_h"]:.4f} CG {cg}{" so asa" if so_asa else ""}\n'
           f'{av["M"]:.4f}          # Mach\n0 0 0.0       # iYsym iZsym Zsym\n'
           f'{av["Sref"]:.4f} {av["Cref"]:.4f} {av["Bref"]:.4f}   # Sref Cref Bref\n'
           f'{xref:.4f} 0.0 0.0   # Xref Yref Zref\n{cdp:.5f}       # CDp\n')
    txt += _asa(secoes, malha, controles=not so_asa)
    if not so_asa:
        txt += _winglet(secoes, malha) + _eh(av, malha) + _ev(av, malha)
        txt += _corpo(av, malha, fuselagem) + _nacele(av, nacele_componente)

    destino = os.path.join(AQUI, caminho)
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    with open(destino, 'w', encoding='utf-8') as f:
        f.write(txt)
    return caminho
