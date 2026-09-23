'''
Escreve os arquivos do AVL a partir do dicionario de aeronave.py.

Deslocamentos em Z (exigencia do professor, para manter distancia entre os
paineis): asa -1,20 m, EH +1,85 m, EV +0,85 m, nacele acompanhando a asa.
Nao sao alteracao de projeto e existem so aqui.

Winglet: superficie propria no COMPONENT 1 (mesmo componente da asa), com
toe zero; geometria do designTool (Raymer, Fig. 7.34): vertical, altura e
corda de raiz iguais a corda da ponta, afilamento 0,21, bordo de fuga reto.

Nacele: anel sustentador no COMPONENT 1, como no 737.avl da disciplina.
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


def _corpo(av, malha):
    f = av['fuselagem']
    return ('#' + '-'*64 + f'\nBODY\nFuselage\n# Nbody Bspace\n{malha["nbody"]} 1.0\n'
            f'SCALE\n{f["L"]:.4f} {f["D"]:.4f} {f["D"]:.4f}\n'
            'BFILE\nfuselage_nondim.dat\n')


def _nacele(av):
    n = av['nacele']
    txt = ('#' + '-'*64 + '\nSURFACE\nNacelle\n#Nchordwise  Cspace   Nspanwise  Sspace\n'
           '6            1.0      12          0.0\nCOMPONENT\n1\nYDUPLICATE\n0.0\n'
           f'SCALE\n{n["L"]:.4f}  {n["D"]/2:.4f}  {n["D"]/2:.4f}\n'
           f'TRANSLATE\n{n["x"]:.4f}  {n["y"]:.4f}  {n["z"] + DZ["asa"]:.4f}\n')
    for y, z in _ANEL:
        txt += ('\nSECTION\n#Xle   Yle    Zle      Chord   Ainc  Nspanwise  Sspace\n'
                f' 0.00  {y:.3f}  {z:.3f}  1.0  0.  1  0.\n')
    return txt


def escreve_avl(av, caminho, cg='aft', torcao=None, malha=None, etas=ETAS,
                so_asa=False, diedro=True, cdp=None):
    '''
    Escreve o .avl em avl/<caminho> e devolve `caminho` (relativo a avl/).
    so_asa=True: so a asa, sem winglet, controles, corpos ou CDp
    (verificacao eliptica).
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
        txt += _corpo(av, malha) + _nacele(av)

    destino = os.path.join(AQUI, caminho)
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    with open(destino, 'w', encoding='utf-8') as f:
        f.write(txt)
    return caminho
