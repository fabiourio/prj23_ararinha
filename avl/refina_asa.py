'''
Refina a discretizacao da ASA em estacoes de torcao, sem mexer em mais
nada do modelo.

Por que: a torcao no AVL e definida por secao (Ainc), de modo que o numero
de secoes da asa e o numero de pontos de controle da otimizacao. Com seis
secoes (cinco torcoes livres) o arrasto ja atinge o piso analitico, mas a
restricao de estol fica no limite do que a parametrizacao alcanca. Mais
secoes dao controle mais fino do carregamento local.

Como os perfis ficam: os TRES perfis otimizados no Lab 03 continuam sendo
os mestres, em airfoils/. As secoes intermediarias recebem perfis
INTERPOLADOS entre eles, gravados em airfoils/gerado/ e regerados a cada
execucao. Nao sao perfis novos nem copias mantidas a mao: sao derivados.

A mistura segue as estacoes de referencia do Lab 03:
  ate eta = 0,1011          -> raiz
  de 0,1011 a 0,398         -> mistura raiz/centro
  de 0,398 a 0,90           -> mistura centro/ponta
  de 0,90 em diante         -> ponta

Preserva do arquivo atual: cabecalho, empenagens, fuselagem, naceles,
malha, deslocamentos em Z e a posicao do aileron.

Rodar de dentro da pasta avl/:   python refina_asa.py
'''

# IMPORTS
import os
import re
import shutil

import numpy as np

#=========================================

BASE = 'aft.avl'
ORIG = os.path.join('..', 'otimizacao_aerofolio', 'ENTREGA',
                    '02_perfis_otimizados')
MESTRES = {'raiz': 'raiz_roteiro.dat', 'centro': 'meio_roteiro.dat',
           'ponta': 'ponta_roteiro.dat'}

# Estacoes de referencia dos perfis do Lab 03
ETA_RAIZ, ETA_CENTRO, ETA_PONTA = 0.1011, 0.398, 0.90

# Bordas do aileron: precisam continuar sendo bordas de secao
ETA_AIL_INI, ETA_AIL_FIM = 0.56, 0.90

# Estacoes da asa refinada. Inclui as de referencia dos perfis e as bordas
# do aileron, mais pontos intermediarios para a otimizacao de torcao.
ETAS = [0.0, 0.1011, 0.22, 0.32, 0.398, 0.48, 0.56, 0.70, 0.82, 0.90, 1.0]

N_FACE = 99          # pontos por face (limite IBX do executavel)


def reamostra(caminho):
    '''Le um .dat Selig e reamostra com espacamento cosseno por face.'''
    linhas = open(caminho).read().split('\n')
    pts = np.array([[float(v) for v in ln.split()]
                    for ln in linhas[1:] if ln.strip()])
    i_le = int(np.argmin(pts[:, 0]))
    faces = []
    for lado in (pts[:i_le + 1], pts[i_le:]):
        x, y = lado[:, 0], lado[:, 1]
        if x[0] > x[-1]:
            x, y = x[::-1], y[::-1]
        xn = x[0] + (x[-1] - x[0])*(1 - np.cos(
            np.linspace(0, np.pi, N_FACE)))/2
        faces.append(np.column_stack([xn, np.interp(xn, x, y)]))
    return np.vstack([faces[0][::-1], faces[1][1:]])


def mistura(eta, perfis):
    '''Perfil na estacao eta, misturando os mestres linearmente.'''
    if eta <= ETA_RAIZ:
        return perfis['raiz'], 'raiz'
    if eta >= ETA_PONTA:
        return perfis['ponta'], 'ponta'
    if eta <= ETA_CENTRO:
        w = (eta - ETA_RAIZ)/(ETA_CENTRO - ETA_RAIZ)
        return (1 - w)*perfis['raiz'] + w*perfis['centro'], f'raiz/centro {w:.2f}'
    w = (eta - ETA_CENTRO)/(ETA_PONTA - ETA_CENTRO)
    return (1 - w)*perfis['centro'] + w*perfis['ponta'], f'centro/ponta {w:.2f}'


def escreve_dat(caminho, curva, nome):
    with open(caminho, 'w', encoding='ascii') as f:
        f.write(nome + '\n')
        for x, y in curva:
            f.write(f'{x:.7f} {y:.7f}\n')


#=========================================
# 1. PERFIS

perfis = {k: reamostra(os.path.join(ORIG, v)) for k, v in MESTRES.items()}
for k, curva in perfis.items():
    escreve_dat(os.path.join('airfoils', f'{k}.dat'), curva, k.upper())

destino = os.path.join('airfoils', 'gerado')
shutil.rmtree(destino, ignore_errors=True)
os.makedirs(destino, exist_ok=True)

arquivo_de = {}
print('perfis por estacao:')
for eta in ETAS:
    curva, origem = mistura(eta, perfis)
    if origem in ('raiz', 'ponta'):
        arquivo_de[eta] = f'airfoils/{origem}.dat'
    else:
        nome = f'eta_{int(round(eta*1000)):04d}.dat'
        escreve_dat(os.path.join(destino, nome), curva, f'ETA {eta}')
        arquivo_de[eta] = f'airfoils/gerado/{nome}'
    print(f'  eta = {eta:.4f}  ->  {arquivo_de[eta]:34s} ({origem})')

#=========================================
# 2. GEOMETRIA DA ASA, LIDA DO ARQUIVO ATUAL
#
# A raiz e a ponta atuais ja carregam o deslocamento em Z e a geometria do
# designTool, entao as novas secoes saem por interpolacao entre elas.

linhas = open(BASE, encoding='utf-8').read().split('\n')
n = len(linhas)
inicios = [i for i in range(n) if linhas[i].strip() in ('SURFACE', 'BODY')]
i_asa = None
for i in inicios:
    if linhas[i].strip() != 'SURFACE':
        continue
    j = i + 1
    while j < n and (not linhas[j].strip() or linhas[j].strip().startswith('#')):
        j += 1
    if j < n and linhas[j].strip() == 'Wing':
        i_asa = i
if i_asa is None:
    raise SystemExit('nao achei a superficie Wing')
prox = [i for i in inicios if i > i_asa]
fim_asa = prox[0] if prox else n

# secoes atuais da asa
secoes = []
espera = False
for i in range(i_asa, fim_asa):
    t = linhas[i].split('#')[0].strip()
    if linhas[i].strip() == 'SECTION':
        espera = True
    elif espera and t:
        secoes.append([float(v) for v in t.split()[:5]])
        espera = False
raiz, ponta = np.array(secoes[0]), np.array(secoes[-1])
semi = ponta[1]
print(f'\nasa atual: {len(secoes)} secoes, semi-envergadura {semi:.4f} m')
print(f'asa nova:  {len(ETAS)} secoes ({len(ETAS)-1} torcoes livres)')

#=========================================
# 3. NOVO BLOCO DE SECOES

# tudo do bloco da asa antes da primeira SECTION (cabecalho da superficie)
i_prim = next(i for i in range(i_asa, fim_asa)
              if linhas[i].strip() == 'SECTION')
cabecalho_asa = linhas[i_asa:i_prim]

blocos = []
for eta in ETAS:
    p = raiz + eta*(ponta - raiz)          # Xle, Yle, Zle, Chord
    b = ['SECTION',
         '#Xle      Yle      Zle      Chord    Ainc',
         f'{p[0]:.4f}  {p[1]:.4f}  {p[2]:.4f}  {p[3]:.4f}  0.0000',
         'AFILE',
         arquivo_de[eta],
         'CLAF',
         '1.0000']
    if ETA_AIL_INI - 1e-9 <= eta <= ETA_AIL_FIM + 1e-9:
        b += ['CONTROL',
              '#name    gain  Xhinge  XYZhvec      SgnDup',
              'aileron  1.0   0.73    0. 0. 0.    -1.0']
    blocos.append('\n'.join(b))

novo = linhas[:i_asa] + cabecalho_asa + [''] + \
       ['\n\n'.join(blocos)] + [''] + linhas[fim_asa:]
texto = '\n'.join(novo)
open(BASE, 'w', encoding='utf-8').write(texto)

# fwd a partir do aft: so o CG muda
fwd = texto.replace('CG traseiro', 'CG dianteiro')
fwd = re.sub(r'^27\.7190 0\.0 0\.0.*$',
             '26.0772 0.0 0.0   # Xref Yref Zref (CG dianteiro)',
             fwd, flags=re.M)
open('fwd.avl', 'w', encoding='utf-8').write(fwd)

print(f'\ngravados: {BASE} e fwd.avl')
print(f'aileron nas secoes de eta {ETA_AIL_INI} a {ETA_AIL_FIM}')
