'''
Lab 04 -- gera os arquivos de entrada do AVL a partir do designTool.

Produz em avl/:
  - fwd.avl e aft.avl (item 1 do roteiro): mesma geometria, mudando apenas
    o Xref (CG dianteiro ou traseiro do designTool). CDp recebe o CD0 do
    designTool no ponto de projeto e as asas usam os perfis otimizados do
    Lab 03 (familia roteiro, ja copiados em avl/aerofolios/).
  - fwd.mass e aft.mass: massa, CG e inercias do ponto de projeto do Lab 03
    (peso medio de cruzeiro, fuel_frac ~0,44, 100%% de carga paga). So sao
    necessarios para analise de modos (comando mode), nao para o roteiro.

Escolhas de modelagem (nao vem do designTool):
  - arrasto viscoso apenas pelo CDp do cabecalho, que recebe o CD0 por
    areas molhadas do designTool (orientacao do professor). Sem CDCL por
    secao: o CD0 ja cobre o arrasto parasita da asa, e o metodo da secao
    critica e as polares do roteiro nao dependem dele;
  - CLAF = 1,0 em todas as secoes por enquanto (sem a correcao de
    espessura 1 + 0,77 t/c do manual do AVL);
  - fuselagem nao modelada (VLM classico, como o b737simple.avl da aula);
  - asa com secoes nas estacoes do Lab 03: eta 0,1011 (raiz), 0,398 (meio)
    e 0,90 (ponta), com o perfil da raiz estendido ate o plano de simetria
    e o da ponta ate eta 1,0;
  - aileron com 27%% da corda (c_ail_c_wing) de eta 0,56 a 0,90, fechando
    b_ail_b_wing = 0,34 da semi-envergadura na estacao da ponta;
  - profundor e leme com charneira em 70%% da corda (valor tipico; o
    designTool nao define);
  - empenagens com NACA 0010 (t/c = 0,10 dos inputs; simetrico);
  - winglet nao modelado.

Numeracao que o AVL vai atribuir (ordem de aparicao no arquivo):
  controles: d1 = aileron, d2 = elevator, d3 = rudder
  variaveis de projeto: 1 = it
  Nos comandos do roteiro, troque "d4 pm 0" por "d2 pm 0" e, no menu de,
  "2 <valor>" por "1 <valor>".

Rodar da raiz do repo:  python lab04_gera_avl.py
Depois rode o AVL de dentro da pasta avl/ (os AFILE sao relativos a ela).
'''

# IMPORTS
import os

import numpy as np

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.aerodynamics import aerodynamics
from designTool.auxiliary import atmosphere
from designTool.moment_of_inertia import moment_of_inertia
from designTool.constants import gravity

#=========================================

# Ponto de projeto (o mesmo de lab04_dados.py)
W_DESIGN_KGF = 229669.3

# SETUP
airplane = standard_airplane('my_airplane')
analyze(airplane, print_log=False, plot=False)

inputs = airplane['inputs']
geom = airplane['geometry']
tm = airplane['thrust_matching']
bal = airplane['balance']
ew = airplane['empty_weight']

h = inputs['altitude_cruise']
M = inputs['Mach_cruise']
atm = atmosphere(h)
rho = atm['density']
V = M*atm['speed_of_sound']
q_inf = 0.5*rho*V**2
S_w = inputs['S_w']

W = W_DESIGN_KGF*gravity
CL = W/(q_inf*S_w)
fuel_frac = (W - tm['W_empty'] - inputs['W_payload'] - inputs['W_crew'])/tm['W_fuel']

# CD0 para o cabecalho CDp
_, _, dragDict = aerodynamics(airplane, Mach=M, altitude=h, CL=CL)
CD0 = dragDict['CD0']

# Inercias e CG do carregamento do ponto de projeto
moment_of_inertia(airplane, fuel_frac=fuel_frac, payload_frac=1.0)
moi = airplane['moment_of_inertia']
xcg_ponto = (tm['W_empty']*ew['xcg_empty'] + fuel_frac*tm['W_fuel']*bal['xcg_fuel']
             + inputs['W_payload']*inputs['xcg_payload']
             + inputs['W_crew']*inputs['xcg_crew'])/W
massa = W/gravity

#=========================================
# PERFIS DA ASA

# O avl.exe da disciplina limita o numero de pontos de aerofolio (IBX < 321,
# e os .dat do Lab 03 tem 321). Reamostramos para 99 pontos por face com
# espacamento cosseno (o AVL so usa a linha de arqueamento, entao nao ha
# perda pratica de fidelidade).
ORIGEM_PERFIS = 'otimizacao_aerofolio/ENTREGA/02_perfis_otimizados'
PERFIS = {'raiz.dat': 'raiz_roteiro.dat',
          'meio.dat': 'meio_roteiro.dat',
          'ponta.dat': 'ponta_roteiro.dat'}


def reamostra_dat(origem, destino, n_face=99):
    '''Reamostra um .dat Selig (TE -> LE -> TE) com cosseno em x por face.'''
    with open(origem) as f:
        linhas = f.read().split('\n')
    nome = linhas[0].strip()
    pts = np.array([[float(v) for v in ln.split()]
                    for ln in linhas[1:] if ln.strip()])
    i_le = int(np.argmin(pts[:, 0]))
    faces = []
    for lado in (pts[:i_le + 1], pts[i_le:]):
        x, y = lado[:, 0], lado[:, 1]
        if x[0] > x[-1]:                      # np.interp exige x crescente
            x, y = x[::-1], y[::-1]
        x_novo = x[0] + (x[-1] - x[0])*(1 - np.cos(
            np.linspace(0, np.pi, n_face)))/2
        faces.append(np.column_stack([x_novo, np.interp(x_novo, x, y)]))
    # Remonta no sentido original: TE -> LE (face 1) + LE -> TE (face 2)
    face1 = faces[0][::-1]
    curva = np.vstack([face1, faces[1][1:]])
    with open(destino, 'w', encoding='ascii') as f:
        f.write(nome + '\n')
        for xx, yy in curva:
            f.write(f'{xx:.7f} {yy:.7f}\n')


for destino, origem in PERFIS.items():
    reamostra_dat(os.path.join(ORIGEM_PERFIS, origem),
                  os.path.join('avl/aerofolios', destino))
    print(f'avl/aerofolios/{destino} reamostrado ({2*99 - 1} pontos)')

#=========================================
# GEOMETRIA DAS SUPERFICIES

# Asa: bordo de ataque interpolado linearmente entre a raiz (plano de
# simetria) e a ponta, como no trapezio do designTool
wing_root = np.array([inputs['xr_w'], 0.0, inputs['zr_w']])
wing_tip = np.array([geom['xt_w'], geom['yt_w'], geom['zt_w']])
cr_w, ct_w = geom['cr_w'], geom['ct_w']

ht_root = np.array([geom['xr_h'], 0.0, inputs['zr_h']])
ht_tip = np.array([geom['xt_h'], geom['yt_h'], geom['zt_h']])
cr_h, ct_h = geom['cr_h'], geom['ct_h']

vt_root = np.array([geom['xr_v'], 0.0, inputs['zr_v']])
vt_tip = np.array([geom['xt_v'], 0.0, geom['zt_v']])
cr_v, ct_v = geom['cr_v'], geom['ct_v']

def secao_asa(eta, estacao, extras=()):
    '''Bloco SECTION da asa na estacao adimensional eta.'''
    le = wing_root + eta*(wing_tip - wing_root)
    c = cr_w + eta*(ct_w - cr_w)
    linhas = ['SECTION',
              '#Xle      Yle      Zle      Chord    Ainc',
              f'{le[0]:.4f}  {le[1]:.4f}  {le[2]:.4f}  {c:.4f}  0.0000',
              'AFILE',
              f'aerofolios/{estacao}.dat',
              'CLAF',
              '1.0000']
    linhas += list(extras)
    return '\n'.join(linhas)


def secao_emp(root, tip, cr, ct, frac, extras=()):
    '''Bloco SECTION de empenagem (frac = 0 na raiz, 1 na ponta).'''
    le = root + frac*(tip - root)
    c = cr + frac*(ct - cr)
    linhas = ['SECTION',
              '#Xle      Yle      Zle      Chord    Ainc',
              f'{le[0]:.4f}  {le[1]:.4f}  {le[2]:.4f}  {c:.4f}  0.0000',
              'NACA',
              '0010',
              'CLAF',
              '1.0000']
    linhas += list(extras)
    return '\n'.join(linhas)


AILERON = ('CONTROL',
           '#name    gain  Xhinge  XYZhvec      SgnDup',
           'aileron  1.0   0.73    0. 0. 0.    -1.0')
ELEVATOR = ('DESIGN',
            '#DName  Wdes',
            'it      1.0',
            'CONTROL',
            '#name    gain  Xhinge  XYZhvec      SgnDup',
            'elevator 1.0   0.70    0. 0. 0.     1.0')
RUDDER = ('CONTROL',
          '#name    gain  Xhinge  XYZhvec      SgnDup',
          'rudder   1.0   0.70    0. 0. 0.     0.0')


def monta_avl(nome_cg, xref):
    '''Monta o texto completo de um arquivo .avl para o CG dado.'''
    partes = []
    partes.append(f'''Ararinha
#
# Lab 04 -- gerado por lab04_gera_avl.py a partir do designTool
# CG {nome_cg} do designTool. Perfis otimizados do Lab 03 (roteiro).
# Controles: d1 = aileron, d2 = elevator, d3 = rudder. Design: 1 = it.
#
0.85          # Mach
0 0 0.0       # iYsym iZsym Zsym
{S_w:.3f} {geom['cm_w']:.4f} {geom['b_w']:.4f}   # Sref Cref Bref
{xref:.4f} 0.0 0.0   # Xref Yref Zref (CG {nome_cg})
{CD0:.5f}       # CDp (CD0 do designTool no ponto de projeto)
#----------------------------------------------------------------
SURFACE
Wing
#Nchordwise Cspace Nspanwise Sspace
12 1.0 40 1.0
YDUPLICATE
0.0
ANGLE
0.0''')
    partes.append(secao_asa(0.0, 'raiz'))
    partes.append(secao_asa(0.1011, 'raiz'))
    partes.append(secao_asa(0.398, 'meio'))
    partes.append(secao_asa(0.56, 'meio', AILERON))
    partes.append(secao_asa(0.90, 'ponta', AILERON))
    partes.append(secao_asa(1.0, 'ponta'))
    partes.append('''#----------------------------------------------------------------
SURFACE
Horizontal tail
#Nchordwise Cspace Nspanwise Sspace
8 1.0 16 1.0
YDUPLICATE
0.0
ANGLE
0.0''')
    partes.append(secao_emp(ht_root, ht_tip, cr_h, ct_h, 0.0, ELEVATOR))
    partes.append(secao_emp(ht_root, ht_tip, cr_h, ct_h, 1.0, ELEVATOR))
    partes.append('''#----------------------------------------------------------------
SURFACE
Vertical tail
#Nchordwise Cspace Nspanwise Sspace
8 1.0 12 1.0
ANGLE
0.0''')
    partes.append(secao_emp(vt_root, vt_tip, cr_v, ct_v, 0.0, RUDDER))
    partes.append(secao_emp(vt_root, vt_tip, cr_v, ct_v, 1.0, RUDDER))
    return '\n\n'.join(partes) + '\n'


def monta_mass(nome_cg, xref):
    '''Arquivo .mass: aeronave inteira como um unico item, com as inercias
    do ponto de projeto. As inercias foram calculadas em torno do CG do
    carregamento do ponto de projeto (x = %.4f m); usar o mesmo valor nos
    dois CGs e uma aproximacao (diferenca da ordem de 0,5%% em Iyy).''' % xcg_ponto
    return f'''#  Ararinha -- CG {nome_cg} (Lab 04, ponto de projeto do Lab 03)
#  Massa e inercias do designTool (moment_of_inertia, fuel_frac = {fuel_frac:.4f},
#  payload_frac = 1.0). Inercias calculadas em torno do CG do carregamento
#  (x = {xcg_ponto:.4f} m); mantidas iguais nos dois arquivos por aproximacao.
Lunit = 1.0 m
Munit = 1.0 kg
Tunit = 1.0 s
g   = {gravity}
rho = {rho:.5f}
#  mass        x         y      z      Ixx         Iyy         Izz         Ixy    Ixz         Iyz
{massa:.1f}  {xref:.4f}  0.0  0.0  {moi['Ixx']:.4e}  {moi['Iyy']:.4e}  {moi['Izz']:.4e}  0.0  {moi['Ixz']:.4e}  0.0
'''


# EXECUTION
cgs = {'dianteiro': ('fwd', bal['xcg_fwd']),
       'traseiro': ('aft', bal['xcg_aft'])}

for nome_cg, (prefixo, xref) in cgs.items():
    with open(f'avl/{prefixo}.avl', 'w', encoding='ascii') as f:
        f.write(monta_avl(nome_cg, xref))
    with open(f'avl/{prefixo}.mass', 'w', encoding='ascii') as f:
        f.write(monta_mass(nome_cg, xref))
    print(f'avl/{prefixo}.avl e avl/{prefixo}.mass gravados '
          f'(Xref = {xref:.4f} m)')

print(f'\nCDp (CD0 designTool) = {CD0:.5f}')
print(f'CL de projeto = {CL:.4f} (meta dos comandos "a c" no AVL)')
print('Controles: d1 = aileron, d2 = elevator, d3 = rudder; design 1 = it')
print('Rode o AVL de dentro de avl/ para os AFILE serem encontrados.')
