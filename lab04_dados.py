'''
Lab 04 -- dados de entrada para as analises no AVL.

Extrai do designTool tudo que o roteiro pede como insumo:
  - Tabela 1 (item 2): ponto de projeto -- o mesmo do Lab 03, peso medio
    de cruzeiro de 229.669,3 kgf, que da o CL = 0,5053 usado na otimizacao
    dos perfis;
  - posicoes de CG dianteiro/traseiro e ponto neutro, para o Xref dos
    arquivos fwd.avl e aft.avl (item 1) e para a margem estatica (item 8);
  - CD0 da polar do designTool no ponto de projeto, para o cabecalho CDp
    dos arquivos do AVL (item 1) e para a comparacao de arrasto (item 5);
  - momentos de inercia com fuel_frac consistente com o peso do ponto de
    projeto, para a tabela de derivadas do MVO (Tabelas 4 e 9).

Rodar da raiz do repo:  python lab04_dados.py
Imprime em tela e grava avl/dados_lab04.md
'''

# IMPORTS
import numpy as np

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.aerodynamics import aerodynamics
from designTool.auxiliary import atmosphere
from designTool.moment_of_inertia import moment_of_inertia
from designTool.constants import gravity

#=========================================

# Peso do ponto de projeto [kgf] -- peso medio de cruzeiro do Lab 03
# (otimizacao_aerofolio/ENTREGA/01_ponto_de_projeto/tabelas_1_e_2.md)
W_DESIGN_KGF = 229669.3

# SETUP
airplane = standard_airplane('my_airplane')
analyze(airplane, print_log=False, plot=False)

inputs = airplane['inputs']
geom = airplane['geometry']
tm = airplane['thrust_matching']
bal = airplane['balance']
ew = airplane['empty_weight']

# Condicao de voo (cruzeiro de projeto)
h = inputs['altitude_cruise']
M = inputs['Mach_cruise']
atm = atmosphere(h)
rho = atm['density']
a_inf = atm['speed_of_sound']
V = M*a_inf
q_inf = 0.5*rho*V**2

# Pesos
S_w = inputs['S_w']
W0 = tm['W0']
W_empty = tm['W_empty']
W_fuel = tm['W_fuel']
W_payload = inputs['W_payload']
W_crew = inputs['W_crew']

W = W_DESIGN_KGF*gravity
CL = W/(q_inf*S_w)

# Fracao de combustivel que reproduz o peso do ponto de projeto
fuel_frac = (W - W_empty - W_payload - W_crew)/W_fuel

# Momentos de inercia nessa condicao de carregamento (referidos ao CG)
moment_of_inertia(airplane, fuel_frac=fuel_frac, payload_frac=1.0)
moi = airplane['moment_of_inertia']

# CG do carregamento do ponto de projeto (mesma conta feita dentro de
# moment_of_inertia, que nao exporta o valor)
xcg_ponto = (W_empty*ew['xcg_empty'] + fuel_frac*W_fuel*bal['xcg_fuel']
             + W_payload*inputs['xcg_payload'] + W_crew*inputs['xcg_crew'])/W

# Polar do designTool no ponto de projeto (configuracao limpa)
CD, CLmax, dragDict = aerodynamics(airplane, Mach=M, altitude=h, CL=CL)

# Geometria de referencia e CGs
xm_w = geom['xm_w']
cm_w = geom['cm_w']
b_w = geom['b_w']
xcg_fwd = bal['xcg_fwd']
xcg_aft = bal['xcg_aft']
xnp = bal['xnp']


def pmac(x):
    '''Converte posicao longitudinal [m] para %MAC.'''
    return (x - xm_w)/cm_w*100


# RELATORIO
L = []
L.append('# Dados do Lab 04 -- gerados por lab04_dados.py')
L.append('')
L.append('Ponto de projeto: o mesmo do Lab 03 (peso medio de cruzeiro).')
L.append('')
L.append('## Tabela 1 -- ponto de projeto (item 2)')
L.append('')
L.append('| Parametro | Valor |')
L.append('|---|---|')
L.append(f'| W0 [N] | {W0:.1f} ({W0/gravity:.1f} kgf) |')
L.append(f'| W [N] | {W:.1f} ({W_DESIGN_KGF:.1f} kgf) |')
L.append(f'| h [m] | {h:.1f} |')
L.append(f'| rho_inf [kg/m3] | {rho:.5f} |')
L.append(f'| a_inf [m/s] | {a_inf:.3f} |')
L.append(f'| M | {M:.2f} |')
L.append(f'| V [m/s] | {V:.1f} |')
L.append(f'| CL | {CL:.4f} |')
L.append(f'| Sref [m2] | {S_w:.3f} |')
L.append('')
L.append('## Cabecalho dos arquivos .avl (item 1)')
L.append('')
L.append('| Parametro | Valor | Uso no .avl |')
L.append('|---|---|---|')
L.append(f'| Sref [m2] | {S_w:.3f} | linha Sref Cref Bref |')
L.append(f'| Cref = MAC [m] | {cm_w:.4f} | linha Sref Cref Bref |')
L.append(f'| Bref [m] | {b_w:.4f} | linha Sref Cref Bref |')
L.append(f'| Xref (CG dianteiro) [m] | {xcg_fwd:.4f} | fwd.avl ({pmac(xcg_fwd):.1f} %MAC) |')
L.append(f'| Xref (CG traseiro) [m] | {xcg_aft:.4f} | aft.avl ({pmac(xcg_aft):.1f} %MAC) |')
L.append(f'| CDp = CD0 designTool | {dragDict["CD0"]:.5f} | linha CDp |')
L.append('')
L.append('## CG, ponto neutro e margem estatica (designTool, item 8)')
L.append('')
L.append('| Parametro | [m] | [%MAC] |')
L.append('|---|---|---|')
L.append(f'| xcg_fwd | {xcg_fwd:.4f} | {pmac(xcg_fwd):.1f} |')
L.append(f'| xcg_aft | {xcg_aft:.4f} | {pmac(xcg_aft):.1f} |')
L.append(f'| xcg no ponto de projeto | {xcg_ponto:.4f} | {pmac(xcg_ponto):.1f} |')
L.append(f'| xnp | {xnp:.4f} | {pmac(xnp):.1f} |')
L.append('')
L.append(f"SM_fwd = {bal['SM_fwd']*100:.1f} %MAC, SM_aft = {bal['SM_aft']*100:.1f} %MAC")
L.append('(valores do designTool -- o item 8 pede os equivalentes via AVL)')
L.append('')
L.append('## Polar do designTool no ponto de projeto (itens 1 e 5)')
L.append('')
L.append('| Parametro | Valor |')
L.append('|---|---|')
L.append(f'| CD0 (parasita, config. limpa) | {dragDict["CD0"]:.5f} |')
L.append(f'| CDind | {dragDict["CDind"]:.5f} |')
L.append(f'| CDwave | {dragDict["CDwave"]:.5f} |')
L.append(f'| CD total | {dragDict["CD"]:.5f} |')
L.append(f'| K (fator de arrasto induzido) | {dragDict["K"]:.5f} |')
L.append(f'| e (Oswald) | {dragDict["e"]:.4f} |')
L.append(f'| CLmax limpa (estimativa designTool) | {dragDict["CLmax_clean"]:.3f} |')
L.append('')
L.append('## Momentos de inercia -- Tabelas 4 e 9 do MVO')
L.append('')
L.append(f'Carregamento: fuel_frac = {fuel_frac:.4f} (reproduz W = {W_DESIGN_KGF:.1f} kgf), '
         'payload_frac = 1.0. Referencia: CG deste carregamento.')
L.append('')
L.append('| Parametro | Valor |')
L.append('|---|---|')
L.append(f'| m [kg] | {W/gravity:.1f} |')
L.append(f'| Ixx [kg m2] | {moi["Ixx"]:.4e} |')
L.append(f'| Iyy [kg m2] | {moi["Iyy"]:.4e} |')
L.append(f'| Izz [kg m2] | {moi["Izz"]:.4e} |')
L.append(f'| Ixz [kg m2] | {moi["Ixz"]:.4e} |')
L.append(f'| Ixy [kg m2] | {moi["Ixy"]:.4e} |')
L.append(f'| Iyz [kg m2] | {moi["Iyz"]:.4e} |')
L.append('')
L.append('## Demais itens da Tabela 4 (motor)')
L.append('')
L.append('| Parametro | Valor | Observacao |')
L.append('|---|---|---|')
L.append(f'| ip [deg] | 0.0 | designTool nao define incidencia de motor |')
L.append(f"| xp [m] | {inputs['x_n']:.3f} | face frontal da nacele (x_n) |")
L.append(f"| zp [m] | {inputs['z_n']:.3f} | inverter o sinal para o MVO |")
L.append(f"| Tmax por motor [N] | {inputs['engine']['Tmax']:.0f} | n_engines = {inputs['n_engines']} |")
L.append(f"| T0 (thrust matching) [N] | {tm['T0']:.0f} | tracao total de decolagem |")
L.append('')

texto = '\n'.join(L)
print(texto)

with open('avl/dados_lab04.md', 'w', encoding='utf-8') as f:
    f.write(texto + '\n')

print('Arquivo gravado: avl/dados_lab04.md')
