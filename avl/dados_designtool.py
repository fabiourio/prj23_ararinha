'''
Ponte unica entre o designTool e o Lab 04.

Faz tres coisas:
  1. calcula o PONTO DE PROJETO (Tabela 1 do roteiro);
  2. reune a geometria e os dados de referencia da aeronave;
  3. VERIFICA os arquivos fwd.avl e aft.avl contra o designTool, rodando o
     proprio AVL para ler as areas que ele integra sobre os paineis.

Grava dados_designtool.json, que os scripts Julia leem. Assim o designTool
e chamado em um unico lugar.

Rodar de dentro da pasta avl/:   python dados_designtool.py
'''

# IMPORTS
import json
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.abspath('..'))

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.aerodynamics import aerodynamics
from designTool.auxiliary import atmosphere
from designTool.moment_of_inertia import moment_of_inertia
from designTool.constants import gravity

#=========================================

# Peso do ponto de projeto [kgf]. E o peso medio de cruzeiro usado no
# Lab 03, o mesmo para o qual os perfis foram otimizados, de modo que o
# CL de projeto da aeronave coincide com o CL de projeto das secoes.
W_DESIGN_KGF = 229669.3

# Deslocamentos em Z aplicados no modelo do AVL para tirar as superficies
# de dentro da fuselagem (sao intencionais e verificados aqui).
DZ = {'asa': -1.20, 'EH': +1.85, 'EV': +0.85}

TOL_AREA = 0.01          # 1% nas areas
TOL_GEOM = 0.002         # 2 mm nas dimensoes lineares


def barra(titulo):
    print('\n' + '=' * 74)
    print(titulo)
    print('=' * 74)


def confere(nome, avl, ref, tol, unidade='', rel=False):
    '''Imprime uma linha de comparacao e devolve se passou.'''
    d = abs(avl - ref)/abs(ref) if rel else abs(avl - ref)
    ok = d <= tol
    marca = 'ok' if ok else 'DIVERGE'
    if rel:
        print(f'  {nome:34s} {ref:12.4f} {avl:12.4f}   {100*d:6.2f}%   {marca}')
    else:
        print(f'  {nome:34s} {ref:12.4f} {avl:12.4f}   {d:7.4f}{unidade:2s}  {marca}')
    return ok


#=========================================
# 1. AERONAVE E PONTO DE PROJETO

airplane = standard_airplane('my_airplane')
analyze(airplane, print_log=False, plot=False)

inp = airplane['inputs']
geo = airplane['geometry']
tm = airplane['thrust_matching']
bal = airplane['balance']
ew = airplane['empty_weight']

h = inp['altitude_cruise']
M = inp['Mach_cruise']
atm = atmosphere(h)
rho, a_inf = atm['density'], atm['speed_of_sound']
V = M*a_inf
q = 0.5*rho*V**2
S_w = inp['S_w']

W = W_DESIGN_KGF*gravity
CL_proj = W/(q*S_w)

# fracao de combustivel que reproduz esse peso, para as inercias
fuel_frac = (W - tm['W_empty'] - inp['W_payload'] - inp['W_crew'])/tm['W_fuel']
moment_of_inertia(airplane, fuel_frac=fuel_frac, payload_frac=1.0)
moi = airplane['moment_of_inertia']

# polar do designTool no ponto de projeto
CD_dt, CLmax_dt, drag = aerodynamics(airplane, Mach=M, altitude=h, CL=CL_proj)

barra('PONTO DE PROJETO (Tabela 1 do roteiro)')
print(f'  {"W0 (peso maximo de decolagem)":34s} {tm["W0"]:12.1f} N'
      f'   ({tm["W0"]/gravity:.1f} kgf)')
print(f'  {"W (peso no ponto de projeto)":34s} {W:12.1f} N'
      f'   ({W_DESIGN_KGF:.1f} kgf)')
print(f'  {"h (altitude)":34s} {h:12.1f} m')
print(f'  {"rho_inf":34s} {rho:12.5f} kg/m3')
print(f'  {"a_inf":34s} {a_inf:12.3f} m/s')
print(f'  {"M":34s} {M:12.2f}')
print(f'  {"V":34s} {V:12.1f} m/s')
print(f'  {"CL de projeto":34s} {CL_proj:12.4f}')
print(f'  {"Sref":34s} {S_w:12.3f} m2')
print(f'\n  fracao de combustivel equivalente: {fuel_frac:.4f} '
      f'(com 100% de carga paga)')
print(f'  polar do designTool nesse ponto: CD0 = {drag["CD0"]:.5f}, '
      f'CDind = {drag["CDind"]:.5f}, CD = {CD_dt:.5f}')

#=========================================
# 2. GEOMETRIA DE REFERENCIA

ref = {
    'Sref': S_w, 'Cref': geo['cm_w'], 'Bref': geo['b_w'],
    'asa': {'b': geo['b_w'], 'cr': geo['cr_w'], 'ct': geo['ct_w'],
            'xr': inp['xr_w'], 'zr': inp['zr_w'],
            'xt': geo['xt_w'], 'yt': geo['yt_w'], 'zt': geo['zt_w'],
            'S_proj': S_w},
    'EH': {'S': geo['S_h'], 'b': geo['b_h'], 'cr': geo['cr_h'],
           'ct': geo['ct_h'], 'xr': geo['xr_h'], 'zr': inp['zr_h'],
           'xt': geo['xt_h'], 'yt': geo['yt_h'], 'zt': geo['zt_h']},
    'EV': {'S': geo['S_v'], 'b': geo['b_v'], 'cr': geo['cr_v'],
           'ct': geo['ct_v'], 'xr': geo['xr_v'], 'zr': inp['zr_v'],
           'xt': geo['xt_v'], 'zt': geo['zt_v']},
    'nacele': {'L': inp['L_n'], 'D': inp['D_n'], 'x': inp['x_n'],
               'y': inp['y_n'], 'z': inp['z_n']},
    'fuselagem': {'L': inp['L_f'], 'D': inp['D_f']},
    'cg': {'fwd': bal['xcg_fwd'], 'aft': bal['xcg_aft'], 'np': bal['xnp']},
}

#=========================================
# 3. VERIFICACAO DOS ARQUIVOS DO AVL


def le_avl(arquivo):
    '''Le cabecalho, secoes por superficie e blocos de corpo/nacele.'''
    linhas = open(arquivo, encoding='utf-8').read().split('\n')
    dados = {'secoes': {}, 'malha': {}}
    # cabecalho: primeiras linhas nao comentadas
    vals = []
    for ln in linhas[1:]:
        t = ln.split('#')[0].strip()
        if t:
            vals.append(t)
        if len(vals) == 5:
            break
    dados['mach'] = float(vals[0])
    dados['sref'], dados['cref'], dados['bref'] = map(float, vals[2].split())
    dados['xref'] = float(vals[3].split()[0])
    dados['cdp'] = float(vals[4])

    sup, esp_nome, esp_malha, esp_sec = '', False, False, False
    for ln in linhas:
        t = ln.split('#')[0].strip()
        bruto = ln.strip()
        if bruto == 'SURFACE':
            esp_nome = True
        elif esp_nome and bruto and not bruto.startswith('#'):
            sup = bruto
            dados['secoes'].setdefault(sup, [])
            esp_nome, esp_malha = False, True
        elif esp_malha and t:
            p = t.split()
            dados['malha'][sup] = (int(p[0]), int(p[2]))
            esp_malha = False
        elif bruto == 'SECTION':
            esp_sec = True
        elif esp_sec and t:
            dados['secoes'][sup].append([float(v) for v in t.split()[:5]])
            esp_sec = False
        elif bruto == 'SCALE' and sup == 'Nacelle':
            esp_sec = False
    # SCALE e TRANSLATE da nacele
    txt = '\n'.join(linhas)
    bloco = txt.split('Nacelle')[-1] if 'Nacelle' in txt else ''
    for chave in ('SCALE', 'TRANSLATE'):
        if chave in bloco:
            for ln in bloco.split(chave)[1].split('\n')[1:]:
                t = ln.split('#')[0].strip()
                if t:
                    dados[chave.lower()] = [float(v) for v in t.split()]
                    break
    return dados


def areas_do_avl(arquivo):
    '''Areas que o proprio AVL integra sobre os paineis.'''
    cmds = f'load {arquivo}\noper\nx\nfn\n\n\nquit\n'
    with open('_v.txt', 'w') as f:
        f.write(cmds)
    with open('_v.txt') as f:
        out = subprocess.run(['./avl.exe'], stdin=f, capture_output=True,
                             text=True).stdout
    os.remove('_v.txt')
    areas = {}
    for m in re.finditer(r'^\s*\d+\s+([\d.]+)\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+'
                         r'\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+'
                         r'\s+(.+?)\s*$', out, re.M):
        nome = m.group(2).replace('(YDUP)', '').strip()
        areas[nome] = areas.get(nome, 0.0) + float(m.group(1))
    return areas


avl = le_avl('aft.avl')
areas = areas_do_avl('aft.avl')
tudo_ok = True

barra('VERIFICACAO: CABECALHO')
print(f'  {"grandeza":34s} {"designTool":>12s} {"AVL":>12s}   desvio')
tudo_ok &= confere('Sref [m2]', avl['sref'], ref['Sref'], TOL_GEOM, ' m2')
tudo_ok &= confere('Cref (MAC) [m]', avl['cref'], ref['Cref'], TOL_GEOM, ' m')
tudo_ok &= confere('Bref (envergadura) [m]', avl['bref'], ref['Bref'], TOL_GEOM, ' m')
tudo_ok &= confere('CDp (= CD0 do designTool)', avl['cdp'], drag['CD0'],
                   1e-4, '')
tudo_ok &= confere('Xref (CG traseiro) [m]', avl['xref'], ref['cg']['aft'],
                   TOL_GEOM, ' m')
xref_fwd = le_avl('fwd.avl')['xref']
tudo_ok &= confere('Xref (CG dianteiro) [m]', xref_fwd, ref['cg']['fwd'],
                   TOL_GEOM, ' m')

barra('VERIFICACAO: GEOMETRIA DA ASA')
sec = avl['secoes']['Wing']
raiz, ponta = sec[0], sec[-1]
print(f'  {"grandeza":34s} {"designTool":>12s} {"AVL":>12s}   desvio')
tudo_ok &= confere('envergadura (2 x Yle ponta) [m]', 2*ponta[1],
                   ref['asa']['b'], TOL_GEOM, ' m')
tudo_ok &= confere('corda de raiz [m]', raiz[3], ref['asa']['cr'], TOL_GEOM, ' m')
tudo_ok &= confere('corda de ponta [m]', ponta[3], ref['asa']['ct'], TOL_GEOM, ' m')
tudo_ok &= confere('Xle da raiz [m]', raiz[0], ref['asa']['xr'], TOL_GEOM, ' m')
tudo_ok &= confere('Xle da ponta [m]', ponta[0], ref['asa']['xt'], TOL_GEOM, ' m')
# area projetada do trapezio (o AVL integra a area real, com diedro)
S_trap = (raiz[3] + ponta[3])*ponta[1]
tudo_ok &= confere('area projetada [m2]', S_trap, ref['asa']['S_proj'],
                   TOL_AREA, rel=True)
print(f'  {"area integrada pelo AVL (diedro)":34s} {"":12s} '
      f'{areas.get("Wing", 0):12.3f}   '
      f'{100*(areas.get("Wing", 0)/ref["asa"]["S_proj"] - 1):+6.2f}% (esperado)')
# deslocamento em Z intencional
dz_asa = raiz[2] - ref['asa']['zr']
tudo_ok &= confere('deslocamento em Z (intencional) [m]', dz_asa, DZ['asa'],
                   TOL_GEOM, ' m')

barra('VERIFICACAO: EMPENAGENS')
print(f'  {"grandeza":34s} {"designTool":>12s} {"AVL":>12s}   desvio')
sh = avl['secoes']['Horizontal tail']
tudo_ok &= confere('EH: area [m2]', areas.get('Horizontal tail', 0),
                   ref['EH']['S'], TOL_AREA, rel=True)
tudo_ok &= confere('EH: envergadura [m]', 2*sh[-1][1], ref['EH']['b'],
                   TOL_GEOM, ' m')
tudo_ok &= confere('EH: corda de raiz [m]', sh[0][3], ref['EH']['cr'],
                   TOL_GEOM, ' m')
tudo_ok &= confere('EH: deslocamento em Z [m]', sh[0][2] - ref['EH']['zr'],
                   DZ['EH'], TOL_GEOM, ' m')
sv = avl['secoes']['Vertical tail']
tudo_ok &= confere('EV: area [m2]', areas.get('Vertical tail', 0),
                   ref['EV']['S'], TOL_AREA, rel=True)
tudo_ok &= confere('EV: altura [m]', sv[-1][2] - sv[0][2], ref['EV']['b'],
                   TOL_GEOM, ' m')
tudo_ok &= confere('EV: corda de raiz [m]', sv[0][3], ref['EV']['cr'],
                   TOL_GEOM, ' m')
tudo_ok &= confere('EV: deslocamento em Z [m]', sv[0][2] - ref['EV']['zr'],
                   DZ['EV'], TOL_GEOM, ' m')

barra('VERIFICACAO: NACELE E FUSELAGEM')
print(f'  {"grandeza":34s} {"designTool":>12s} {"AVL":>12s}   desvio')
esc, tra = avl.get('scale', [0, 0, 0]), avl.get('translate', [0, 0, 0])
tudo_ok &= confere('nacele: comprimento [m]', esc[0], ref['nacele']['L'],
                   TOL_GEOM, ' m')
tudo_ok &= confere('nacele: diametro [m]', 2*esc[1], ref['nacele']['D'],
                   TOL_GEOM, ' m')
tudo_ok &= confere('nacele: x da face frontal [m]', tra[0],
                   ref['nacele']['x'], TOL_GEOM, ' m')
tudo_ok &= confere('nacele: y [m]', tra[1], ref['nacele']['y'], TOL_GEOM, ' m')
tudo_ok &= confere('nacele: z (com desloc. da asa) [m]', tra[2],
                   ref['nacele']['z'] + DZ['asa'], TOL_GEOM, ' m')

barra('RESULTADO')
if tudo_ok:
    print('  O modelo do AVL reproduz a geometria do designTool.')
else:
    print('  ATENCAO: ha divergencias marcadas acima.')

#=========================================
# 4. ARQUIVO PARA OS SCRIPTS JULIA

saida = {
    'ponto_de_projeto': {
        'W0_N': tm['W0'], 'W_N': W, 'W_kgf': W_DESIGN_KGF,
        'h_m': h, 'rho': rho, 'a_inf': a_inf, 'M': M, 'V': V,
        'CL': CL_proj, 'Sref': S_w, 'q_inf': q,
        'fuel_frac': fuel_frac, 'payload_frac': 1.0,
    },
    'referencia': {'Sref': ref['Sref'], 'Cref': ref['Cref'],
                   'Bref': ref['Bref']},
    'cg': ref['cg'],
    'polar_designtool': {'CD0': drag['CD0'], 'CDind': drag['CDind'],
                         'CDwave': drag['CDwave'], 'CD': CD_dt,
                         'K': drag['K'], 'e': drag['e'],
                         'CLmax_limpa': drag['CLmax_clean']},
    'inercias': {'m_kg': W/gravity, 'Ixx': moi['Ixx'], 'Iyy': moi['Iyy'],
                 'Izz': moi['Izz'], 'Ixz': moi['Ixz']},
    'motor': {'x_n': inp['x_n'], 'y_n': inp['y_n'], 'z_n': inp['z_n'],
              'Tmax_por_motor': inp['engine']['Tmax'],
              'n_motores': inp['n_engines'], 'T0': tm['T0']},
    'geometria': {'b_w': geo['b_w'], 'cr_w': geo['cr_w'], 'ct_w': geo['ct_w'],
                  'xm_w': geo['xm_w'], 'cm_w': geo['cm_w'],
                  'yt_w': geo['yt_w']},
    'dz_avl': DZ,
}
with open('dados_designtool.json', 'w', encoding='ascii') as f:
    json.dump(saida, f, indent=2)
print('\ngravado: dados_designtool.json')
