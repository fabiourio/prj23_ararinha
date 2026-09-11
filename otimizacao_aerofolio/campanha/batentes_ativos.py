'''
Quais batentes estao ativos no otimo? -- verificacao de otimo interior.

Um otimo encostado num batente nao e o otimo do PROBLEMA, e o otimo da
CAIXA. A condicao de otimalidade muda: em vez de gradiente nulo, ha
multiplicador do batente. Se o batente for arbitrario -- e os do roteiro sao
sugestoes, nao fisica -- o resultado esta limitado por uma escolha nossa.

O Lab 02 fez essa mesma verificacao para a aeronave ("o otimo e interior?").
Aqui a resposta importa ainda mais, porque os batentes do CST nao tem
significado fisico obvio.

Rodar com:  python batentes_ativos.py
'''

import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import otimiza_secao as ot
from analisa_otimos import carrega

TOL = 1e-4


def analisa(nome, pasta):
    d = carrega(pasta)
    if d is None:
        return None
    Al, Au, alpha = ot.desmonta(d['xx_otimo'])
    lin = []
    for i, v in enumerate(Al):
        if abs(v - ot.AL_UPPER[i]) < TOL:
            lin.append((f'Al{i+1}', v, 'superior', ot.AL_UPPER[i]))
        elif abs(v - ot.AL_LOWER[i]) < TOL:
            lin.append((f'Al{i+1}', v, 'inferior', ot.AL_LOWER[i]))
    for i, v in enumerate(Au):
        if abs(v - ot.AU_LOWER[i]) < TOL:
            lin.append((f'Au{i+1}', v, 'inferior', ot.AU_LOWER[i]))
        elif abs(v - ot.AU_UPPER[i]) < TOL:
            lin.append((f'Au{i+1}', v, 'superior', ot.AU_UPPER[i]))
    if abs(alpha - ot.ALPHA_MIN) < TOL or abs(alpha - ot.ALPHA_MAX) < TOL:
        lin.append(('alpha', alpha * 180 / np.pi, 'batente', np.nan))
    return d, Al, Au, alpha, lin


def main():
    print('Batentes do roteiro:  Al em [-1,00, -0,05]   Au em [+0,05, +1,00]')
    print(f'Tolerancia de deteccao: {TOL}\n')

    total = 0
    for nome in ot.ESTACOES:
        for sufixo, rot in [('', 'batentes do roteiro'),
                            ('_cusp', 'com Al4 solto')]:
            r = analisa(nome, f'otim_{nome}{sufixo}')
            if r is None:
                continue
            d, Al, Au, alpha, lin = r
            cd = d['hist']['CD'][d['i_otimo']]
            print(f'--- {nome} ({rot}) --- c_d = {cd:.6f}')
            print(f'    Al = {np.array2string(Al, precision=4)}')
            print(f'    Au = {np.array2string(Au, precision=4)}')
            if lin:
                total += len(lin)
                for var, v, lado, lim in lin:
                    print(f'    >> {var} = {v:+.5f} NO BATENTE {lado} '
                          f'({lim:+.3f})')
            else:
                print('    (nenhum batente ativo -- otimo interior)')
            print()

    print('=' * 70)
    if total:
        print(f'{total} batente(s) ativo(s) no total.')
        print('\nConsequencia: estes otimos NAO sao interiores. A condicao de')
        print('otimalidade tem multiplicador de batente, e o valor de c_d')
        print('esta limitado por uma escolha de caixa, nao pela fisica.')
        print('\nOs batentes do roteiro nao sao fisica: Au >= 0,05 e')
        print('Al <= -0,05 garantiam um raio de bordo de ataque minimo e')
        print('impediam as superficies de se cruzarem no bordo. Mas nos ja')
        print('temos restricoes EXPLICITAS para as duas coisas -- a bluntez')
        print('(substituta de cl_max) e a espessura minima -- e elas fazem o')
        print('servico melhor, porque sao calibradas em vez de arbitradas.')
        print('\nEntao vale alargar os batentes e deixar as restricoes')
        print('fisicas trabalharem. E o que a rodada --batentes-largos faz.')
    else:
        print('Nenhum batente ativo: os otimos sao interiores.')


if __name__ == '__main__':
    main()
