'''
Quanto o NOSSO XFoil erra em cl_max? -- afericao contra dados experimentais.

Motivo: o alvo clmax_w = 1,80 do designTool e um valor de HANDBOOK, e o
cl_max que medimos vem do XFoil a Re alto com Ncrit = 9 (tunel limpo, sem
rugosidade de fabricacao, sem efeito 3D). Se a escala do XFoil estiver
deslocada, comparar os dois numeros e comparar regua com metro -- e a
calibracao do limiar da restricao fica errada na mesma medida.

Metodo: rodar perfis NACA com cl_max experimental publicado (Abbott & von
Doenhoff, "Theory of Wing Sections", dados da NACA TR-824, tunel de
densidade variavel) pela nossa configuracao exata, no mesmo Reynolds dos
ensaios, e medir o desvio.

Rodar com:  python afere_xfoil.py
'''

import os
import sys
import tempfile

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import airfoil_mod as am
import xfoil_runner as xr

# cl_max experimental de Abbott & von Doenhoff, Apendice IV (TR-824),
# perfis lisos, a Re = 8,9e6 -- o maior Reynolds tabelado.
REFERENCIA = [
    ('0009', 8.9e6, 1.40),
    ('0012', 8.9e6, 1.60),
    ('2412', 8.9e6, 1.70),
    ('4412', 8.9e6, 1.85),
    ('23012', 8.9e6, 1.80),
]

ALPHA_SEQ = (0.0, 30.0, 0.5)


def roda_naca(naca, Re, n_crit=9.0):
    '''
    Usa o gerador NACA NATIVO do XFoil, em vez de alimentar coordenadas.

    Alimentar as coordenadas do am.nacafoil nao convergia em nenhum dos cinco
    perfis, provavelmente por causa do bordo de fuga rombudo que o nacafoil
    gera por padrao (sharp_te=False). O comando NACA do proprio XFoil ja
    estava comprovado no diagnostico inicial da interface, entao usamos ele:
    tira uma etapa do caminho e uma fonte de erro junto.
    '''
    import shutil
    import subprocess

    tmp = tempfile.mkdtemp(prefix='afere_')
    try:
        a0, a1, da = ALPHA_SEQ
        cmds = ['PLOP', 'G F', '', f'NACA {naca}',
                'PPAR', 'T', f'{xr.TE_BUNCH}', 'N', f'{xr.N_PANEL}', '', '',
                'OPER', 'ITER 200', 'VPAR', f'N {n_crit}', '',
                f'VISC {Re:.6g}', 'PACC', 'polar.txt', '',
                f'ASEQ {a0} {a1} {da}', 'PACC', '', 'QUIT']
        subprocess.run([xr.XFOIL_EXE], input='\n'.join(cmds) + '\n', cwd=tmp,
                       text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, timeout=600)
        pol = xr._read_polar(os.path.join(tmp, 'polar.txt'))
        return xr.clmax_from_polar(pol, ALPHA_SEQ)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    print('Afericao do XFoil contra dados experimentais (Abbott & von')
    print('Doenhoff, TR-824, perfis lisos, Re = 8,9e6)\n')
    print(f'{"perfil":>8s} {"exp.":>6s} {"XFoil":>7s} {"desvio":>8s} '
          f'{"a_estol":>8s} {"situacao":>14s}')

    desvios = []
    for naca, Re, cl_exp in REFERENCIA:
        r = roda_naca(naca, Re)
        if not np.isfinite(r['clmax']):
            print(f'{naca:>8s} {cl_exp:6.2f} {"--":>7s}   nao convergiu')
            continue
        d = r['clmax'] - cl_exp
        desvios.append(d)
        print(f'{naca:>8s} {cl_exp:6.2f} {r["clmax"]:7.3f} {d:+8.3f} '
              f'{r["alpha_clmax"]:7.1f}o {r["situacao"]:>14s}')

    if not desvios:
        return 1
    dv = np.array(desvios)
    print(f'\ndesvio medio      : {dv.mean():+.3f}')
    print(f'desvio mediano    : {np.median(dv):+.3f}')
    print(f'faixa             : {dv.min():+.3f} a {dv.max():+.3f}')

    print(f'\n{"="*66}\nO QUE ISSO SIGNIFICA PARA A RESTRICAO\n{"="*66}')
    vies = float(np.median(dv))
    alvo_xfoil = 1.80 + vies
    print(f'O nosso XFoil le cl_max cerca de {vies:+.2f} acima do experimental.')
    print(f'Logo, um perfil com cl_max REAL de 1,80 (o alvo do designTool)')
    print(f'deve marcar aproximadamente {alvo_xfoil:.2f} na nossa medicao.')
    print(f'\nO perfil otimizado da MAC marcou 2,058.')
    if 2.058 >= alvo_xfoil:
        print(f'-> ACIMA do alvo corrigido ({alvo_xfoil:.2f}): ha folga, e o')
        print(f'   limiar pode ser afrouxado para recuperar arrasto.')
    else:
        print(f'-> ABAIXO do alvo corrigido ({alvo_xfoil:.2f}): o limiar NAO')
        print(f'   estava apertado demais; estava quase justo.')
    print('\nRessalva: os perfis de referencia sao NACA de 4 e 5 digitos, mais')
    print('finos e nao supercriticos. O vies pode diferir para a nossa')
    print('familia. Isso limita a precisao da correcao, mas o SINAL e o que')
    print('decide se faz sentido afrouxar o limiar.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
