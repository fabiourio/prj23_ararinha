'''
Lab 04, item 3 -- incidencia da empenagem que anula o profundor em cruzeiro.

Para cada arquivo (fwd.avl e aft.avl), fixa a meta de sustentacao do ponto
de projeto (M = 0,85, CL = 0,5053) e a trimagem de arfagem pelo profundor
(d2 pm 0), e resolve a incidencia it que zera a deflexao de profundor.
Como a resposta delta_e(it) e linear no VLM, bastam duas rodadas para achar
o zero por secante, mais uma de verificacao.

Gera:
  - avl/saidas/trim.json (usado pelos demais scripts do Lab 04)
  - avl/saidas/trim_fwd.txt e trim_aft.txt (bloco de forcas da rodada
    trimada, evidencia para o relatorio)

Rodar da raiz do repo:  python lab04_trim.py
'''

# IMPORTS
import json
import os

from avl_batch import roda_avl, pega

#=========================================

CL_PROJ = 0.5053
MACH = 0.85


def rodada(cg, it):
    '''Roda um caso trimado em arfagem com a incidencia dada e devolve a
    deflexao de profundor resultante e a saida completa.'''
    cmds = (f'load {cg}.avl\n'
            'oper\n'
            f'm\nmn {MACH}\n\n'
            f'a c {CL_PROJ}\n'
            'd2 pm 0\n'
            f'de\n1 {it}\n\n'
            'x\n'
            '\nquit\n')
    saida = roda_avl(cmds)
    return pega(saida, r'elevator\s+=\s+([-\d.]+)'), saida


os.makedirs('avl/saidas', exist_ok=True)
resultado = {}
for cg in ('fwd', 'aft'):
    de0, _ = rodada(cg, 0.0)
    de1, _ = rodada(cg, -3.0)
    it_zero = 0.0 - de0*(-3.0 - 0.0)/(de1 - de0)
    de_res, saida = rodada(cg, it_zero)

    alpha = pega(saida, r'Alpha =\s+([-\d.]+)')
    cm = pega(saida, r'Cmtot =\s+([-\d.]+)')
    cl = pega(saida, r'CLtot =\s+([-\d.]+)')
    cd = pega(saida, r'CDtot =\s+([-\d.]+)')

    # Bloco de forcas totais da rodada final, como evidencia
    ini = saida.rfind('Vortex Lattice Output')
    fim = saida.find('---', ini + 100)
    with open(f'avl/saidas/trim_{cg}.txt', 'w', encoding='ascii') as f:
        f.write(saida[ini:fim])

    resultado[cg] = {'it_deg': round(it_zero, 4),
                     'delta_e_residual_deg': round(de_res, 5),
                     'alpha_deg': round(alpha, 4),
                     'CL': round(cl, 4), 'CD': round(cd, 5),
                     'Cm': round(cm, 5)}
    print(f"{cg}: it = {it_zero:.4f} graus, delta_e residual = "
          f"{de_res:.5f} graus, alpha = {alpha:.4f}, CL = {cl:.4f}, "
          f"Cm = {cm:.5f}")

with open('avl/saidas/trim.json', 'w', encoding='ascii') as f:
    json.dump(resultado, f, indent=2)
print('gravado: avl/saidas/trim.json')
