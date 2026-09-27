'''
D: faixas da ponta/juncao e da raiz. Para uma torcao (arg: armazenado ou
caminho de um json com 'twist'), CG traseiro, it retrimado:
  - alfa_i de todas as faixas (inclusive excluidas) em 4 malhas asa_ns;
  - estol (alfa, eta) com exclusao ct/8, ct/4, ct/2 em cada malha;
  - cl_norm/clmax na ultima faixa da asa e largura da faixa, vs malha.
Saida: runs/ponta_<rotulo>.json
'''

import json
import os
import sys

import numpy as np

from comum import aeronave, LC, RUNS, AVL_DIR, gera_avl, escreve_avl, it_para_de_zero, caso
from estol import modelo_estol, clmax_local, MACH_BAIXO

rotulo = sys.argv[1]
fonte = sys.argv[2]
if fonte == 'armazenado':
    with open(os.path.join(AVL_DIR, 'resultados', 'torcao_otimizada.json')) as f:
        tw = np.concatenate([[0.0], json.load(f)['twist']])
elif fonte == 'zero':
    tw = np.zeros(11)
else:
    with open(fonte) as f:
        tw = np.array(json.load(f)['twist'])
        if len(tw) == 10:
            tw = np.concatenate([[0.0], tw])

av = aeronave(LC)
semi, ct = av['asa']['yt'], av['asa']['ct']
out = {'twist': tw.tolist(), 'malhas': {}}
for ns in (30, 60, 90, 120):
    malha = gera_avl.malha_adotada(); malha['asa_ns'] = ns
    arq = escreve_avl(av, f'resultados/_tmp/inv_ponta_{rotulo}_{ns}.avl', cg='aft',
                      torcao=tw, malha=malha)
    it, _ = it_para_de_zero(arq, av['M'], av['CL'])
    m = modelo_estol(arq, it, trim=True, av=av)
    eta, al = m['eta'], m['alfa_i']
    y = eta*semi
    ent = {'it': it, 'eta': eta.tolist(), 'alfa_i': al.tolist(), 'excl': {}}
    for nome, fr in (('ct/8', 1/8), ('ct/4', 1/4), ('ct/2', 1/2), ('nenhuma', 0.0)):
        d = fr*ct
        exc = ((semi - y) <= d) | (y <= d) if fr > 0 else np.zeros_like(eta, bool)
        i = int(np.where(~exc)[0][np.argmin(al[~exc])])
        ent['excl'][nome] = {'alfa': float(al[i]), 'eta': float(eta[i]),
                             'n_excluidas': int(exc.sum())}
    a_st = ent['excl']['ct/4']['alfa']
    r = caso(arq, MACH_BAIXO, alfa=a_st, it=it, trim=True, faixas=True)
    fs = r['faixas']
    ent['a_estol'] = a_st
    ent['cl_norm'] = fs['cl_norm'].tolist()
    ent['clmax'] = clmax_local(fs['y']/semi).tolist()
    ent['largura_ultima_m'] = float(semi - fs['y'][-1])*2
    ent['cl_ultima'] = float(fs['cl_norm'][-1])
    ent['cl_primeira'] = float(fs['cl_norm'][0])
    out['malhas'][ns] = ent
    print(f'ns {ns:3d}: it {it:.3f}  ' + '  '.join(
        f'{k}: α {v["alfa"]:.2f} η {v["eta"]:.3f} (n_exc {v["n_excluidas"]})'
        for k, v in ent['excl'].items())
        + f'  | a α_estol: cl_norm ultima faixa {ent["cl_ultima"]:.3f}, primeira {ent["cl_primeira"]:.3f}, '
          f'η ultima {fs["y"][-1]/semi:.4f}', flush=True)

with open(os.path.join(RUNS, f'ponta_{rotulo}.json'), 'w') as f:
    json.dump(out, f, indent=1, default=float)
print('ok')
