import numpy as np
import pytest

from aeronave import aeronave
from analises import it_para_de_zero
from estol import clmax_local, estol, ETA_CLMAX, CLMAX
from gera_avl import escreve_avl, MALHA_PADRAO


def test_clmax_local_interpola_e_satura():
    assert clmax_local(0.0) == pytest.approx(CLMAX[0])
    assert clmax_local(1.0) == pytest.approx(CLMAX[-1])
    meio = 0.5*(ETA_CLMAX[0] + ETA_CLMAX[1])
    assert clmax_local(meio) == pytest.approx(0.5*(CLMAX[0] + CLMAX[1]))


def test_it_zera_o_profundor():
    av = aeronave(4.6)
    arq = escreve_avl(av, 'resultados/_tmp/teste_an.avl', cg='aft', malha=MALHA_PADRAO)
    it, r = it_para_de_zero(arq, av['M'], av['CL'])
    assert abs(r['de']) < 0.01
    assert r['CL'] == pytest.approx(av['CL'], abs=1e-4)
    assert np.isfinite(r['xnp'])


def test_estol_modelo_linear_bate_com_a_conferencia_direta():
    # malha=MALHA_PADRAO (nao a adotada pela convergencia): tolerancia
    # calibrada contra essa malha, para o teste nao depender do resultado
    # do estudo de convergencia_malha.py.
    av = aeronave(4.6)
    arq = escreve_avl(av, 'resultados/_tmp/teste_an.avl', cg='fwd', malha=MALHA_PADRAO)
    e = estol(arq, it=0.0, trim=True, av=av)
    assert abs(e['excesso_max']) < 0.02     # a faixa critica chega ao clmax
    assert 5.0 < e['alfa'] < 25.0


def test_estol_ignora_a_faixa_da_juncao_com_o_winglet():
    # Sem a exclusao da(s) faixa(s) de canto, o alfa de estol cai de ~10 para
    # ~5 graus so por refinar a malha da envergadura (singularidade de
    # vórtice de canto na ultima faixa, ver docstring de estol.py).
    av = aeronave(4.6)
    arq20 = escreve_avl(av, 'resultados/_tmp/teste_estol_ns20.avl', cg='fwd',
                        malha=dict(MALHA_PADRAO, asa_ns=20))
    arq60 = escreve_avl(av, 'resultados/_tmp/teste_estol_ns60.avl', cg='fwd',
                        malha=dict(MALHA_PADRAO, asa_ns=60))
    e20 = estol(arq20, it=0.0, trim=True, av=av)
    e60 = estol(arq60, it=0.0, trim=True, av=av)
    assert abs(e60['alfa'] - e20['alfa']) < 0.3
