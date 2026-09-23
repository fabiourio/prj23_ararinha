import numpy as np
import pytest

from aeronave import aeronave
from analises import it_para_de_zero
from estol import clmax_local, estol, ETA_CLMAX, CLMAX
from gera_avl import escreve_avl


def test_clmax_local_interpola_e_satura():
    assert clmax_local(0.0) == pytest.approx(CLMAX[0])
    assert clmax_local(1.0) == pytest.approx(CLMAX[-1])
    meio = 0.5*(ETA_CLMAX[0] + ETA_CLMAX[1])
    assert clmax_local(meio) == pytest.approx(0.5*(CLMAX[0] + CLMAX[1]))


def test_it_zera_o_profundor():
    av = aeronave(4.6)
    arq = escreve_avl(av, 'resultados/_tmp/teste_an.avl', cg='aft')
    it, r = it_para_de_zero(arq, av['M'], av['CL'])
    assert abs(r['de']) < 0.01
    assert r['CL'] == pytest.approx(av['CL'], abs=1e-4)
    assert np.isfinite(r['xnp'])


def test_estol_modelo_linear_bate_com_a_conferencia_direta():
    av = aeronave(4.6)
    arq = escreve_avl(av, 'resultados/_tmp/teste_an.avl', cg='fwd')
    e = estol(arq, it=0.0, trim=True, semi=av['asa']['yt'])
    assert abs(e['excesso_max']) < 0.02     # a faixa critica chega ao clmax
    assert 5.0 < e['alfa'] < 25.0
