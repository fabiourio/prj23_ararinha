import numpy as np
import pytest

from aeronave import aeronave
from otimizacao_torcao import Avaliador, elipse_normalizada


def test_elipse_normalizada_integra_o_CL():
    eta = np.linspace(0, 1, 20001)
    assert np.trapezoid(elipse_normalizada(eta), eta) == pytest.approx(np.pi/4, rel=1e-4)


def test_avaliador_asa_limpa_tem_cache_e_bate_o_CL():
    av = aeronave(4.6)
    a = Avaliador(av, modelo='asa_limpa', com_estol=False, nome='teste_ot')
    x = np.zeros(10)
    f1 = a.objetivo(x)
    n = a.n_rodadas
    f2 = a.objetivo(x.copy())
    assert f1 == f2 and a.n_rodadas == n          # segunda chamada vem do cache
    assert a.ultimo['r']['CL'] == pytest.approx(av['CL'], abs=1e-4)
