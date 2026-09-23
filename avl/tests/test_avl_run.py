import os

import numpy as np
import pytest

from avl_run import le_resultados, le_faixas, caso

AMOSTRA = os.path.join(os.path.dirname(__file__), 'dados', 'amostra_avl.txt')


def texto():
    with open(AMOSTRA) as f:
        return f.read()


def test_le_forcas_totais():
    r = le_resultados(texto())
    assert r['alfa'] == pytest.approx(4.22269)
    assert r['CL'] == pytest.approx(0.5)
    assert r['CDff'] == pytest.approx(0.0113659)
    assert r['CDvis'] == pytest.approx(0.01473)
    assert r['e'] == pytest.approx(0.6928)
    assert r['Cm'] == pytest.approx(0.0, abs=1e-5)
    assert r['de'] == pytest.approx(-10.58309)   # nao confunde com a lista de restricoes
    assert r['it'] == pytest.approx(0.0)
    assert r['xnp'] == pytest.approx(27.427111)
    assert r['CLa'] == pytest.approx(6.967041)


def test_le_faixas_da_asa_direita():
    f = le_faixas(texto(), superficie=1)
    assert len(f['y']) == 20
    assert f['y'][0] == pytest.approx(0.0491)
    assert f['cl_norm'][0] == pytest.approx(0.4497)
    assert f['ccl'][-1] == pytest.approx(2.2154)
    assert np.all(np.diff(f['y']) > 0)


def test_caso_roda_de_verdade():
    r = caso('fwd.avl', mach=0.85, alfa=3.0, it=0.0, de=0.0)
    assert np.isfinite(r['CL']) and r['alfa'] == pytest.approx(3.0)
