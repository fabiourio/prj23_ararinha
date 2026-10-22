import numpy as np
import pytest

import esforcos as es


def test_carga_uniforme():
    L, q = 10.0, 3.0
    y = np.linspace(0, L, 1001)
    V, M = es.integra(y, np.full(1000, q), [])
    assert V[0] == pytest.approx(q*L)
    assert M[0] == pytest.approx(q*L**2/2)
    assert V[-1] == pytest.approx(0.0)
    assert M[-1] == pytest.approx(0.0)
    # M(y) = q (L - y)^2 / 2
    assert M[500] == pytest.approx(q*(L - y[500])**2/2)


def test_carga_pontual():
    y = np.linspace(0, 10.0, 1001)
    V, M = es.integra(y, np.zeros(1000), [(4.0, -2.0)])
    assert V[0] == pytest.approx(-2.0)
    assert M[0] == pytest.approx(-8.0)
    i = np.searchsorted(y, 4.0) + 1        # logo para fora da carga
    assert V[i] == pytest.approx(0.0)
    assert M[i] == pytest.approx(0.0)


def test_torcao():
    y = np.linspace(0, 5.0, 501)
    T = es.integra_torcao(y, np.full(500, 2.0), [(3.0, 7.0)])
    assert T[0] == pytest.approx(2.0*5.0 + 7.0)
    assert T[-1] == pytest.approx(0.0)


def test_bandas_das_faixas():
    # tres faixas de largura 1, 2 e 3 (ao longo da superficie, sem diedro)
    corda = np.array([2.0, 2.0, 1.0])
    area = np.array([2.0, 4.0, 3.0])
    bordas = es.bordas_faixas(area, corda, diedro=0.0)
    assert bordas == pytest.approx([0.0, 1.0, 3.0, 6.0])
    # com diedro, a largura projetada encolhe por cos(diedro)
    b2 = es.bordas_faixas(area, corda, diedro=np.radians(6))
    assert b2[-1] == pytest.approx(6.0*np.cos(np.radians(6)))
