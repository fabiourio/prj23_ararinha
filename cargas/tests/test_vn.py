import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import vn  # noqa: E402
from designTool.constants import ft2m, lb2N  # noqa: E402

FT_S = ft2m   # 1 ft/s em m/s


def test_u_ref_pontos_da_norma():
    assert vn.u_ref(0.0) == pytest.approx(56.0*FT_S)
    assert vn.u_ref(15000*ft2m) == pytest.approx(44.0*FT_S)
    assert vn.u_ref(60000*ft2m) == pytest.approx(20.86*FT_S)
    # linear entre 15000 e 60000 ft
    assert vn.u_ref(37500*ft2m) == pytest.approx(0.5*(44.0 + 20.86)*FT_S)


def test_n_max_manobra():
    assert vn.n_max_manobra(50000*lb2N) == pytest.approx(2.5)
    assert vn.n_max_manobra(640000*lb2N) == pytest.approx(2.5)
    assert vn.n_max_manobra(4000*lb2N) == pytest.approx(3.8)
    assert vn.n_max_manobra(20000*lb2N) == pytest.approx(2.1 + 24000/30000)


def test_fator_de_alivio_limites():
    assert vn.fator_alivio(1e9) == pytest.approx(0.88, rel=1e-6)
    assert vn.fator_alivio(5.3) == pytest.approx(0.44)


def test_dn_rajada_linear_em_V_e_U():
    a = vn.dn_rajada(100.0, 10.0, 5000.0, 6.0, 0.8)
    assert vn.dn_rajada(200.0, 10.0, 5000.0, 6.0, 0.8) == pytest.approx(2*a)
    assert vn.dn_rajada(100.0, 5.0, 5000.0, 6.0, 0.8) == pytest.approx(a/2)
    assert a == pytest.approx(0.8*vn.RHO0*10*100*6/(2*5000))


def test_envelope_coerente():
    W, S, c, h = 2.0e6, 370.0, 7.0, 10668.0
    e = vn.envelope(W, S, c, h, CLmax=1.2, CLmin=-1.0, CLa=6.5, M_C=0.85)
    assert e['V_A'] == pytest.approx(e['V_S1']*np.sqrt(e['n_max']))
    assert e['M']['V_C'] == pytest.approx(0.85)
    assert e['M']['V_D'] == pytest.approx(0.92)
    assert e['V_B'] <= e['V_C'] + 1e-9
    # V_B: curva de estol cruza a reta de rajada de V_C
    n_g = 1 + vn.dn_rajada(e['V_C'], e['U_ref'], W/S, 6.5, e['Kg'])
    assert e['V_B'] == pytest.approx(min(e['V_S1']*np.sqrt(n_g), e['V_C']))
    # rajada em V_D com metade da referencia
    assert e['U']['V_D'] == pytest.approx(0.5*e['U_ref'])
    # velocidade verdadeira: EAS/sqrt(sigma)
    sig = e['rho']/vn.RHO0
    assert e['TAS']['V_C'] == pytest.approx(e['V_C']/np.sqrt(sig))


def test_contorno_manobra_respeita_estol_e_n_max():
    e = vn.envelope(2.0e6, 370.0, 7.0, 0.0, CLmax=1.2, CLmin=-1.0, CLa=6.5, M_C=0.85)
    V, n = vn.contorno_manobra(e)
    pos = n > 0
    assert np.all(n[pos] <= e['n_max'] + 1e-9)
    assert np.all(n[pos] <= (V[pos]/e['V_S1'])**2 + 1e-9)
    assert np.all(n >= -1.0 - 1e-9)
