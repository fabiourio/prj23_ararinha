import pytest

from aeronave import aeronave, lc_h_maximo, W_PROJETO_BASE_KGF, FOLGA_FUSELAGEM
from designTool.constants import gravity


def test_base_reproduz_ponto_de_projeto_do_lab03():
    av = aeronave(4.6)
    assert av['W'] == pytest.approx(W_PROJETO_BASE_KGF*gravity, abs=1.0)
    assert av['CL'] == pytest.approx(0.5053, abs=1e-3)
    assert av['xcg_fwd'] == pytest.approx(26.0797, abs=1e-3)
    assert 0.0 < av['deda'] < 1.0


def test_alturas_sao_as_reais_do_designtool():
    av = aeronave(4.6)
    assert av['asa']['zr'] == pytest.approx(-1.3)
    assert av['EH']['zr'] == pytest.approx(2.0)
    assert av['EV']['zr'] == pytest.approx(3.0)


def test_eh_mais_para_tras_reduz_W0_e_area():
    a, b = aeronave(4.6), aeronave(4.8)
    assert b['W0'] < a['W0']
    assert b['EH']['S'] < a['EH']['S']


def test_limite_da_fuselagem():
    lc = lc_h_maximo()
    av = aeronave(lc)
    fim_eh = av['EH']['xr'] + av['EH']['cr']
    assert fim_eh == pytest.approx(av['fuselagem']['L'] - FOLGA_FUSELAGEM, abs=1e-3)


def test_retorno_nao_tem_aliasing_com_o_cache():
    av1 = aeronave(4.6)
    original = av1['EH']['xr']
    av1['EH']['xr'] = -999.0
    av2 = aeronave(4.6)
    assert av2['EH']['xr'] == pytest.approx(original)
