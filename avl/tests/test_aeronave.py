import pytest

from aeronave import (aeronave, aeronave_mod, lc_h_maximo, W_PROJETO_BASE_KGF,
                      FOLGA_FUSELAGEM)
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


def test_aeronave_mod_sem_overrides_e_a_aeronave_entregue():
    base = aeronave_mod()
    ref = aeronave(base['Lc_h'])
    for k in ('W0', 'W', 'CL', 'CD0', 'xcg_fwd', 'xcg_aft', 'xnp_dt', 'deda'):
        assert base[k] == pytest.approx(ref[k], rel=1e-12)
    assert base['EH'] == pytest.approx(ref['EH'])


def test_aeronave_mod_com_lc_h_igual_a_aeronave():
    a = aeronave_mod({'Lc_h': 4.6})
    b = aeronave(4.6)
    assert a['W0'] == pytest.approx(b['W0'], rel=1e-12)
    assert a['EH']['xr'] == pytest.approx(b['EH']['xr'])


def test_aeronave_mod_aplica_overrides_sem_afetar_o_cache():
    ref = aeronave(4.6)
    mod = aeronave_mod({'Lc_h': 4.6, 'Cht': 0.8, 'xr_w': 16.5})
    assert mod['EH']['S'] == pytest.approx(ref['EH']['S']*0.8/0.7, rel=1e-9)
    assert mod['asa']['xr'] == pytest.approx(16.5)
    assert mod['W0'] > ref['W0']
    assert aeronave(4.6)['EH']['S'] == pytest.approx(ref['EH']['S'])


def test_lc_h_maximo_com_overrides():
    lc = lc_h_maximo(lo=3.0, overrides={'xr_w': 16.9966})
    av = aeronave_mod({'xr_w': 16.9966, 'Lc_h': lc})
    fim_eh = av['EH']['xr'] + av['EH']['cr']
    assert fim_eh == pytest.approx(av['fuselagem']['L'] - FOLGA_FUSELAGEM, abs=1e-3)
    assert lc < lc_h_maximo()
