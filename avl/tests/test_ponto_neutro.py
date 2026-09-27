import pytest

from aeronave import aeronave
from ponto_neutro import np_designtool


def test_np_designtool_reproduz_o_ponto_neutro_do_designtool():
    # np_designtool replica designTool/balance.py, mas com o Mach como
    # parametro livre. No Mach_cruise (M da aeronave), tem de bater
    # exatamente com o xnp que o proprio designTool calculou (bal['xnp']).
    av = aeronave(4.6)
    d = np_designtool(4.6, av['M'])
    assert d['xnp'] == pytest.approx(av['xnp_dt'], abs=1e-6)
