import os

import numpy as np
import pytest

from avl_saida import derivadas, casos_ft, por_rad, converte, sessao
from derivadas_estabilidade import dif_finita_q, tabela9

AMOSTRA_SB = os.path.join(os.path.dirname(__file__), 'dados', 'amostra_sb.txt')


def texto():
    with open(AMOSTRA_SB) as f:
        return f.read()


def test_le_st():
    st = derivadas(texto(), 'st')
    assert st['CLa'] == pytest.approx(6.950467)
    assert st['Cma'] == pytest.approx(0.213541)
    assert st['CLq'] == pytest.approx(6.822587)
    assert st['Cmq'] == pytest.approx(-36.561115)
    assert st['CYb'] == pytest.approx(-0.767324)
    assert st['Clp'] == pytest.approx(-0.595339)      # eixos de estabilidade
    assert st['CLd2'] == pytest.approx(0.006095)
    assert st['CDffd2'] == pytest.approx(0.000003)
    assert st['CLg1'] == pytest.approx(0.011207)
    assert st['Cmg1'] == pytest.approx(-0.052034)
    assert st['CDffg1'] == pytest.approx(-0.000012)
    assert st['Xnp'] == pytest.approx(27.508125)


def test_st_nao_confunde_cnb_com_razao_espiral():
    # linha "Clb Cnr / Clr Cnb  =  -2.947581" logo depois do Xnp
    assert derivadas(texto(), 'st')['Cnb'] == pytest.approx(-0.081488)


def test_le_sb():
    sb = derivadas(texto(), 'sb')
    assert sb['Clp'] == pytest.approx(-0.605934)      # eixos do corpo (difere do st)
    assert sb['Clr'] == pytest.approx(0.193310)
    assert sb['Cnp'] == pytest.approx(-0.059073)
    assert sb['Cnr'] == pytest.approx(-0.208283)
    assert sb['Cld1'] == pytest.approx(-0.004183)
    assert sb['Cld3'] == pytest.approx(-0.000530)
    assert sb['Cnd1'] == pytest.approx(-0.000032)
    assert sb['Cnd3'] == pytest.approx(0.002192)
    assert sb['CXw'] == pytest.approx(0.735076)       # X para a frente: CXw > 0
    assert sb['CZg1'] == pytest.approx(-0.011215)
    assert 'Xnp' not in sb


def test_sb_ausente_gera_erro():
    s = texto()
    with pytest.raises(RuntimeError):
        derivadas(s[:s.find('Geometry-axis')], 'sb')


def test_casos_ft_um_por_execucao():
    rs = casos_ft(texto())
    assert len(rs) == 3
    assert all(r['alfa'] == pytest.approx(3.79107) for r in rs)
    assert rs[0]['it'] == pytest.approx(-2.0)
    assert rs[0]['CDind'] == pytest.approx(0.0059038)


def test_conversoes():
    assert por_rad(1.0) == pytest.approx(180/np.pi)
    assert converte(0.006095, por_grau=True) == pytest.approx(0.349218, rel=1e-5)
    assert converte(-0.767324, inverte=True) == pytest.approx(0.767324)
    assert converte(-0.004183, inverte=True, por_grau=True) == pytest.approx(0.239668, rel=1e-5)


def test_tabela9_convencoes():
    s = texto()
    st, sb = derivadas(s, 'st'), derivadas(s, 'sb')
    t4 = [{'parametro': 'zp', 'explicacao': '', 'valor': 3.0, 'fonte': 'designTool'}]
    t6 = [{'parametro': 'CL0', 'valor': 0.1}, {'parametro': 'CM0', 'valor': -0.2}]
    t7 = [{'parametro': k, 'valor': v} for k, v in (('CD0', 0.02), ('CDa', 0.01), ('CDa2', 0.5))]
    t9 = {l['parametro']: l for l in tabela9(t4, t6, t7, st, sb, {'CDq': 0.3})}
    g = 180/np.pi
    assert t9['CLde']['valor_final'] == pytest.approx(st['CLd2']*g)
    assert t9['CDit']['valor_final'] == pytest.approx(st['CDffg1']*g)
    assert t9['CYb']['valor_final'] == pytest.approx(-st['CYb'])
    assert t9['CYdr']['valor_final'] == pytest.approx(st['CYd3']*g)
    assert t9['Clb']['valor_final'] == pytest.approx(st['Clb'])
    assert t9['Clp']['valor_final'] == pytest.approx(sb['Clp'])
    assert t9['Clda']['valor_final'] == pytest.approx(-sb['Cld1']*g)
    assert t9['Cndr']['valor_final'] == pytest.approx(-sb['Cnd3']*g)
    assert t9['CDq']['valor_final'] == pytest.approx(0.3)
    assert t9['zp']['conversao'] == 'sinal invertido'


def test_diferenca_finita_em_q_reproduz_clq_do_st():
    av = {'M': 0.85}
    alfa, it = 3.0, -4.0
    s = sessao('aft.avl', av['M'], it, ['d2 d2 0', f'a a {alfa}', 'x', 'st', ''])
    st = derivadas(s, 'st')
    fd, _ = dif_finita_q('aft.avl', av, it, alfa)
    assert fd['CLq'] == pytest.approx(st['CLq'], rel=0.02)
    assert fd['Cmq'] == pytest.approx(st['Cmq'], rel=0.02)
    assert np.isfinite(fd['CDq'])
