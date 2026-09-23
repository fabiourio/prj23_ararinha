import numpy as np
import pytest

from aeronave import aeronave
from avl_run import caso
from gera_avl import escreve_avl, secoes_asa, DZ, MALHA_PADRAO

REF = 'tests/dados/fwd_referencia_11c6802.avl'


def test_secoes_batem_com_o_arquivo_de_referencia():
    s = secoes_asa(aeronave(4.6))
    assert len(s) == 11
    assert s[0]['z'] == pytest.approx(-2.5, abs=1e-4)          # -1,3 + DZ asa
    assert s[-1]['x'] == pytest.approx(38.6635, abs=1e-3)
    assert s[-1]['z'] == pytest.approx(0.6611, abs=1e-3)
    assert s[1]['y'] == pytest.approx(3.0407, abs=1e-3)
    assert [x['aileron'] for x in s] == [False]*6 + [True]*4 + [False]


def test_regenera_o_modelo_de_referencia():
    escreve_avl(aeronave(4.6), 'resultados/_tmp/teste_fwd.avl', cg='fwd',
                malha=MALHA_PADRAO)
    for a in (0.0, 4.0):
        novo = caso('resultados/_tmp/teste_fwd.avl', 0.85, alfa=a)
        ref = caso(REF, 0.85, alfa=a)
        # O modelo congelado (commit 11c6802) tem 3 secoes da asa arredondadas
        # de forma diferente na 4a casa decimal (0,1 mm: eta 0,22, 0,56, 0,70,
        # 0,82). Isso sozinho explica |dCL| = 7e-5 e |dCm| = 4e-5 em alfa =
        # 4 graus (corrigindo so essas linhas no .avl gerado, o resultado bate
        # exatamente com a referencia). 1e-4 fica bem abaixo de qualquer
        # diferenca fisicamente relevante. CDff fica com verificacao mais
        # apertada (2e-6) porque nao e afetado por esse arredondamento.
        assert novo['CDff'] == pytest.approx(ref['CDff'], abs=2e-6)
        for k in ('CL', 'Cm', 'CD'):
            assert novo[k] == pytest.approx(ref[k], abs=1e-4), k


def test_modelo_so_asa_roda():
    escreve_avl(aeronave(4.6), 'resultados/_tmp/teste_asa.avl', so_asa=True,
                diedro=False)
    r = caso('resultados/_tmp/teste_asa.avl', 0.85, cl=0.5, empenagem=False)
    assert r['CL'] == pytest.approx(0.5, abs=1e-4)
    assert 0.5 < r['e'] < 1.05


def test_torcao_entra_no_ainc():
    import os
    from gera_avl import AQUI
    tw = np.linspace(0, -4, 11)
    escreve_avl(aeronave(4.6), 'resultados/_tmp/teste_tw.avl', torcao=tw)
    with open(os.path.join(AQUI, 'resultados', '_tmp', 'teste_tw.avl')) as f:
        assert '-4.0000' in f.read()
