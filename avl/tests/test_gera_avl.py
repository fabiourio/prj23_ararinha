import os

import numpy as np
import pytest

from aeronave import aeronave
from avl_run import caso
from gera_avl import (escreve_avl, secoes_asa, reamostra_fuselagem, AQUI, DZ,
                      MALHA_PADRAO)

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
    # Modelo congelado (commit 11c6802): nacele no mesmo COMPONENT da asa e
    # fuselagem original (sem a reamostragem), como era escrito na epoca.
    escreve_avl(aeronave(4.6), 'resultados/_tmp/teste_fwd.avl', cg='fwd',
                malha=MALHA_PADRAO, nacele_componente=1, fuselagem='original')
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
    tw = np.linspace(0, -4, 11)
    escreve_avl(aeronave(4.6), 'resultados/_tmp/teste_tw.avl', torcao=tw)
    with open(os.path.join(AQUI, 'resultados', '_tmp', 'teste_tw.avl')) as f:
        assert '-4.0000' in f.read()


def test_fuselagem_reamostrada_e_fechada_e_fiel_a_original():
    destino = reamostra_fuselagem(destino='resultados/_tmp/teste_fus_reamostrada.dat')
    caminho = os.path.join(AQUI, destino)
    with open(caminho) as f:
        linhas = f.readlines()
    pts = np.array([[float(v) for v in l.split()] for l in linhas[1:]])

    n = 81
    assert len(pts) == 2*n - 1                       # fecha TE->cima->LE->baixo->TE
    assert pts[0, 0] == pytest.approx(1.0)            # TE
    assert pts[n - 1, 0] == pytest.approx(0.0)        # LE
    assert pts[-1, 0] == pytest.approx(1.0)           # TE de novo (fechado)

    cima = pts[:n]
    baixo = pts[n - 1:]
    assert np.all(np.diff(cima[:, 0]) < 0)            # x decresce TE->LE
    assert np.all(np.diff(baixo[:, 0]) > 0)            # x cresce LE->TE

    orig = np.loadtxt(os.path.join(AQUI, 'fuselage_nondim.dat'), skiprows=1)
    i0 = int(np.argmin(orig[:, 0]))
    orig_cima = orig[:i0 + 1][::-1]
    orig_baixo = orig[i0:]
    z_cima_lin = np.interp(cima[:, 0], orig_cima[:, 0], orig_cima[:, 1])
    z_baixo_lin = np.interp(baixo[:, 0], orig_baixo[:, 0], orig_baixo[:, 1])
    # Tolerancia = arredondamento do formato texto (%.6f) propagado pela
    # interpolacao linear, nao erro de metodo.
    assert np.max(np.abs(cima[:, 1] - z_cima_lin)) == pytest.approx(0.0, abs=1e-4)
    assert np.max(np.abs(baixo[:, 1] - z_baixo_lin)) == pytest.approx(0.0, abs=1e-4)


def test_modelo_padrao_usa_nacele_em_componente_2_e_fuselagem_reamostrada():
    escreve_avl(aeronave(4.6), 'resultados/_tmp/teste_padrao.avl', malha=MALHA_PADRAO)
    with open(os.path.join(AQUI, 'resultados', '_tmp', 'teste_padrao.avl')) as f:
        txt = f.read()
    assert 'COMPONENT\n2' in txt
    assert 'fuselage_reamostrada.dat' in txt
