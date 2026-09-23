import os

import numpy as np
import pytest

from avl_run import le_resultados, le_faixas, caso

AMOSTRA = os.path.join(os.path.dirname(__file__), 'dados', 'amostra_avl.txt')


def texto():
    with open(AMOSTRA) as f:
        return f.read()


def texto_faixas_overflow():
    '''Bloco sintetico de uma superficie com uma faixa marcada com "******".'''
    return '''
  Surface # 1     Wing
     # Chordwise =  8   # Spanwise =  3     First strip =  1

 Strip Forces referred to Strip Area, Chord
    j      Yle    Chord     Area     c cl      ai      cl_norm  cl       cd       cdv    cm_c/4    cm_LE  C.P.x/c
     1   0.0491   2.2637   0.2224   1.0179   0.0740   0.4497   0.4494   0.0091   0.0000  -0.0982  -0.3527    0.218
     2   0.2145   2.2601   0.4425   ******   0.0813   0.4581   0.4577   0.0087   0.0000  -0.1013  -0.3600    0.221
     3   0.4127   2.2493   0.4395   1.0349   0.0895   0.4680   0.4675   0.0084   0.0000  -0.1054  -0.3691    0.226

  Surface # 2     Wing (YDUP)
     # Chordwise =  8   # Spanwise =  3     First strip =  4
'''



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


def test_le_faixas_conta_bate_com_spanwise_mesmo_com_coluna_em_branco():
    # nacele (YDUP), superficie 9: strips 88 e 91 tem C.P.x/c em branco
    # (cl_norm ~ 0); # Spanwise = 12 mas antes so 10 linhas eram lidas.
    f = le_faixas(texto(), superficie=9)
    assert len(f['y']) == 12
    assert np.isnan(f['cpx'][5])   # strip 88: 6a linha do bloco (j=88)


def test_le_faixas_overflow_gera_erro():
    with pytest.raises(RuntimeError):
        le_faixas(texto_faixas_overflow(), superficie=1)


def test_caso_trim_impossivel_levanta_erro():
    with pytest.raises(RuntimeError):
        caso('fwd.avl', mach=0.85, cl=5.0, trim=True)
