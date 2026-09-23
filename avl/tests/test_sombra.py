import numpy as np
import pytest

from sombra import alfa_saida, fracao_leme_encoberto


def av_teste(x_eh=40.0, tc=0.0, deda=0.5):
    return {'deda': deda,
            'asa': {'xr': 0.0, 'zr': 0.0, 'cr': 10.0, 'ct': 10.0, 'xt': 0.0,
                    'yt': 20.0, 'zt': 0.0},
            'EH': {'xr': x_eh, 'zr': 3.0, 'cr': 4.0, 'ct': 4.0, 'xt': x_eh,
                   'yt': 5.0, 'zt': 3.0, 'tc': tc},
            'EV': {'xr': 0.0, 'zr': 3.0, 'cr': 10.0, 'ct': 10.0, 'xt': 0.0,
                   'zt': 13.0}}


def test_alfa_de_saida_caso_simples():
    # EH 3 m acima e 30 m atras do bordo de fuga: gama = atan(3/30)
    esperado = np.degrees(np.arctan(3/30))/(1 - 0.5)
    assert alfa_saida(av_teste()) == pytest.approx(esperado, rel=1e-6)


def test_espessura_da_eh_conta():
    esperado = np.degrees(np.arctan(3.2/30))/0.5      # extradorso: 3 + 0,1*4/2
    assert alfa_saida(av_teste(tc=0.1)) == pytest.approx(esperado, rel=1e-6)


def test_eh_mais_para_tras_sai_da_sombra_antes():
    assert alfa_saida(av_teste(x_eh=50.0)) < alfa_saida(av_teste(x_eh=40.0))


def test_leme_livre_quando_eh_longe_para_frente():
    av = av_teste()
    av['EH'].update({'xr': -100.0, 'cr': 1.0, 'zr': 3.0})
    assert fracao_leme_encoberto(av) == pytest.approx(0.0)


def test_leme_todo_encoberto():
    # EV de 0 a 10 em x, charneira em 7; linhas de 60 graus (da raiz em x=0)
    # e de 30 graus (do bordo de fuga em x=10) cobrem o leme inteiro
    av = av_teste()
    av['EH'].update({'xr': 0.0, 'cr': 10.0, 'zr': 3.0})
    assert fracao_leme_encoberto(av) == pytest.approx(1.0)
