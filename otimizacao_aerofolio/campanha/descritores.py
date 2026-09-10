'''
Descritores geometricos de um perfil CST -- candidatos a preditor de cl_max.

A ideia e nao pre-comprometer a restricao do otimizador com uma unica
variavel. Registramos varios descritores, todos analiticos e suaves (logo
utilizaveis como g(x) num otimizador de gradiente), e deixamos o DOE decidir
qual deles -- ou qual combinacao -- realmente prediz cl_max.

Descritores implementados:

  r_LE            raio de bordo de ataque, r/c = A[0]^2/2. Sai da funcao de
                  classe do CST (N1=0.5): perto do bordo, y -> A[0]*sqrt(x),
                  e a parabola y^2 = A^2 x tem raio A^2/2 no vertice.

  delta_y         parametro de afiamento de bordo de ataque de Abbott & von
                  Doenhoff, usado no DATCOM/Roskam para correlacionar cl_max e
                  tipo de estol: diferenca de ordenada do extradorso entre
                  x/c = 6% e x/c = 0.15%, em porcentagem da corda.

  t_01, t_05      espessura relativa em x/c = 1% e 5% (bluntez do nariz).

  t_max, x_tmax   espessura maxima e posicao.
  c_max, x_cmax   arqueamento maximo e posicao.
  t_min           espessura minima (negativa = superficies cruzadas).

  camber_085      arqueamento em x/c = 0.85 e
  dcamber_te      inclinacao da linha media no bordo de fuga -- os dois medem
                  carregamento traseiro, caracteristico de perfis
                  supercriticos, que sobe o cl mas antecipa a separacao de
                  bordo de fuga e portanto pode DERRUBAR cl_max.
'''

import numpy as np


def _superficies(Au, Al, x):
    '''Ordenadas do extradorso e do intradorso nas estacoes x (x/c em [0,1]).'''
    from scipy.special import comb

    x = np.atleast_1d(np.asarray(x, dtype=float))
    C = x ** 0.5 * (1.0 - x) ** 1.0

    def _skin(A):
        A = np.asarray(A, dtype=float)
        N = len(A) - 1
        S = np.zeros_like(x)
        for ii in range(N + 1):
            S = S + A[ii] * comb(N, ii) * x ** ii * (1.0 - x) ** (N - ii)
        return C * S

    return _skin(Au), _skin(Al)


def le_radius(A0):
    '''Raio de bordo de ataque adimensional: r/c = A[0]^2 / 2.'''
    return float(A0) ** 2 / 2.0


def delta_y(Au):
    '''
    Parametro de afiamento de bordo de ataque (Abbott & von Doenhoff),
    em porcentagem da corda. Quanto maior, mais "cheio" o nariz e maior
    tende a ser cl_max.
    '''
    yu, _ = _superficies(Au, Au, np.array([0.0015, 0.06]))
    return float((yu[1] - yu[0]) * 100.0)


def descritores(Au, Al):
    '''
    Calcula todos os descritores de um perfil CST. Devolve um dicionario.

    ATENCAO: espessura maxima e sua posicao sao calculadas aqui, e nao lidas
    do cstfoil do airfoil_mod. O cstfoil filtra o vetor de espessuras para
    x > 0.05 (linha 216) mas depois indexa o vetor de abscissas COMPLETO com
    o indice do vetor filtrado (linha 219), entao x_max_thickness sai
    deslocado -- para o NACA 1411 ele devolve 0.117 quando o valor correto e
    0.30. O valor de max_thickness em si esta certo; so a posicao esta errada.
    O arqueamento nao sofre do problema, porque o vetor cc nao e filtrado.
    Isso afeta a coluna x_t/c,max da Tab. 2 do roteiro.
    '''
    Au = np.asarray(Au, dtype=float)
    Al = np.asarray(Al, dtype=float)

    # geometria com adensamento nos bordos
    x = (1 - np.cos(np.linspace(0, 1, 801) * np.pi)) / 2
    yu, yl = _superficies(Au, Al, x)
    esp = yu - yl
    camber = 0.5 * (yu + yl)

    # espessura maxima e posicao, na mesma janela que o cstfoil usa
    janela = (x > 0.05) & (x < 0.97)
    xj, ej = x[janela], esp[janela]
    it = int(np.argmax(ej))

    ic = int(np.argmax(np.abs(camber)))

    def _em(xq):
        '''
        Avalia a espessura DIRETAMENTE na estacao pedida, sem interpolar.
        Precisa bater bit a bit com o t01 analitico do otimizador, senao o
        limiar calibrado aqui nao e o mesmo que a restricao aplica la.
        '''
        yu_q, yl_q = _superficies(Au, Al, np.array([xq]))
        return float(yu_q[0] - yl_q[0])

    # inclinacao da linha media perto do bordo de fuga (carregamento traseiro)
    m = x > 0.90
    dcamber_te = float(np.polyfit(x[m], camber[m], 1)[0])

    return {
        'r_LE_sup': le_radius(Au[0]),
        'r_LE_inf': le_radius(Al[0]),
        'delta_y': delta_y(Au),
        't_01': _em(0.01),
        't_05': _em(0.05),
        't_max': float(ej[it]),
        'x_tmax': float(xj[it]),
        't_min': float(np.min(ej)),
        'c_max': float(camber[ic]),
        'x_cmax': float(x[ic]),
        'camber_085': float(np.interp(0.85, x, camber)),
        'dcamber_te': dcamber_te,
    }


NOMES = ['r_LE_sup', 'r_LE_inf', 'delta_y', 't_01', 't_05',
         't_max', 'x_tmax', 't_min', 'c_max', 'x_cmax',
         'camber_085', 'dcamber_te']
