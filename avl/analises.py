'''Compensacao no cruzeiro: incidencia da EH que zera o profundor.'''

from avl_run import caso


def it_para_de_zero(arquivo, mach, cl, it0=0.0, passo=-2.0, tol=0.01, max_iter=8):
    '''
    Secante sobre it ate o profundor de compensacao (d2 pm 0) ficar abaixo de
    `tol` graus. Devolve (it, resultado do caso com derivadas).
    '''
    it_a = it0
    r_a = caso(arquivo, mach, cl=cl, it=it_a, trim=True, derivadas=True)
    if abs(r_a['de']) < tol:
        return it_a, r_a
    it_b = it0 + passo
    r_b = caso(arquivo, mach, cl=cl, it=it_b, trim=True, derivadas=True)
    for _ in range(max_iter):
        if abs(r_b['de']) < tol:
            return it_b, r_b
        it_c = it_b - r_b['de']*(it_b - it_a)/(r_b['de'] - r_a['de'])
        it_a, r_a = it_b, r_b
        it_b = it_c
        r_b = caso(arquivo, mach, cl=cl, it=it_b, trim=True, derivadas=True)
    raise RuntimeError(f'it nao convergiu: de = {r_b["de"]:.4f} com it = {it_b:.4f}')
