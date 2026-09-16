'''
Helper comum dos scripts do Lab 04: roda o avl.exe em lote (stdin) e
extrai numeros da saida por expressao regular.
'''

# IMPORTS
import os
import re
import subprocess

AVL_EXE = os.path.abspath(os.path.join('avl', 'avl.exe'))


def roda_avl(comandos, timeout=580):
    '''Executa o AVL na pasta avl/ com a sequencia de comandos dada e
    devolve o texto completo da saida.'''
    r = subprocess.run([AVL_EXE], cwd='avl', input=comandos,
                       capture_output=True, text=True, timeout=timeout)
    return r.stdout


def pega(texto, padrao, indice=-1):
    '''Ultimo (ou i-esimo) numero que casa com o padrao, como float.'''
    achados = re.findall(padrao, texto)
    if not achados:
        return float('nan')
    return float(achados[indice])


def pega_todos(texto, padrao):
    '''Todos os numeros que casam com o padrao, como floats.'''
    return [float(v) for v in re.findall(padrao, texto)]
