'''
Execucao do AVL por stdin e leitura da saida.

O AVL roda sempre com cwd = pasta avl/, e os caminhos passados sao
relativos a ela: o executavel nao lida bem com o "º" do caminho absoluto.

Numeracao do modelo: d1 aileron, d2 profundor, d3 leme; variavel de
projeto 1 = it (incidencia da EH).
'''

import os
import re
import subprocess

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
AVL = os.path.join(AQUI, 'avl.exe')

COLUNAS_FS = ['y', 'corda', 'area', 'ccl', 'ai', 'cl_norm', 'cl', 'cd',
              'cdv', 'cm_c4', 'cm_le', 'cpx']


def roda(comandos, timeout=300):
    entrada = '\n'.join(comandos) + '\n'
    r = subprocess.run([AVL], input=entrada, capture_output=True, text=True,
                       cwd=AQUI, timeout=timeout)
    return r.stdout


def _num(txt, padrao):
    m = re.findall(padrao, txt, flags=re.M)
    return float(m[-1]) if m else float('nan')


def le_resultados(s):
    '''Forcas totais (ft), derivadas (st) e ponto neutro do ultimo caso.'''
    return {
        'alfa': _num(s, r'Alpha\s*=\s*([-\d.]+)'),
        'CL': _num(s, r'CLtot\s*=\s*([-\d.]+)'),
        'CD': _num(s, r'CDtot\s*=\s*([-\d.]+)'),
        'CDvis': _num(s, r'CDvis\s*=\s*([-\d.Ee+]+)'),
        'CDind': _num(s, r'CDind\s*=\s*([-\d.Ee+]+)'),
        'CLff': _num(s, r'CLff\s*=\s*([-\d.]+)'),
        'CDff': _num(s, r'CDff\s*=\s*([-\d.Ee+]+)'),
        'e': _num(s, r'\se =\s+([-\d.]+)'),
        'Cm': _num(s, r'Cmtot\s*=\s*([-\d.]+)'),
        # ancorado no inicio da linha: exclui "D2  elevator -> elevator = ..."
        'de': _num(s, r'^\s+elevator\s+=\s+([-\d.]+)'),
        'it': _num(s, r'^\s+it\s+=\s+([-\d.]+)'),
        'xnp': _num(s, r'Xnp\s*=\s*([-\d.]+)'),
        'CLa': _num(s, r'CLa\s*=\s*([-\d.]+)'),
        'Cma': _num(s, r'Cma\s*=\s*([-\d.]+)'),
    }


def le_faixas(s, superficie=1):
    '''Faixas de uma superficie (bloco fs). superficie=1 e a asa direita.'''
    inicios = [m.end() for m in re.finditer(rf'Surface #\s+{superficie}\s', s)]
    if not inicios:
        raise RuntimeError(f'superficie {superficie} nao encontrada na saida do AVL')
    trecho = s[inicios[-1]:]
    fim = re.search(r'Surface #\s+\d+', trecho)
    bloco = trecho[:fim.start()] if fim else trecho
    linhas = re.findall(r'^\s*\d+((?:\s+[-\d.Ee+]+){12})\s*$', bloco, flags=re.M)
    if not linhas:
        raise RuntimeError('faixas nao encontradas na saida do AVL')
    d = np.array([[float(v) for v in ln.split()] for ln in linhas])
    return {c: d[:, i] for i, c in enumerate(COLUNAS_FS)}


def caso(arquivo, mach, alfa=None, cl=None, it=0.0, trim=False, de=0.0,
         empenagem=True, faixas=False, derivadas=False):
    '''
    Roda um caso. Informe alfa [graus] OU cl (alvo de CL).
    trim=True compensa a arfagem pelo profundor (d2 pm 0); senao o
    profundor fica em `de`. empenagem=False para o modelo so de asa.
    '''
    if (alfa is None) == (cl is None):
        raise ValueError('informe alfa ou cl')
    cmds = [f'load {arquivo}', 'oper', 'm', f'mn {mach}', '']
    if empenagem:
        cmds += ['de', f'1 {it}', '', 'd2 pm 0' if trim else f'd2 d2 {de}']
    cmds += [f'a c {cl}' if cl is not None else f'a a {alfa}', 'x']
    if faixas:
        cmds += ['fs', '']
    if derivadas:
        cmds += ['st', '']
    cmds += ['', 'quit']
    s = roda(cmds)
    r = le_resultados(s)
    if not np.isfinite(r['alfa']) or not np.isfinite(r['CL']):
        raise RuntimeError('AVL nao convergiu:\n' + s[-2000:])
    if faixas:
        r['faixas'] = le_faixas(s, superficie=1)
    r['saida'] = s
    return r
