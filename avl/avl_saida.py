'''
Leitura de blocos da saida do AVL 3.37 que o avl_run.py nao cobre, usados
pelo roteiro do Lab 04 (lab04_roteiro.py) e pelas derivadas de estabilidade
(derivadas_estabilidade.py):

  - varios casos (varios "x") na mesma sessao: um bloco de forcas totais
    por execucao (casos_ft);
  - derivadas do "st" (eixos de estabilidade) e do "sb" (eixos do corpo)
    como dicionario nome -> valor (derivadas);
  - lista de superficies do "fs" com numero, nome, malha e area integrada
    pelo AVL (superficies), e as faixas de todas elas (faixas_todas);
  - conversoes de unidade/sinal pedidas pelo roteiro (por_rad, etc.).

Sobre o "sb" do AVL 3.37: o cabecalho impresso e "Geometry-axis
derivatives...", mas os valores estao nos eixos do corpo PADRAO (X para a
frente, Z para baixo), os mesmos do "Standard axis orientation" do ft:
CXw > 0 (a sustentacao inclina para a frente com alfa), CZw < 0, Clv com o
mesmo sinal de Clb do st. O leitor aceita os dois cabecalhos ("Geometry-axis"
e "Body-axis").

Comando de arfagem no menu OPER (verificado no AVL 3.37): a variavel e "P"
(P itch rate) e a restricao tambem "P" (qc/2V), ou seja, "p p 0.005" fixa
qc/2V = 0,005. "q q ..." NAO e reconhecido ("Option not recognized").
'''

import re

import numpy as np

from avl_run import le_resultados, le_faixas, roda

GRAUS_POR_RAD = 180.0/np.pi

_NUM = r'-?\d+\.?\d*(?:[Ee][-+]?\d+)?'
_PAR = re.compile(rf'([A-Za-z][A-Za-z0-9]*)\s*=\s*({_NUM})')

CABECALHOS = {
    'st': [r'Stability-axis derivatives\.\.\.'],
    'sb': [r'Geometry-axis derivatives\.\.\.', r'Body-axis derivatives\.\.\.'],
}


def por_rad(valor_por_grau):
    '''Derivada de controle do AVL (1/grau) -> 1/rad (x 180/pi).'''
    return valor_por_grau*GRAUS_POR_RAD


def converte(valor, inverte=False, por_grau=False):
    '''Conversao do roteiro: inverte o sinal e/ou multiplica por 180/pi.'''
    v = -valor if inverte else valor
    return por_rad(v) if por_grau else v


def derivadas(s, bloco='st'):
    '''
    Derivadas do ULTIMO bloco "st" ou "sb" da saida `s`: dicionario
    nome -> valor (ex.: 'CLa', 'Cmq', 'CLd2', 'CDffd2', 'CLg1', 'Clp').
    No st o bloco termina antes da linha "Neutral point" (a linha seguinte,
    "Clb Cnr / Clr Cnb = ...", seria lida como 'Cnb'); o Xnp e lido a parte.
    '''
    inicios = []
    for padrao in CABECALHOS[bloco]:
        inicios += [m.end() for m in re.finditer(padrao, s)]
    if not inicios:
        raise RuntimeError(f'bloco "{bloco}" nao encontrado na saida do AVL')
    trecho = s[max(inicios):]
    fim = re.search(r'Neutral point|Operation of run case|Clb Cnr', trecho)
    corpo = trecho[:fim.start()] if fim else trecho
    d = {k: float(v) for k, v in _PAR.findall(corpo)}
    if bloco == 'st':
        m = re.search(rf'Neutral point\s+Xnp\s*=\s*({_NUM})', trecho)
        if m:
            d['Xnp'] = float(m.group(1))
    if not d:
        raise RuntimeError(f'bloco "{bloco}" vazio na saida do AVL')
    return d


def casos_ft(s):
    '''Um dicionario de le_resultados por bloco de forcas totais (um por "x").'''
    partes = re.split(r'Vortex Lattice Output -- Total Forces', s)[1:]
    out = []
    for p in partes:
        fim = p.find('Operation of run case')
        out.append(le_resultados(p[:fim] if fim > 0 else p))
    return out


def superficies(s):
    '''
    Superficies do ultimo bloco "fs": lista de dicionarios com numero,
    nome, n_corda, n_env (malha) e area (integrada pelo AVL, m^2).
    '''
    inicio = s.rfind('Surface and Strip Forces by surface')
    if inicio < 0:
        raise RuntimeError('bloco fs nao encontrado na saida do AVL')
    trecho = s[inicio:]
    padrao = re.compile(
        r'Surface #\s*(\d+)\s+(.+?)\s*\n\s*# Chordwise =\s*(\d+)\s+# Spanwise =\s*(\d+)'
        rf'.*?\n\s*Surface area =\s*({_NUM})', flags=re.S)
    out = []
    for m in padrao.finditer(trecho):
        out.append({'numero': int(m.group(1)), 'nome': m.group(2).strip(),
                    'n_corda': int(m.group(3)), 'n_env': int(m.group(4)),
                    'area': float(m.group(5))})
    if not out:
        raise RuntimeError('nenhuma superficie lida do bloco fs')
    return out


def faixas_todas(s):
    '''Faixas de todas as superficies do ultimo fs: {nome: faixas}.'''
    inicio = s.rfind('Surface and Strip Forces by surface')
    trecho = s[inicio:]
    return {sp['nome']: le_faixas(trecho, sp['numero']) for sp in superficies(s)}


def cabecalho_avl(caminho_abs):
    '''Mach, Sref, Cref, Bref, Xref, Yref, Zref e CDp do cabecalho do .avl.'''
    with open(caminho_abs, encoding='utf-8', errors='replace') as f:
        linhas = [ln.split('#')[0].split('!')[0].strip() for ln in f]
    uteis = [ln for ln in linhas if ln]
    # uteis[0] = titulo; depois Mach; iYsym iZsym Zsym; Sref Cref Bref;
    # Xref Yref Zref; [CDp]
    num = [[float(v) for v in ln.split()] for ln in uteis[1:6]]
    return {'Mach': num[0][0], 'Sref': num[2][0], 'Cref': num[2][1],
            'Bref': num[2][2], 'Xref': num[3][0], 'Yref': num[3][1],
            'Zref': num[3][2], 'CDp': num[4][0] if len(num[4]) == 1 else 0.0}


def secoes_superficie(caminho_abs, nome='Wing'):
    '''Secoes (Xle, Yle, Zle, corda, Ainc) de uma SURFACE do .avl.'''
    with open(caminho_abs, encoding='utf-8', errors='replace') as f:
        txt = f.read()
    m = re.search(rf'SURFACE\s*\n{re.escape(nome)}\s*\n(.*?)(?=\n#-{{10,}}|\nSURFACE|\nBODY|\Z)',
                  txt, flags=re.S)
    if not m:
        raise RuntimeError(f'superficie {nome} nao encontrada em {caminho_abs}')
    secs = re.findall(rf'SECTION\s*\n#[^\n]*\n\s*({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})',
                      m.group(1))
    return np.array([[float(v) for v in sec] for sec in secs])


def sessao(arquivo, mach, it, comandos_oper):
    '''
    Roda uma sessao do AVL: carrega `arquivo`, fixa Mach e it e executa os
    comandos do menu OPER dados. Devolve o texto da saida.
    '''
    cmds = [f'load {arquivo}', 'oper', 'm', f'mn {mach}', '',
            'de', f'1 {it}', ''] + list(comandos_oper) + ['', 'quit']
    return roda(cmds)
