# Lab 04 — posição da EH e torção da asa: plano de implementação

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Escolher a posição longitudinal da EH (`Lc_h`) por uma varredura com
restrições de fuselagem, sombreamento, profundor e margem estática. Depois,
otimizar a torção da asa nessa posição. As duas coisas vêm precedidas de
convergência de malha e de uma verificação do otimizador (asa limpa →
distribuição elíptica).

**Architecture:** Módulos Python pequenos dentro de `avl/`, de baixo para cima:
- `aeronave.py`: designTool → dicionário;
- `gera_avl.py`: dicionário → `.avl`;
- `avl_run.py`: roda o `avl.exe` e lê a saída;
- `analises.py` e `estol.py`: compensação e método da seção crítica;
- `sombra.py`: geometria pura;
- scripts de etapa: `convergencia_malha.py`, `verifica_eliptica.py`,
  `varredura_eh.py`, `otimizacao_torcao.py`.

O AVL sempre roda com `cwd = avl/` e caminhos relativos, porque a pasta do
repositório tem "º" no caminho.

**Tech Stack:** Python 3.13, numpy, scipy (SLSQP), matplotlib, pytest,
`avl/avl.exe` (3.37), designTool da disciplina.

**Especificação:** `docs/superpowers/specs/2026-09-23-lab04-posicao-eh-torcao-design.md`

**Convenções:**
- Todos os comandos rodam a partir da raiz do repositório.
- Testes: `python -m pytest avl/tests -q`.
- Commits direto na `main` (fluxo da equipe). Não fazer push sem pedir.

---

## Mapa de arquivos

| arquivo | ação | responsabilidade |
|---|---|---|
| `designTool/balance.py` | modificar (1 linha) | expor `deda` em `airplane['balance']` |
| `avl/aeronave.py` | criar | designTool → dicionário da aeronave para um `Lc_h`; limite de `Lc_h` pela fuselagem |
| `avl/gera_avl.py` | criar | escreve `.avl` (completo ou só asa) com deslocamento em Z, torção e malha |
| `avl/avl_run.py` | criar | roda o AVL por stdin e lê `ft`, `fs` e `st` |
| `avl/analises.py` | criar | `it` que zera o δe no cruzeiro |
| `avl/estol.py` | criar | método da seção crítica (modelo linear por faixa + conferência direta) |
| `avl/sombra.py` | criar | α de saída da EH da sombra da asa; fração do leme encoberta |
| `avl/estilo.py` | criar | paleta e rcParams das figuras |
| `avl/convergencia_malha.py` | criar | etapa 0a |
| `avl/otimizacao_torcao.py` | criar | avaliador e otimizador de torção (etapas 0b e 2) |
| `avl/verifica_eliptica.py` | criar | etapa 0b |
| `avl/varredura_eh.py` | criar | etapa 1 |
| `avl/tests/…` | criar | testes |
| `avl/README.md` | criar | como rodar |
| `.gitignore` | modificar | ignorar `avl/resultados/_tmp/` |
| `designTool/standard_airplane.py:171` | modificar (última tarefa) | `Lc_h` escolhido |
| `avl/fwd.avl`, `avl/aft.avl` | regerar (última tarefa) | geometria final |

---

### Task 1: expor dε/dα no designTool e criar a infraestrutura de testes

**Files:**
- Modify: `designTool/balance.py:194`
- Create: `avl/tests/conftest.py`
- Modify: `.gitignore`

- [ ] **Step 1: Adicionar a linha no `balance.py`**

Logo depois de `airplane['balance']['CLv'] = CLv` (linha 194), acrescentar:

```python
    airplane['balance']['deda'] = deda
```

- [ ] **Step 2: Criar `avl/tests/conftest.py`**

```python
import os
import sys

AVL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if AVL_DIR not in sys.path:
    sys.path.insert(0, AVL_DIR)
```

- [ ] **Step 3: Acrescentar ao `.gitignore`**

```
# Arquivos temporarios do AVL gerados pelos scripts do Lab 04
avl/resultados/_tmp/
```

- [ ] **Step 4: Conferir que o designTool continua rodando**

Run: `python -c "from designTool.standard_airplane import standard_airplane; from designTool.analyze import analyze; a=standard_airplane('my_airplane'); analyze(a); print(round(a['balance']['deda'],4))"`
Expected: um número entre 0 e 1 (algo perto de 0,3–0,5), sem erro.

- [ ] **Step 5: Commit**

```bash
git add designTool/balance.py avl/tests/conftest.py .gitignore
git commit -m "designTool expoe deda; infraestrutura de testes do Lab 04"
```

---

### Task 2: `aeronave.py` — ponte com o designTool

**Files:**
- Create: `avl/aeronave.py`
- Test: `avl/tests/test_aeronave.py`

- [ ] **Step 1: Escrever o teste**

```python
import pytest

from aeronave import aeronave, lc_h_maximo, W_PROJETO_BASE_KGF, FOLGA_FUSELAGEM
from designTool.constants import gravity


def test_base_reproduz_ponto_de_projeto_do_lab03():
    av = aeronave(4.6)
    assert av['W'] == pytest.approx(W_PROJETO_BASE_KGF*gravity, abs=1.0)
    assert av['CL'] == pytest.approx(0.5053, abs=1e-3)
    assert av['xcg_fwd'] == pytest.approx(26.0797, abs=1e-3)
    assert 0.0 < av['deda'] < 1.0


def test_alturas_sao_as_reais_do_designtool():
    av = aeronave(4.6)
    assert av['asa']['zr'] == pytest.approx(-1.3)
    assert av['EH']['zr'] == pytest.approx(2.0)
    assert av['EV']['zr'] == pytest.approx(3.0)


def test_eh_mais_para_tras_reduz_W0_e_area():
    a, b = aeronave(4.6), aeronave(4.8)
    assert b['W0'] < a['W0']
    assert b['EH']['S'] < a['EH']['S']


def test_limite_da_fuselagem():
    lc = lc_h_maximo()
    av = aeronave(lc)
    fim_eh = av['EH']['xr'] + av['EH']['cr']
    assert fim_eh == pytest.approx(av['fuselagem']['L'] - FOLGA_FUSELAGEM, abs=1e-3)
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `python -m pytest avl/tests/test_aeronave.py -q`
Expected: FAIL com `ModuleNotFoundError: No module named 'aeronave'`

- [ ] **Step 3: Implementar `avl/aeronave.py`**

```python
'''
Ponte entre o designTool e o Lab 04.

aeronave(Lc_h) roda o designTool com o braco da EH pedido e devolve um
dicionario com tudo o que os scripts do AVL usam: ponto de projeto,
referencias, CGs, dε/dα e a geometria de asa, EH, EV, nacele e fuselagem.
As alturas (z) sao as REAIS do designTool; o deslocamento em Z exigido
no AVL e aplicado so no gera_avl.py.

O peso no ponto de projeto segue a definicao do Lab 03 (peso medio de
cruzeiro, 229.669,3 kgf com Lc_h = 4,6). Para outros Lc_h mantem-se a mesma
fracao de combustivel, com 100% de carga paga, e o peso e recalculado com o
W_empty e o W_fuel do designTool.
'''

import functools
import os
import sys

import numpy as np

RAIZ_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if RAIZ_REPO not in sys.path:
    sys.path.insert(0, RAIZ_REPO)

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.geometry import geometry
from designTool.aerodynamics import aerodynamics
from designTool.auxiliary import atmosphere
from designTool.constants import gravity

LC_H_BASE = 4.6
W_PROJETO_BASE_KGF = 229669.3
FOLGA_FUSELAGEM = 0.5      # [m] bordo de fuga da raiz da EH ate o fim da fuselagem


def _analisa(Lc_h):
    ap = standard_airplane('my_airplane')
    ap['inputs']['Lc_h'] = float(Lc_h)
    analyze(ap, print_log=False, plot=False)
    return ap


@functools.lru_cache(maxsize=None)
def fracao_combustivel():
    '''Fracao de combustivel que reproduz o peso de projeto do Lab 03.'''
    ap = _analisa(LC_H_BASE)
    tm, inp = ap['thrust_matching'], ap['inputs']
    W = W_PROJETO_BASE_KGF*gravity
    return (W - tm['W_empty'] - inp['W_payload'] - inp['W_crew'])/tm['W_fuel']


@functools.lru_cache(maxsize=None)
def aeronave(Lc_h=LC_H_BASE):
    '''Dicionario da aeronave para um Lc_h. Nao modifique o retorno (cache).'''
    ap = _analisa(Lc_h)
    inp, geo = ap['inputs'], ap['geometry']
    tm, bal = ap['thrust_matching'], ap['balance']

    h, M = inp['altitude_cruise'], inp['Mach_cruise']
    atm = atmosphere(h)
    rho, a_inf = atm['density'], atm['speed_of_sound']
    V = M*a_inf
    q = 0.5*rho*V**2
    S_w = inp['S_w']

    W = (tm['W_empty'] + inp['W_payload'] + inp['W_crew']
         + fracao_combustivel()*tm['W_fuel'])
    CL = W/(q*S_w)
    _, _, drag = aerodynamics(ap, Mach=M, altitude=h, CL=CL)

    return {
        'Lc_h': float(Lc_h),
        'W0': tm['W0'], 'W': W, 'W_kgf': W/gravity, 'W0_kgf': tm['W0']/gravity,
        'fuel_frac': fracao_combustivel(),
        'M': M, 'h': h, 'rho': rho, 'a': a_inf, 'V': V, 'CL': CL,
        'CD0': drag['CD0'],
        'Sref': S_w, 'Cref': geo['cm_w'], 'Bref': geo['b_w'],
        'xcg_fwd': bal['xcg_fwd'], 'xcg_aft': bal['xcg_aft'],
        'xnp_dt': bal['xnp'], 'deda': bal['deda'],
        'asa': {'xr': inp['xr_w'], 'zr': inp['zr_w'], 'cr': geo['cr_w'],
                'ct': geo['ct_w'], 'xt': geo['xt_w'], 'yt': geo['yt_w'],
                'zt': geo['zt_w']},
        'EH': {'S': geo['S_h'], 'xr': geo['xr_h'], 'zr': inp['zr_h'],
               'cr': geo['cr_h'], 'ct': geo['ct_h'], 'xt': geo['xt_h'],
               'yt': geo['yt_h'], 'zt': geo['zt_h'], 'tc': inp['tcr_h']},
        'EV': {'xr': geo['xr_v'], 'zr': inp['zr_v'], 'cr': geo['cr_v'],
               'ct': geo['ct_v'], 'xt': geo['xt_v'], 'zt': geo['zt_v']},
        'nacele': {'x': inp['x_n'], 'y': inp['y_n'], 'z': inp['z_n'],
                   'L': inp['L_n'], 'D': inp['D_n']},
        'fuselagem': {'L': inp['L_f'], 'D': inp['D_f']},
    }


def _fim_eh(Lc_h):
    ap = standard_airplane('my_airplane')
    ap['inputs']['Lc_h'] = float(Lc_h)
    geometry(ap)
    return ap['geometry']['xr_h'] + ap['geometry']['cr_h'], ap['inputs']['L_f']


def lc_h_maximo(folga=FOLGA_FUSELAGEM, lo=4.0, hi=6.0, tol=1e-5):
    '''Maior Lc_h com o bordo de fuga da raiz da EH a `folga` do fim da fuselagem.'''
    while hi - lo > tol:
        mid = 0.5*(lo + hi)
        fim, L_f = _fim_eh(mid)
        if fim <= L_f - folga:
            lo = mid
        else:
            hi = mid
    return lo
```

- [ ] **Step 4: Rodar e ver passar**

Run: `python -m pytest avl/tests/test_aeronave.py -q`
Expected: `4 passed`. Se `test_base_reproduz_ponto_de_projeto_do_lab03`
falhar em `xcg_fwd`, é porque o `balance.py` mudou: confira com
`git diff designTool/`.

- [ ] **Step 5: Commit**

```bash
git add avl/aeronave.py avl/tests/test_aeronave.py
git commit -m "Lab 04: ponte com o designTool parametrizada por Lc_h"
```

---

### Task 3: `avl_run.py` — rodar o AVL e ler a saída

**Files:**
- Create: `avl/avl_run.py`
- Create: `avl/tests/dados/amostra_avl.txt` (saída real capturada)
- Test: `avl/tests/test_avl_run.py`

- [ ] **Step 1: Capturar a amostra de saída**

```bash
mkdir -p avl/tests/dados
cd avl && printf 'load fwd.avl\noper\nm\nmn 0.85\n\nde\n1 0\n\nd2 pm 0\na c 0.5\nx\nfs\n\nst\n\n\nquit\n' | ./avl.exe > tests/dados/amostra_avl.txt; cd ..
```

Conferir: `grep -n "CDff\|Xnp" avl/tests/dados/amostra_avl.txt` mostra
`CDff  = 0.0113659` e `Xnp =  27.427111`.

- [ ] **Step 2: Escrever o teste**

```python
import os

import numpy as np
import pytest

from avl_run import le_resultados, le_faixas, caso

AMOSTRA = os.path.join(os.path.dirname(__file__), 'dados', 'amostra_avl.txt')


def texto():
    with open(AMOSTRA) as f:
        return f.read()


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
```

- [ ] **Step 3: Rodar e ver falhar**

Run: `python -m pytest avl/tests/test_avl_run.py -q`
Expected: FAIL com `No module named 'avl_run'`

- [ ] **Step 4: Implementar `avl/avl_run.py`**

```python
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
```

- [ ] **Step 5: Rodar e ver passar**

Run: `python -m pytest avl/tests/test_avl_run.py -q`
Expected: `3 passed`

- [ ] **Step 6: Commit**

```bash
git add avl/avl_run.py avl/tests/test_avl_run.py avl/tests/dados/amostra_avl.txt
git commit -m "Lab 04: execucao do AVL e leitura de ft, fs e st"
```

---

### Task 4: `gera_avl.py` — escrever o `.avl` a partir do designTool

**Files:**
- Create: `avl/gera_avl.py`
- Create: `avl/tests/dados/fwd_referencia_11c6802.avl`
- Test: `avl/tests/test_gera_avl.py`

- [ ] **Step 1: Congelar a referência**

```bash
git show 11c6802:avl/fwd.avl > avl/tests/dados/fwd_referencia_11c6802.avl
```

- [ ] **Step 2: Escrever o teste**

```python
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


def test_regenera_o_modelo_de_referencia(tmp_path):
    escreve_avl(aeronave(4.6), 'resultados/_tmp/teste_fwd.avl', cg='fwd',
                malha=MALHA_PADRAO)
    for a in (0.0, 4.0):
        novo = caso('resultados/_tmp/teste_fwd.avl', 0.85, alfa=a)
        ref = caso(REF, 0.85, alfa=a)
        for k in ('CL', 'Cm', 'CDff', 'CD'):
            assert novo[k] == pytest.approx(ref[k], abs=2e-5), k


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
```

- [ ] **Step 3: Rodar e ver falhar**

Run: `python -m pytest avl/tests/test_gera_avl.py -q`
Expected: FAIL com `No module named 'gera_avl'`

- [ ] **Step 4: Implementar `avl/gera_avl.py`**

```python
'''
Escreve os arquivos do AVL a partir do dicionario de aeronave.py.

Deslocamentos em Z (exigencia do professor, para manter distancia entre os
paineis): asa -1,20 m, EH +1,85 m, EV +0,85 m, nacele acompanhando a asa.
Nao sao alteracao de projeto e existem so aqui.

Winglet: superficie propria no COMPONENT 1 (mesmo componente da asa), com
toe zero; geometria do designTool (Raymer, Fig. 7.34): vertical, altura e
corda de raiz iguais a corda da ponta, afilamento 0,21, bordo de fuga reto.

Nacele: anel sustentador no COMPONENT 1, como no 737.avl da disciplina.
'''

import json
import os

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))

DZ = {'asa': -1.20, 'EH': +1.85, 'EV': +0.85}
TAPER_WINGLET = 0.21

# Estacoes da asa (fracao da semienvergadura) e perfis do Lab 03
ETAS = [0.0, 0.1011, 0.22, 0.32, 0.398, 0.48, 0.56, 0.70, 0.82, 0.90, 1.0]
PERFIS = ['raiz.dat', 'raiz.dat', 'gerado/eta_0220.dat', 'gerado/eta_0320.dat',
          'gerado/eta_0398.dat', 'gerado/eta_0480.dat', 'gerado/eta_0560.dat',
          'gerado/eta_0700.dat', 'gerado/eta_0820.dat', 'ponta.dat',
          'ponta.dat']
ETA_AILERON = (0.56, 0.90)

MALHA_PADRAO = {'asa_nc': 8, 'asa_ns': 20, 'winglet_ns': 8, 'eh_nc': 16,
                'eh_ns': 10, 'ev_nc': 8, 'ev_ns': 6, 'nbody': 20}
ARQ_MALHA = os.path.join(AQUI, 'resultados', 'malha_adotada.json')

_ANEL = [(0.0, 1.0), (0.5, 0.866), (0.866, 0.5), (1.0, 0.0), (0.866, -0.5),
         (0.5, -0.866), (0.0, -1.0), (-0.5, -0.866), (-0.866, -0.5),
         (-1.0, 0.0), (-0.866, 0.5), (-0.5, 0.866), (0.0, 1.0)]


def malha_adotada():
    '''Malha da convergencia (resultados/malha_adotada.json) ou a padrao.'''
    m = dict(MALHA_PADRAO)
    if os.path.exists(ARQ_MALHA):
        with open(ARQ_MALHA) as f:
            m.update(json.load(f))
    return m


def perfil(eta):
    i = int(np.argmin([abs(eta - e) for e in ETAS]))
    return 'airfoils/' + PERFIS[i]


def secoes_asa(av, torcao=None, etas=ETAS, diedro=True):
    w = av['asa']
    etas = list(etas)
    torcao = np.zeros(len(etas)) if torcao is None else np.asarray(torcao, float)
    if len(torcao) != len(etas):
        raise ValueError('torcao precisa ter um valor por estacao')
    out = []
    for eta, tw in zip(etas, torcao):
        z_real = w['zr'] + eta*(w['zt'] - w['zr']) if diedro else w['zr']
        out.append({'eta': eta, 'y': eta*w['yt'],
                    'x': w['xr'] + eta*(w['xt'] - w['xr']),
                    'z': z_real + DZ['asa'],
                    'c': w['cr'] + eta*(w['ct'] - w['cr']),
                    'ainc': float(tw), 'perfil': perfil(eta),
                    'aileron': ETA_AILERON[0] - 1e-9 <= eta <= ETA_AILERON[1] + 1e-9})
    return out


def _secao(x, y, z, c, ainc=0.0):
    return ('\nSECTION\n#Xle      Yle      Zle      Chord    Ainc\n'
            f'{x:.4f}  {y:.4f}  {z:.4f}  {c:.4f}  {ainc:.4f}\n')


def _superficie(nome, nc, ns, extra=''):
    return ('#' + '-'*64 + f'\nSURFACE\n{nome}\n#Nchordwise Cspace Nspanwise Sspace\n'
            f'{nc} 1.0 {ns} 1.0\n{extra}')


def _asa(secoes, malha, controles):
    txt = _superficie('Wing', malha['asa_nc'], malha['asa_ns'],
                      'YDUPLICATE\n0.0\nANGLE\n0.0\n')
    for s in secoes:
        txt += _secao(s['x'], s['y'], s['z'], s['c'], s['ainc'])
        txt += f'AFILE\n{s["perfil"]}\nCLAF\n1.0000\n'
        if controles and s['aileron']:
            txt += 'CONTROL\naileron  1.0   0.73    0. 0. 0.    -1.0\n'
    return txt


def _winglet(secoes, malha):
    p = secoes[-1]
    ct = p['c']
    txt = _superficie('Winglet', malha['asa_nc'], malha['winglet_ns'],
                      'COMPONENT\n1\nYDUPLICATE\n0.0\nANGLE\n0.0\n')
    txt += _secao(p['x'], p['y'], p['z'], ct) + 'AFILE\nairfoils/ponta.dat\nCLAF\n1.0000\n'
    txt += (_secao(p['x'] + ct - TAPER_WINGLET*ct, p['y'], p['z'] + ct,
                   TAPER_WINGLET*ct) + 'AFILE\nairfoils/ponta.dat\nCLAF\n1.0000\n')
    return txt


def _eh(av, malha):
    h = av['EH']
    txt = _superficie('Horizontal tail', malha['eh_nc'], malha['eh_ns'],
                      'YDUPLICATE\n0.0\nANGLE\n0.0\n')
    for x, y, z, c in ((h['xr'], 0.0, h['zr'], h['cr']),
                       (h['xt'], h['yt'], h['zt'], h['ct'])):
        txt += _secao(x, y, z + DZ['EH'], c)
        txt += ('NACA\n0010\nCLAF\n1.0000\nDESIGN\nit      1.0\n'
                'CONTROL\nelevator 1.0   0.70    0. 0. 0.     1.0\n')
    return txt


def _ev(av, malha):
    v = av['EV']
    txt = _superficie('Vertical tail', malha['ev_nc'], malha['ev_ns'], 'ANGLE\n0.0\n')
    for x, z, c in ((v['xr'], v['zr'], v['cr']), (v['xt'], v['zt'], v['ct'])):
        txt += _secao(x, 0.0, z + DZ['EV'], c)
        txt += ('NACA\n0010\nCLAF\n1.0000\n'
                'CONTROL\nrudder   1.0   0.70    0. 0. 0.     0.0\n')
    return txt


def _corpo(av, malha):
    f = av['fuselagem']
    return ('#' + '-'*64 + f'\nBODY\nFuselage\n# Nbody Bspace\n{malha["nbody"]} 1.0\n'
            f'SCALE\n{f["L"]:.4f} {f["D"]:.4f} {f["D"]:.4f}\n'
            'BFILE\nfuselage_nondim.dat\n')


def _nacele(av):
    n = av['nacele']
    txt = ('#' + '-'*64 + '\nSURFACE\nNacelle\n#Nchordwise  Cspace   Nspanwise  Sspace\n'
           '6            1.0      12          0.0\nCOMPONENT\n1\nYDUPLICATE\n0.0\n'
           f'SCALE\n{n["L"]:.4f}  {n["D"]/2:.4f}  {n["D"]/2:.4f}\n'
           f'TRANSLATE\n{n["x"]:.4f}  {n["y"]:.4f}  {n["z"] + DZ["asa"]:.4f}\n')
    for y, z in _ANEL:
        txt += ('\nSECTION\n#Xle   Yle    Zle      Chord   Ainc  Nspanwise  Sspace\n'
                f' 0.00  {y:.3f}  {z:.3f}  1.0  0.  1  0.\n')
    return txt


def escreve_avl(av, caminho, cg='aft', torcao=None, malha=None, etas=ETAS,
                so_asa=False, diedro=True, cdp=None):
    '''
    Escreve o .avl em avl/<caminho> e devolve `caminho` (relativo a avl/).
    so_asa=True: so a asa, sem winglet, controles, corpos ou CDp
    (verificacao eliptica).
    '''
    malha = malha_adotada() if malha is None else malha
    secoes = secoes_asa(av, torcao, etas, diedro)
    xref = av['xcg_fwd'] if cg == 'fwd' else av['xcg_aft']
    cdp = (0.0 if so_asa else av['CD0']) if cdp is None else cdp

    txt = (f'Ararinha Lc_h={av["Lc_h"]:.4f} CG {cg}{" so asa" if so_asa else ""}\n'
           f'{av["M"]:.4f}          # Mach\n0 0 0.0       # iYsym iZsym Zsym\n'
           f'{av["Sref"]:.4f} {av["Cref"]:.4f} {av["Bref"]:.4f}   # Sref Cref Bref\n'
           f'{xref:.4f} 0.0 0.0   # Xref Yref Zref\n{cdp:.5f}       # CDp\n')
    txt += _asa(secoes, malha, controles=not so_asa)
    if not so_asa:
        txt += _winglet(secoes, malha) + _eh(av, malha) + _ev(av, malha)
        txt += _corpo(av, malha) + _nacele(av)

    destino = os.path.join(AQUI, caminho)
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    with open(destino, 'w', encoding='utf-8') as f:
        f.write(txt)
    return caminho
```

- [ ] **Step 5: Rodar e ver passar**

Run: `python -m pytest avl/tests/test_gera_avl.py -q`
Expected: `4 passed`. Se a regressão falhar por mais de 2e-5, compare
seção por seção com `diff` entre `avl/resultados/_tmp/teste_fwd.avl` e a
referência. A causa mais provável é o CDp (CD0 do designTool com CL
ligeiramente diferente), que só afeta `CD`.

- [ ] **Step 6: Commit**

```bash
git add avl/gera_avl.py avl/tests/test_gera_avl.py avl/tests/dados/fwd_referencia_11c6802.avl
git commit -m "Lab 04: gerador do .avl a partir do designTool, conferido contra o modelo de referencia"
```

---

### Task 5: `analises.py` e `estol.py` — compensação e seção crítica

**Files:**
- Create: `avl/analises.py`
- Create: `avl/estol.py`
- Test: `avl/tests/test_analises.py`

- [ ] **Step 1: Escrever o teste**

```python
import numpy as np
import pytest

from aeronave import aeronave
from analises import it_para_de_zero
from estol import clmax_local, estol, ETA_CLMAX, CLMAX
from gera_avl import escreve_avl


def test_clmax_local_interpola_e_satura():
    assert clmax_local(0.0) == pytest.approx(CLMAX[0])
    assert clmax_local(1.0) == pytest.approx(CLMAX[-1])
    meio = 0.5*(ETA_CLMAX[0] + ETA_CLMAX[1])
    assert clmax_local(meio) == pytest.approx(0.5*(CLMAX[0] + CLMAX[1]))


def test_it_zera_o_profundor():
    av = aeronave(4.6)
    arq = escreve_avl(av, 'resultados/_tmp/teste_an.avl', cg='aft')
    it, r = it_para_de_zero(arq, av['M'], av['CL'])
    assert abs(r['de']) < 0.01
    assert r['CL'] == pytest.approx(av['CL'], abs=1e-4)
    assert np.isfinite(r['xnp'])


def test_estol_modelo_linear_bate_com_a_conferencia_direta():
    av = aeronave(4.6)
    arq = escreve_avl(av, 'resultados/_tmp/teste_an.avl', cg='fwd')
    e = estol(arq, it=0.0, trim=True, semi=av['asa']['yt'])
    assert abs(e['excesso_max']) < 0.02     # a faixa critica chega ao clmax
    assert 5.0 < e['alfa'] < 25.0
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `python -m pytest avl/tests/test_analises.py -q`
Expected: FAIL com `No module named 'analises'`

- [ ] **Step 3: Implementar `avl/analises.py`**

```python
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
```

- [ ] **Step 4: Implementar `avl/estol.py`**

```python
'''
Metodo da secao critica.

O cl_norm de cada faixa da asa e linear em alfa no VLM (inclusive com o
profundor de compensacao, que tambem varia linearmente). Duas rodadas dao o
alfa em que cada faixa alcanca o clmax do seu perfil; a menor define o
estol. Uma terceira rodada, no alfa encontrado, confere o resultado.

clmax dos perfis do Lab 03 (XFoil corrigido, otimizacao_aerofolio/campanha/
resultados/revisao_plano_clmax.csv, plano "normal"), interpolado entre as
estacoes de referencia e constante fora delas.
'''

import numpy as np

from avl_run import caso

ETA_CLMAX = [0.1011, 0.398, 0.90]
CLMAX = [1.774, 1.7985, 1.7338]
MACH_BAIXO = 0.2
ALFAS_BASE = (8.0, 14.0)


def clmax_local(eta):
    return np.interp(eta, ETA_CLMAX, CLMAX)


def modelo_estol(arquivo, it, trim, semi, empenagem=True):
    '''alfa de estol de cada faixa da asa direita (modelo linear).'''
    a0, a1 = ALFAS_BASE
    r0, r1 = (caso(arquivo, MACH_BAIXO, alfa=a, it=it, trim=trim,
                   empenagem=empenagem, faixas=True) for a in ALFAS_BASE)
    eta = r0['faixas']['y']/semi
    c0, c1 = r0['faixas']['cl_norm'], r1['faixas']['cl_norm']
    b = (c1 - c0)/(a1 - a0)
    a = c0 - b*a0
    return {'eta': eta, 'alfa_i': (clmax_local(eta) - a)/b}


def estol(arquivo, it, trim, semi, empenagem=True):
    m = modelo_estol(arquivo, it, trim, semi, empenagem)
    i = int(np.argmin(m['alfa_i']))
    alfa = float(m['alfa_i'][i])
    r = caso(arquivo, MACH_BAIXO, alfa=alfa, it=it, trim=trim,
             empenagem=empenagem, faixas=True)
    eta = r['faixas']['y']/semi
    lim = clmax_local(eta)
    return {'alfa': alfa, 'CL': r['CL'], 'de': r['de'] if trim else 0.0,
            'eta_crit': float(m['eta'][i]), 'eta': eta,
            'cl_norm': r['faixas']['cl_norm'], 'clmax': lim,
            'excesso_max': float(np.max(r['faixas']['cl_norm'] - lim)),
            'alfa_i': m['alfa_i']}
```

- [ ] **Step 5: Rodar e ver passar**

Run: `python -m pytest avl/tests/test_analises.py -q`
Expected: `3 passed`

- [ ] **Step 6: Commit**

```bash
git add avl/analises.py avl/estol.py avl/tests/test_analises.py
git commit -m "Lab 04: it de compensacao e metodo da secao critica"
```

---

### Task 6: `sombra.py` — sombreamento da EH e leme encoberto

**Files:**
- Create: `avl/sombra.py`
- Test: `avl/tests/test_sombra.py`

- [ ] **Step 1: Escrever o teste (geometria construída à mão)**

```python
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
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `python -m pytest avl/tests/test_sombra.py -q`
Expected: FAIL com `No module named 'sombra'`

- [ ] **Step 3: Implementar `avl/sombra.py`**

```python
'''
Sombreamento, com geometria REAL do designTool (sem o deslocamento em Z do
AVL, cuja esteira e reta e nao enxerga o fenomeno).

Sombra da asa sobre a EH: faixa entre as retas que partem do bordo de
ataque e do bordo de fuga de cada secao da asa na direcao do escoamento
local, inclinada de gama = alfa - eps = (1 - dε/dα)·alfa em relacao ao eixo
do corpo (raiz sem incidencia, eps0 desprezado). A EH esta fora da sombra
quando o extradorso do seu bordo de ataque fica abaixo da reta do bordo de
fuga. Como gama cresce com alfa, a reta so sobe: depois de sair, a EH nao
volta.

Leme encoberto (recuperacao de parafuso): regiao entre a reta a 60 graus
pelo bordo de ataque da raiz da EH e a reta a 30 graus pelo bordo de fuga,
as mesmas linhas tracejadas do plots.py do designTool. So informativo.
'''

import numpy as np


def alfa_saida(av, n=41):
    '''alfa [graus] a partir do qual toda a EH fica abaixo da sombra da asa.'''
    w, h = av['asa'], av['EH']
    y = np.linspace(0.0, h['yt'], n)
    s_h = y/h['yt']
    x_le_h = h['xr'] + s_h*(h['xt'] - h['xr'])
    c_h = h['cr'] + s_h*(h['ct'] - h['cr'])
    z_topo = h['zr'] + s_h*(h['zt'] - h['zr']) + h['tc']*c_h/2

    s_w = y/w['yt']
    x_te_w = w['xr'] + s_w*(w['xt'] - w['xr']) + w['cr'] + s_w*(w['ct'] - w['cr'])
    z_w = w['zr'] + s_w*(w['zt'] - w['zr'])

    gama = np.arctan2(z_topo - z_w, x_le_h - x_te_w).max()
    return float(np.degrees(gama)/(1.0 - av['deda']))


def fracao_leme_encoberto(av, n=400, charneira=0.70):
    h, v = av['EH'], av['EV']
    z = np.linspace(v['zr'], v['zt'], n)
    s = (z - v['zr'])/(v['zt'] - v['zr'])
    x_le = v['xr'] + s*(v['xt'] - v['xr'])
    c = v['cr'] + s*(v['ct'] - v['cr'])
    x_ch, x_te = x_le + charneira*c, x_le + c

    dz = z - h['zr']
    lim1 = h['xr'] + dz/np.tan(np.radians(60.0))
    lim2 = h['xr'] + h['cr'] + dz/np.tan(np.radians(30.0))
    sobre = np.clip(np.minimum(x_te, lim2) - np.maximum(x_ch, lim1), 0.0, None)
    sobre[dz < 0] = 0.0
    return float(sobre.sum()/(x_te - x_ch).sum())
```

- [ ] **Step 4: Rodar e ver passar**

Run: `python -m pytest avl/tests/test_sombra.py -q`
Expected: `5 passed`

- [ ] **Step 5: Conferir contra a aeronave real**

Run: `cd avl && python -c "from aeronave import aeronave; from sombra import *; av=aeronave(4.6); print(round(alfa_saida(av),2), round(fracao_leme_encoberto(av),3))"`
Expected: α de saída perto de 10–12° (a estimativa da especificação era
≈ 11°) e fração entre 0 e 1.

- [ ] **Step 6: Commit**

```bash
git add avl/sombra.py avl/tests/test_sombra.py
git commit -m "Lab 04: criterio de sombreamento da EH pela asa e leme encoberto"
```

---

### Task 7: `estilo.py` e `convergencia_malha.py` (etapa 0a)

**Files:**
- Create: `avl/estilo.py`
- Create: `avl/convergencia_malha.py`

- [ ] **Step 1: Criar `avl/estilo.py`**

```python
'''Estilo comum das figuras: marcas finas, grade discreta, paleta fixa.'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

AZUL, LARANJA, VERDE = '#2a78d6', '#eb6834', '#1baf7a'
TINTA, TINTA2, GRADE, CINZA = '#0b0b0b', '#52514e', '#e1e0d9', '#b9b8b3'

plt.rcParams.update({
    'figure.dpi': 150, 'savefig.dpi': 200, 'savefig.bbox': 'tight',
    'axes.grid': True, 'grid.color': GRADE, 'grid.linewidth': 0.7,
    'axes.edgecolor': CINZA, 'axes.labelcolor': TINTA2,
    'xtick.color': TINTA2, 'ytick.color': TINTA2, 'text.color': TINTA,
    'axes.spines.top': False, 'axes.spines.right': False,
    'lines.linewidth': 2.0, 'lines.markersize': 5, 'font.size': 9,
    'legend.frameon': False,
})
```

- [ ] **Step 2: Criar `avl/convergencia_malha.py`**

```python
'''
Etapa 0a -- convergencia de malha.

Geometria de base (Lc_h = 4,6, sem torcao, CG traseiro), ponto de projeto
compensado com it que zera o profundor. Refina uma direcao por vez, na
ordem de SEQ, ja usando o que foi adotado nas anteriores. Adota a menor malha
cuja diferenca para a mais fina da sequencia fica abaixo de TOL.

Rodar de dentro de avl/:
    python convergencia_malha.py            # estudo completo
    python convergencia_malha.py --lc 4.8   # confere a malha adotada em outro Lc_h
'''

import argparse
import json
import os

import numpy as np

from estilo import plt, AZUL, LARANJA, TINTA2
from aeronave import aeronave, LC_H_BASE
from analises import it_para_de_zero
from estol import estol
from gera_avl import escreve_avl, MALHA_PADRAO, ARQ_MALHA, malha_adotada, AQUI

SEQ = [('asa_ns', [10, 20, 30, 40]), ('asa_nc', [4, 8, 12, 16]),
       ('winglet_ns', [4, 8, 12]), ('eh_ns', [5, 10, 15, 20]),
       ('eh_nc', [8, 16, 24]), ('nbody', [10, 20, 40])]
TOL = {'CDff': 0.1, 'it': 0.02}          # [count], [graus]
ARQ = 'resultados/_tmp/malha.avl'


def metricas(av, malha):
    escreve_avl(av, ARQ, cg='aft', malha=malha)
    it, r = it_para_de_zero(ARQ, av['M'], av['CL'])
    e = estol(ARQ, it=it, trim=True, semi=av['asa']['yt'])
    return {'CDff': 1e4*r['CDff'], 'it': it, 'alfa_estol': e['alfa']}


def estudo():
    av = aeronave(LC_H_BASE)
    malha = dict(MALHA_PADRAO)
    historico = {}
    for chave, valores in SEQ:
        res = []
        for n in valores:
            m = dict(malha, **{chave: n})
            res.append(metricas(av, m))
            print(f'  {chave:10s} = {n:3d}   CDff = {res[-1]["CDff"]:8.3f} count'
                  f'   it = {res[-1]["it"]:7.3f}   alfa_estol = {res[-1]["alfa_estol"]:6.2f}')
        fino = res[-1]
        adotado = next(n for n, r in zip(valores, res)
                       if abs(r['CDff'] - fino['CDff']) < TOL['CDff']
                       and abs(r['it'] - fino['it']) < TOL['it'])
        malha[chave] = adotado
        historico[chave] = {'valores': valores, 'res': res, 'adotado': adotado}
        print(f'  -> {chave} adotado = {adotado}\n')

    os.makedirs(os.path.dirname(ARQ_MALHA), exist_ok=True)
    with open(ARQ_MALHA, 'w') as f:
        json.dump(malha, f, indent=2)
    with open(os.path.join(AQUI, 'resultados', 'convergencia_malha.json'), 'w') as f:
        json.dump(historico, f, indent=2)
    figura(historico)
    print('malha adotada:', malha)


def figura(hist):
    fig, eixos = plt.subplots(2, len(hist), figsize=(2.3*len(hist), 4.2), sharex='col')
    for j, (chave, h) in enumerate(hist.items()):
        v = h['valores']
        for i, (k, cor, rot, tol) in enumerate((('CDff', AZUL, 'CDff [count]', TOL['CDff']),
                                                 ('it', LARANJA, 'it [°]', TOL['it']))):
            y = [r[k] for r in h['res']]
            ax = eixos[i, j]
            ax.axhspan(y[-1] - tol, y[-1] + tol, color='#e1e0d9', lw=0)
            ax.plot(v, y, '-o', color=cor)
            ax.axvline(h['adotado'], color=TINTA2, lw=0.8, ls='--')
            if j == 0:
                ax.set_ylabel(rot)
        eixos[1, j].set_xlabel(chave)
    fig.suptitle('Convergência de malha (faixa cinza: tolerância em torno da malha mais fina)')
    fig.savefig(os.path.join(AQUI, 'resultados', 'convergencia_malha.png'))


def confere(lc):
    av = aeronave(lc)
    m = malha_adotada()
    base = metricas(av, m)
    fina = metricas(av, dict(m, asa_ns=2*m['asa_ns'], eh_ns=2*m['eh_ns']))
    d = {k: fina[k] - base[k] for k in base}
    ok = abs(d['CDff']) < TOL['CDff'] and abs(d['it']) < TOL['it']
    print(f'Lc_h = {lc}: dCDff = {d["CDff"]:+.3f} count, dit = {d["it"]:+.4f} graus'
          f' -> {"ok" if ok else "REVER MALHA"}')
    return ok


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--lc', type=float, default=None)
    a = p.parse_args()
    confere(a.lc) if a.lc is not None else estudo()
```

- [ ] **Step 3: Rodar o estudo**

Run: `cd avl && python convergencia_malha.py`
Expected: tabela por direção e uma linha `-> <chave> adotado = N` para cada
uma. Gera `resultados/malha_adotada.json`, `resultados/convergencia_malha.json`
e `resultados/convergencia_malha.png`. Abra o PNG e confira que as curvas
entram na faixa cinza. Se `asa_ns` = 40 for adotado (a última), estenda a
sequência para `[10, 20, 30, 40, 60]` e rode de novo: o platô tem de ficar
dentro da sequência.

- [ ] **Step 4: Rodar a suíte inteira (a malha nova não pode quebrar nada)**

Run: `python -m pytest avl/tests -q`
Expected: todos passam. A regressão da Task 4 usa `MALHA_PADRAO`
explicitamente, então não depende da malha adotada.

- [ ] **Step 5: Commit**

```bash
git add avl/estilo.py avl/convergencia_malha.py avl/resultados/malha_adotada.json avl/resultados/convergencia_malha.json avl/resultados/convergencia_malha.png
git commit -m "Lab 04: convergencia de malha e malha adotada"
```

---

### Task 8: `otimizacao_torcao.py` (avaliador) e `verifica_eliptica.py` (etapa 0b)

**Files:**
- Create: `avl/otimizacao_torcao.py` (só o avaliador e `otimiza`; o `main` vem na Task 10)
- Create: `avl/verifica_eliptica.py`
- Test: `avl/tests/test_otimizacao.py`

- [ ] **Step 1: Escrever o teste**

```python
import numpy as np
import pytest

from aeronave import aeronave
from otimizacao_torcao import Avaliador, elipse_normalizada


def test_elipse_normalizada_integra_o_CL():
    # integral de sqrt(1-eta^2) em [0,1] = pi/4 -> area total = Sref*CL
    eta = np.linspace(0, 1, 20001)
    assert np.trapezoid(np.sqrt(1 - eta**2), eta) == pytest.approx(np.pi/4, rel=1e-4)


def test_avaliador_asa_limpa_tem_cache_e_bate_o_CL():
    av = aeronave(4.6)
    a = Avaliador(av, modelo='asa_limpa', com_estol=False, nome='teste_ot')
    x = np.zeros(10)
    f1 = a.objetivo(x)
    n = a.n_rodadas
    f2 = a.objetivo(x.copy())
    assert f1 == f2 and a.n_rodadas == n          # segunda chamada vem do cache
    assert a.ultimo['r']['CL'] == pytest.approx(av['CL'], abs=1e-4)
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `python -m pytest avl/tests/test_otimizacao.py -q`
Expected: FAIL com `No module named 'otimizacao_torcao'`

- [ ] **Step 3: Implementar a parte reutilizável de `avl/otimizacao_torcao.py`**

```python
'''
Otimizacao de torcao da asa.

Variaveis: Ainc das estacoes da asa, com a raiz fixa em 0 (a raiz define o
eixo do corpo). Objetivo: CDff (Trefftz) no CL de projeto, M = 0,85 -- o
CDp e constante, entao minimizar CDff e minimizar o CD.

modelo='completo': aft.avl com compensacao pelo profundor (d2 pm 0) e it
fixo; restricao de estol: toda faixa externa a ETA_ESTOL tem de estolar
pelo menos DESEMPATE graus depois da primeira faixa interna.
modelo='asa_limpa': so a asa, plana, sem winglet, sem compensacao
(verificacao eliptica).

O mesmo codigo serve as duas coisas, de proposito: a verificacao eliptica
testa exatamente o otimizador usado no projeto.
'''

import json
import os

import numpy as np
from scipy.optimize import minimize

from avl_run import caso
from estol import modelo_estol
from gera_avl import escreve_avl, ETAS, AQUI

LIMITES = (-8.0, 3.0)
ETA_ESTOL = 0.56        # raiz do aileron: o estol tem de comecar para dentro
DESEMPATE = 0.2         # [graus]
EPS_FD = 0.2            # [graus] passo das diferencas finitas


def elipse_normalizada(eta):
    return np.sqrt(np.clip(1.0 - np.asarray(eta)**2, 0.0, None))


def carga_normalizada(av, r):
    '''c·cl dividido pela carga elíptica na raiz, 4 Sref CL / (pi b).'''
    f = r['faixas']
    return f['y']/av['asa']['yt'], f['ccl']/(4*av['Sref']*r['CL']/(np.pi*av['Bref']))


class Avaliador:
    def __init__(self, av, modelo='completo', etas=ETAS, it=0.0,
                 com_estol=True, nome='ot'):
        self.av, self.modelo, self.etas = av, modelo, list(etas)
        self.it, self.com_estol = it, com_estol
        self.arquivo = f'resultados/_tmp/{nome}.avl'
        self.cache, self.n_rodadas, self.ultimo = {}, 0, None
        self.historico = []

    def torcao(self, x):
        return np.concatenate([[0.0], np.asarray(x, float)])

    def avalia(self, x):
        chave = tuple(np.round(np.asarray(x, float), 8))
        if chave in self.cache:
            self.ultimo = self.cache[chave]
            return self.ultimo
        limpa = self.modelo == 'asa_limpa'
        escreve_avl(self.av, self.arquivo, cg='aft', torcao=self.torcao(x),
                    etas=self.etas, so_asa=limpa, diedro=not limpa)
        r = caso(self.arquivo, self.av['M'], cl=self.av['CL'], it=self.it,
                 trim=not limpa, empenagem=not limpa, faixas=True)
        self.n_rodadas += 1
        res = {'CDff': r['CDff'], 'e': r['e'], 'r': r}
        if self.com_estol:
            m = modelo_estol(self.arquivo, self.it, trim=True,
                             semi=self.av['asa']['yt'], empenagem=not limpa)
            self.n_rodadas += 2
            dentro = m['eta'] <= ETA_ESTOL
            res['g'] = m['alfa_i'][~dentro] - m['alfa_i'][dentro].min() - DESEMPATE
            res['estol'] = m
        self.cache[chave] = res
        self.ultimo = res
        self.historico.append((np.array(x, float), 1e4*res['CDff']))
        return res

    def objetivo(self, x):
        return 1e4*self.avalia(x)['CDff']          # [count]

    def restricoes(self, x):
        return self.avalia(x)['g']


def otimiza(avaliador, x0, maxiter=100):
    cons = ([{'type': 'ineq', 'fun': avaliador.restricoes}]
            if avaliador.com_estol else [])
    r = minimize(avaliador.objetivo, np.asarray(x0, float), method='SLSQP',
                 bounds=[LIMITES]*len(x0), constraints=cons,
                 options={'ftol': 1e-3, 'eps': EPS_FD, 'maxiter': maxiter})
    return r
```

- [ ] **Step 4: Rodar e ver passar**

Run: `python -m pytest avl/tests/test_otimizacao.py -q`
Expected: `2 passed`

- [ ] **Step 5: Criar `avl/verifica_eliptica.py`**

```python
'''
Etapa 0b -- verificacao do otimizador: asa limpa -> distribuicao eliptica.

So a asa (planta, enflechamento e perfis do projeto), plana (diedro zero),
sem winglet, fuselagem, naceles ou empenagens. O mesmo Avaliador e o mesmo
SLSQP da otimizacao de projeto minimizam o CDff no CL de projeto. Pela
teoria (Prandtl/Munk), o otimo e a carga eliptica, com e -> 1.

Criterios: e >= 0,98; desvio RMS da carga contra a elipse <= 0,02 (a carga
e normalizada pela elipse na raiz); duas partidas chegam ao mesmo CDff com
diferenca de ate 0,1 count. Se o e falhar, repete com o dobro de estacoes
para separar parametrizacao de otimizador.

Rodar de dentro de avl/:   python verifica_eliptica.py
'''

import json
import os

import numpy as np

from estilo import plt, AZUL, LARANJA, CINZA, TINTA2
from aeronave import aeronave, LC_H_BASE
from gera_avl import ETAS, AQUI
from otimizacao_torcao import (Avaliador, otimiza, carga_normalizada,
                               elipse_normalizada)

E_MIN, RMS_MAX, DIF_MAX = 0.98, 0.02, 0.1


def roda(av, etas):
    etas = np.asarray(etas)
    partidas = {'sem torção': np.zeros(len(etas) - 1),
                'washout linear −4°': -4.0*etas[1:]}
    saida = {}
    for rotulo, x0 in partidas.items():
        a = Avaliador(av, modelo='asa_limpa', etas=etas, com_estol=False,
                      nome='eliptica')
        base = a.avalia(np.zeros(len(etas) - 1))
        r = otimiza(a, x0)
        fim = a.avalia(r.x)
        eta, carga = carga_normalizada(av, fim['r'])
        rms = float(np.sqrt(np.mean((carga - elipse_normalizada(eta))**2)))
        saida[rotulo] = {'x': r.x.tolist(), 'CDff': 1e4*fim['CDff'], 'e': fim['e'],
                         'e_base': base['e'], 'CDff_base': 1e4*base['CDff'],
                         'rms': rms, 'n_rodadas': a.n_rodadas,
                         'eta': eta.tolist(), 'carga': carga.tolist(),
                         'carga_base': carga_normalizada(av, base['r'])[1].tolist()}
        print(f'  partida {rotulo:20s}: e {base["e"]:.4f} -> {fim["e"]:.4f}, '
              f'CDff {1e4*base["CDff"]:.2f} -> {1e4*fim["CDff"]:.2f} count, '
              f'RMS {rms:.4f}, {a.n_rodadas} rodadas')
    return saida


def julga(s):
    v = list(s.values())
    crit = {'e >= 0,98': min(x['e'] for x in v) >= E_MIN,
            'RMS <= 0,02': max(x['rms'] for x in v) <= RMS_MAX,
            'partidas concordam': abs(v[0]['CDff'] - v[1]['CDff']) <= DIF_MAX}
    for k, ok in crit.items():
        print(f'    {k:22s} {"ok" if ok else "FALHOU"}')
    return all(crit.values()), crit


def figura(s, nome):
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    v = s['sem torção']
    eta = np.array(v['eta'])
    ax.plot(eta, elipse_normalizada(eta), color=TINTA2, lw=1.2, ls='--', label='elíptica')
    ax.plot(eta, v['carga_base'], color=CINZA, label=f'sem torção (e = {v["e_base"]:.3f})')
    ax.plot(eta, v['carga'], color=AZUL, label=f'otimizada (e = {v["e"]:.3f})')
    ax.set_xlabel('η = 2y/b')
    ax.set_ylabel('c·cl / (4 Sref CL / π b)')
    ax.set_title('Asa limpa: a otimização de torção recupera a carga elíptica')
    ax.legend()
    fig.savefig(os.path.join(AQUI, 'resultados', nome))


def main():
    av = aeronave(LC_H_BASE)
    print('Estacoes do projeto:')
    s = roda(av, ETAS)
    ok, crit = julga(s)
    res = {'etas': ETAS, 'resultado': s, 'criterios': crit, 'aprovado': ok}
    figura(s, 'verifica_eliptica.png')
    if min(x['e'] for x in s.values()) < E_MIN:
        etas2 = sorted(set(ETAS) | set(np.round(0.5*(np.array(ETAS[1:]) + np.array(ETAS[:-1])), 4)))
        print('e abaixo de 0,98: repetindo com o dobro de estacoes')
        s2 = roda(av, etas2)
        ok2, crit2 = julga(s2)
        res['refinado'] = {'etas': etas2, 'resultado': s2, 'criterios': crit2,
                           'aprovado': ok2}
        figura(s2, 'verifica_eliptica_refinada.png')
    with open(os.path.join(AQUI, 'resultados', 'verifica_eliptica.json'), 'w') as f:
        json.dump(res, f, indent=2)


if __name__ == '__main__':
    main()
```

- [ ] **Step 6: Rodar a verificação**

Run: `cd avl && python verifica_eliptica.py`
Expected: as duas partidas saem com e ≥ 0,98 e RMS ≤ 0,02, e os três
critérios dão `ok`. Abra `resultados/verifica_eliptica.png`: a curva azul
tem de estar sobre a tracejada. **Se algum critério falhar, pare e investigue
antes da Task 9.** Primeiro, confira se `EPS_FD` é grande demais (0,2° →
teste 0,1°). Depois, confira se o SLSQP parou por `maxiter`: imprima
`r.message` em `roda`.

- [ ] **Step 7: Commit**

```bash
git add avl/otimizacao_torcao.py avl/verifica_eliptica.py avl/tests/test_otimizacao.py avl/resultados/verifica_eliptica.json avl/resultados/verifica_eliptica*.png
git commit -m "Lab 04: otimizador de torcao verificado na asa limpa contra a carga eliptica"
```

---

### Task 9: `varredura_eh.py` (etapa 1)

**Files:**
- Create: `avl/varredura_eh.py`

- [ ] **Step 1: Criar `avl/varredura_eh.py`**

```python
'''
Etapa 1 -- varredura da posicao longitudinal da EH (Lc_h), asa sem torcao.

Para cada Lc_h: designTool (W0, S_h, CL de projeto, CGs, dε/dα) -> .avl
fwd/aft -> it que zera o profundor no cruzeiro, CD compensado, margem
estatica -> estol (secao critica, compensado) -> sombreamento e leme.
Escolhe o MAIOR Lc_h viavel.

Rodar de dentro de avl/:   python varredura_eh.py
'''

import csv
import json
import os

import numpy as np

from estilo import plt, AZUL, LARANJA, CINZA, TINTA2
from aeronave import aeronave, lc_h_maximo, FOLGA_FUSELAGEM
from analises import it_para_de_zero
from estol import estol
from gera_avl import escreve_avl, AQUI
from sombra import alfa_saida, fracao_leme_encoberto

LC_MIN, PASSO = 4.0, 0.1
FOLGA_SOMBRA = 2.0      # [graus] a EH sai da sombra ao menos isso antes do estol
DE_MAX = 20.0           # [graus] profundor na compensacao no CLmax, CG dianteiro
MS_MIN = 0.05           # margem estatica minima, CG traseiro
SAIDA = os.path.join(AQUI, 'resultados')


def ponto(lc):
    av = aeronave(lc)
    semi = av['asa']['yt']
    r = {'Lc_h': lc, 'S_h': av['EH']['S'], 'W0_kgf': av['W0_kgf'], 'CL': av['CL'],
         'folga_fus': av['fuselagem']['L'] - (av['EH']['xr'] + av['EH']['cr'])}
    for cg in ('fwd', 'aft'):
        arq = escreve_avl(av, f'resultados/varredura/{cg}_{lc:.3f}.avl', cg=cg)
        it, c = it_para_de_zero(arq, av['M'], av['CL'])
        e = estol(arq, it=it, trim=True, semi=semi)
        xcg = av[f'xcg_{cg}']
        r.update({f'it_{cg}': it, f'CD_{cg}': 1e4*c['CD'], f'CDff_{cg}': 1e4*c['CDff'],
                  f'e_{cg}': c['e'], f'MS_{cg}': (c['xnp'] - xcg)/av['Cref'],
                  f'alfa_estol_{cg}': e['alfa'], f'CLmax_{cg}': e['CL'],
                  f'de_estol_{cg}': e['de'], f'eta_estol_{cg}': e['eta_crit']})
    r['alfa_estol'] = min(r['alfa_estol_fwd'], r['alfa_estol_aft'])
    r['alfa_saida'] = alfa_saida(av)
    r['folga_sombra'] = r['alfa_estol'] - r['alfa_saida']
    r['leme_encoberto'] = fracao_leme_encoberto(av)
    r['ok_fus'] = r['folga_fus'] >= FOLGA_FUSELAGEM - 1e-6
    r['ok_sombra'] = r['folga_sombra'] >= FOLGA_SOMBRA
    r['ok_profundor'] = abs(r['de_estol_fwd']) <= DE_MAX
    r['ok_ms'] = r['MS_aft'] >= MS_MIN
    # a sombra so e exigida depois da torcao (etapa 2): aqui e informativa
    r['viavel'] = r['ok_fus'] and r['ok_profundor'] and r['ok_ms']
    return r


def figura(linhas, escolhido):
    lc = np.array([l['Lc_h'] for l in linhas])
    paineis = [
        ('W0 [kgf]', [('W0_kgf', AZUL, None)], None),
        ('CD compensado [count]', [('CD_fwd', LARANJA, 'CG dianteiro'),
                                   ('CD_aft', AZUL, 'CG traseiro')], None),
        ('δe no CLmax, CG dianteiro [°]', [('de_estol_fwd', AZUL, None)], -DE_MAX),
        ('margem estática, CG traseiro', [('MS_aft', AZUL, None)], MS_MIN),
        ('α_estol − α_saída da sombra [°]', [('folga_sombra', AZUL, None)], FOLGA_SOMBRA),
        ('fração do leme encoberta', [('leme_encoberto', AZUL, None)], None),
    ]
    fig, eixos = plt.subplots(2, 3, figsize=(10.5, 5.6), sharex=True)
    for ax, (titulo, series, limite) in zip(eixos.flat, paineis):
        for chave, cor, rot in series:
            ax.plot(lc, [l[chave] for l in linhas], '-o', color=cor, label=rot)
        if limite is not None:
            ax.axhline(limite, color=TINTA2, lw=0.8, ls='--')
        for l in linhas:
            if not l['viavel']:
                ax.axvspan(l['Lc_h'] - PASSO/2, l['Lc_h'] + PASSO/2, color='#f1f0ea', lw=0, zorder=0)
        if escolhido is not None:
            ax.axvline(escolhido, color=TINTA2, lw=1.0)
        ax.set_title(titulo, fontsize=9)
        if any(rot for _, _, rot in series):
            ax.legend(fontsize=8)
    for ax in eixos[1]:
        ax.set_xlabel('Lc_h (braço da EH / CMA)')
    fig.suptitle('Varredura da posição da EH (fundo claro: inviável; linha vertical: escolhida; tracejado: limite)')
    fig.savefig(os.path.join(SAIDA, 'varredura_eh.png'))


def main():
    lc_max = lc_h_maximo()
    grade = list(np.round(np.arange(LC_MIN, lc_max - 1e-9, PASSO), 3)) + [round(lc_max, 4)]
    print(f'Lc_h de {LC_MIN} a {lc_max:.4f} (limite da fuselagem, folga {FOLGA_FUSELAGEM} m)')
    linhas = []
    for lc in grade:
        l = ponto(float(lc))
        linhas.append(l)
        print(f'  Lc_h {lc:6.3f}  W0 {l["W0_kgf"]:9.0f}  CD_aft {l["CD_aft"]:7.2f}'
              f'  de_estol {l["de_estol_fwd"]:6.2f}  MS_aft {l["MS_aft"]:6.3f}'
              f'  folga_sombra {l["folga_sombra"]:5.2f}  leme {l["leme_encoberto"]:.2f}'
              f'  {"VIAVEL" if l["viavel"] else "inviavel"}')

    viaveis = [l for l in linhas if l['viavel']]
    escolhido = max(l['Lc_h'] for l in viaveis) if viaveis else None
    w0 = [l['W0_kgf'] for l in linhas]
    cd = [l['CD_aft'] for l in linhas]
    monotono = bool(np.all(np.diff(w0) < 0) and np.all(np.diff(cd) < 0))
    if viaveis and not monotono:
        escolhido = min(viaveis, key=lambda l: l['CD_aft'])['Lc_h']
        print('W0 ou CD nao caem monotonamente: escolha pelo menor CD viavel')
    print(f'\nLc_h escolhido: {escolhido}   (W0 e CD monotonos: {monotono})')

    os.makedirs(SAIDA, exist_ok=True)
    with open(os.path.join(SAIDA, 'varredura_eh.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)
    with open(os.path.join(SAIDA, 'lc_h_escolhido.json'), 'w') as f:
        json.dump({'Lc_h': escolhido, 'monotono': monotono, 'lc_max': lc_max,
                   'restricoes': {'folga_fus_m': FOLGA_FUSELAGEM,
                                  'folga_sombra_graus': FOLGA_SOMBRA,
                                  'de_max_graus': DE_MAX, 'ms_min': MS_MIN}}, f, indent=2)
    figura(linhas, escolhido)


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Rodar a varredura**

Run: `cd avl && python varredura_eh.py`
Expected: uma linha por `Lc_h` (≈ 10 pontos) e `Lc_h escolhido: X`. Pela
especificação, W0 e CD devem cair com `Lc_h`. Abra
`resultados/varredura_eh.png` e confira qual restrição fica ativa.
**Se nenhum ponto for viável, pare e leve a tabela ao usuário:** os limites
de 2°, 20° e 5% são premissas a revisar e não devem ser afrouxados por
conta própria.

- [ ] **Step 3: Conferir a malha no `Lc_h` escolhido**

Run: `cd avl && python convergencia_malha.py --lc <Lc_h escolhido>`
Expected: `-> ok`

- [ ] **Step 4: Commit**

```bash
git add avl/varredura_eh.py avl/resultados/varredura_eh.csv avl/resultados/varredura_eh.png avl/resultados/lc_h_escolhido.json
git commit -m "Lab 04: varredura da posicao da EH com restricoes de fuselagem, sombra, profundor e margem"
```

`avl/resultados/varredura/` (os `.avl` de cada ponto) fica fora do commit;
acrescente `avl/resultados/varredura/` ao `.gitignore` neste mesmo commit.

---

### Task 10: otimização de torção na posição escolhida (etapa 2)

**Files:**
- Modify: `avl/otimizacao_torcao.py` (acrescentar `main` ao final)

- [ ] **Step 1: Acrescentar o `main` ao fim de `avl/otimizacao_torcao.py`**

```python
# ======================================================================
# ETAPA 2 -- execucao

def _projeto(lc):
    from aeronave import aeronave
    from analises import it_para_de_zero
    av = aeronave(lc)
    base = escreve_avl(av, f'resultados/_tmp/base_{lc:.3f}.avl', cg='aft')
    it, _ = it_para_de_zero(base, av['M'], av['CL'])
    a = Avaliador(av, modelo='completo', it=it, com_estol=True, nome=f'ot_{lc:.3f}')
    x0 = np.zeros(len(ETAS) - 1)
    f0 = a.objetivo(x0)
    r = otimiza(a, x0)
    fim = a.avalia(r.x)
    viavel = bool(np.all(fim['g'] >= -1e-3))
    print(f'  Lc_h {lc:.3f}: CDff {f0:.2f} -> {1e4*fim["CDff"]:.2f} count, '
          f'e {fim["e"]:.4f}, viavel {viavel}, {a.n_rodadas} rodadas, {r.message}')
    return av, it, a, r, fim, viavel


def main():
    from estilo import plt, AZUL, LARANJA, VERDE, CINZA, TINTA2
    from estol import estol

    with open(os.path.join(AQUI, 'resultados', 'lc_h_escolhido.json')) as f:
        lc = json.load(f)['Lc_h']

    print('Otimizacao de torcao na posicao escolhida e em dois vizinhos')
    casos = {}
    for l in (lc, round(lc - 0.2, 4), round(lc - 0.4, 4)):
        casos[l] = _projeto(l)

    av, it, a, r, fim, viavel = casos[lc]
    tw = a.torcao(r.x)

    # conferencia direta do estol no projeto final
    for cg in ('fwd', 'aft'):
        arq = escreve_avl(av, f'{cg}.avl', cg=cg, torcao=tw)
    e = estol('aft.avl', it=it, trim=True, semi=av['asa']['yt'])
    print(f'  estol conferido direto: comeca em eta = {e["eta_crit"]:.3f} '
          f'(exigido <= {ETA_ESTOL}), alfa = {e["alfa"]:.2f}, CLmax = {e["CL"]:.3f}')
    # sombreamento com o alfa de estol da asa torcida (menor entre os CGs)
    from sombra import alfa_saida
    it_f, _ = it_para_de_zero('fwd.avl', av['M'], av['CL'])
    e_f = estol('fwd.avl', it=it_f, trim=True, semi=av['asa']['yt'])
    a_estol = min(e['alfa'], e_f['alfa'])
    folga_sombra = a_estol - alfa_saida(av)
    print(f'  sombra: alfa_estol {a_estol:.2f} - alfa_saida {alfa_saida(av):.2f} '
          f'= {folga_sombra:.2f} graus (exigido >= 2) -> '
          f'{"ok" if folga_sombra >= 2.0 else "FALHOU: levar a equipe"}')

    dif = {l: float(np.max(np.abs(c[2].torcao(c[3].x) - tw))) for l, c in casos.items()}
    print('  maior diferenca de torcao contra a posicao escolhida:',
          {k: round(v, 2) for k, v in dif.items()})

    res = {'Lc_h': lc, 'it_base': it, 'etas': ETAS, 'torcao': tw.tolist(),
           'CDff_count': 1e4*fim['CDff'], 'e': fim['e'], 'viavel': viavel,
           'estol_direto': {'eta_crit': e['eta_crit'], 'alfa': e['alfa'],
                            'CLmax': e['CL'], 'ok': e['eta_crit'] <= ETA_ESTOL},
           'vizinhos': {str(l): {'torcao': c[2].torcao(c[3].x).tolist(),
                                 'CDff_count': 1e4*c[4]['CDff'],
                                 'dif_max_torcao': dif[l]} for l, c in casos.items()}}
    with open(os.path.join(AQUI, 'resultados', 'torcao_otimizada.json'), 'w') as f:
        json.dump(res, f, indent=2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.6))
    for (l, c), cor in zip(casos.items(), (AZUL, LARANJA, VERDE)):
        ax1.plot(ETAS, c[2].torcao(c[3].x), '-o', color=cor, label=f'Lc_h = {l}')
    ax1.set_xlabel('η')
    ax1.set_ylabel('torção [°]')
    ax1.legend()
    ax1.set_title('Torção ótima')
    ax2.plot(e['eta'], e['cl_norm'], color=AZUL, label='cl_norm no estol')
    ax2.plot(e['eta'], e['clmax'], color=TINTA2, ls='--', lw=1.2, label='clmax do perfil')
    ax2.axvline(ETA_ESTOL, color=CINZA, lw=0.8)
    ax2.set_xlabel('η')
    ax2.legend()
    ax2.set_title(f'Estol começa em η = {e["eta_crit"]:.2f}')
    fig.savefig(os.path.join(AQUI, 'resultados', 'otimizacao_torcao.png'))


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Rodar**

Run: `cd avl && python otimizacao_torcao.py`
Expected, nesta ordem:
- três linhas `Lc_h ...: CDff A -> B count`, com B < A e `viavel True`;
- `estol conferido direto: comeca em eta <= 0.56`;
- as diferenças de torção entre a posição escolhida e os vizinhos.

Se a diferença passar de ~1°, a premissa "varredura antes, torção depois"
não vale. Nesse caso, avise o usuário antes de seguir: a varredura teria de
ser refeita com a torção ótima. O `fwd.avl` e o `aft.avl` da raiz de `avl/`
ficam regerados com a torção ótima.

- [ ] **Step 3: Rodar a suíte inteira**

Run: `python -m pytest avl/tests -q`
Expected: todos passam. A referência da regressão é o arquivo congelado,
não o `fwd.avl` regerado.

- [ ] **Step 4: Commit**

```bash
git add avl/otimizacao_torcao.py avl/fwd.avl avl/aft.avl avl/resultados/torcao_otimizada.json avl/resultados/otimizacao_torcao.png
git commit -m "Lab 04: torcao otimizada na posicao escolhida da EH, com conferencia de independencia"
```

---

### Task 11: registrar o `Lc_h` no designTool e escrever o README

**Files:**
- Modify: `designTool/standard_airplane.py:171`
- Create: `avl/README.md`

- [ ] **Step 1: Atualizar o `Lc_h`**

Na linha 171, trocar `'Lc_h' : 4.6,` pelo valor de
`avl/resultados/lc_h_escolhido.json`, mantendo o histórico no comentário:

```python
                  'Lc_h' : <VALOR>, # Non-dimensional lever of the horizontal tail (lever/wing_mac) (v2: 4.0; v3: 4.6; Lab 04: varredura em avl/varredura_eh.py)
```

(`<VALOR>` é o número do JSON, com 4 casas.)

- [ ] **Step 2: Conferir que o designTool roda e que o W0 caiu**

Run: `python -c "from designTool.standard_airplane import standard_airplane; from designTool.analyze import analyze; a=standard_airplane('my_airplane'); analyze(a); print(a['thrust_matching']['W0']/9.81)"`
Expected: W0 [kgf] igual ao da linha escolhida em `varredura_eh.csv`.

- [ ] **Step 3: Escrever `avl/README.md`**

Conteúdo:
- **Ordem de execução**, de dentro de `avl/`: `python convergencia_malha.py`,
  `python verifica_eliptica.py`, `python varredura_eh.py`,
  `python otimizacao_torcao.py`. Testes: `python -m pytest tests -q`.
- **Deslocamento em Z:** existe por exigência do professor e o sombreamento
  usa as alturas reais.
- **Numeração dos controles:** `d1` aileron, `d2` profundor, `d3` leme; a
  variável de projeto `1` é o `it`. Ela difere da do roteiro (`d4`, `2`).
- **Tabela de resultados**, preenchida com os números dos JSONs:
  - malha adotada;
  - verificação elíptica (e, RMS);
  - `Lc_h` escolhido e a restrição ativa;
  - W0 antes e depois;
  - CDff antes e depois da torção;
  - η e α do estol.
- **Premissas:**
  - limites de 0,5 m, 2°, 20° e 5%;
  - dε/dα pelo Roskam;
  - clmax dos perfis do Lab 03;
  - peso de projeto com a fração de combustível do Lab 03.

- [ ] **Step 4: Commit**

```bash
git add designTool/standard_airplane.py avl/README.md
git commit -m "Lab 04: Lc_h escolhido no designTool e README da pasta avl"
```

- [ ] **Step 5: Perguntar ao usuário se pode fazer o push**

A receita está na memória do projeto: alternar para a conta
`JoseJunior4144`, `git push origin main` e depois voltar para a conta
anterior.

---

## Fora deste plano

Os itens 1–8 do roteiro e as derivadas de estabilidade (etapa 3 da
especificação) terão plano próprio. Eles reaproveitam `avl_run.caso`,
`analises.it_para_de_zero` e `estol.estol` com os `.avl` finais.
