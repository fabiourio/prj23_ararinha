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

from estilo import plt, AZUL, CINZA, TINTA2
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
                         'rms': rms, 'n_rodadas': a.n_rodadas, 'mensagem': r.message,
                         'eta': eta.tolist(), 'carga': carga.tolist(),
                         'carga_base': carga_normalizada(av, base['r'])[1].tolist()}
        print(f'  partida {rotulo:20s}: e {base["e"]:.4f} -> {fim["e"]:.4f}, '
              f'CDff {1e4*base["CDff"]:.2f} -> {1e4*fim["CDff"]:.2f} count, '
              f'RMS {rms:.4f}, {a.n_rodadas} rodadas, {r.message}')
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
        json.dump(res, f, indent=2, default=str)


if __name__ == '__main__':
    main()
