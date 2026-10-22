import os
import sys

AQUI = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RAIZ = os.path.abspath(os.path.join(AQUI, '..'))
for p in (AQUI, RAIZ, os.path.join(RAIZ, 'avl')):
    if p not in sys.path:
        sys.path.insert(0, p)
