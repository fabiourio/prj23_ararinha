import os
import sys

AVL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if AVL_DIR not in sys.path:
    sys.path.insert(0, AVL_DIR)
