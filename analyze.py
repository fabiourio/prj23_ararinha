'''Demo rapida: analyze() no Fokker 100 com log e vista 3D.'''

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze

airplane = standard_airplane('fokker100')
analyze(airplane, print_log=True, plot=True)
