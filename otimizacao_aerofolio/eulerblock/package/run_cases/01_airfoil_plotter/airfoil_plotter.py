'''
Run this script to see a graphical interface to generate airfoils using
the CST parameters. This interface also lets the use run Euler analyses
and export the airfoil coordinates in Xfoil format.

Written by:
Ney Rafael Secco
Instituto Tecnologico de Aeronautica
Sao Jose dos Campos - Brazil
ney@ita.br
Sept 2022
'''

# Modify path to include the Euler solver
import sys
with open('../eulerblock_path.txt') as f:
    exec(f.read())
sys.path.append(eulerblock_path)
from eulerblock import euler_mod as eb
import airfoil_mod as am

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import ctypes

# Initial set of parameters
cref=1.0
mach=0.75
iter=50000
CFL=0.2
res_tol=1e-6
reinitialize=0
Nchord=31

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.29, bottom=0.35)

# Initial airfoil
Al10 = -0.2
Al20 = -0.2
Al30 = -0.2
Al40 = -0.2
Au10 =  0.2
Au20 =  0.2
Au30 =  0.2
Au40 =  0.2
AoA0 =  0.0

# Chord positions
xa = (1-np.cos(np.linspace(0, 1, Nchord)*np.pi))/2

Al0 = [Al10, Al20, Al30, Al40]
Au0 = [Au10, Au20, Au30, Au40]
airfoil = am.cstfoil(Au0, Al0, xa, N1=0.5, N2=1.0, tte=0.0, plot=False)
xf = airfoil['x_coord']
yf = airfoil['y_coord']
maxt = airfoil['max_thickness']
xmaxt = airfoil['x_max_thickness']
mint = airfoil['min_thickness']
maxc = airfoil['max_camber']
xmaxc = airfoil['x_max_camber']
ax.set_title(r'$t/c_{\max}=%.3f \; x_{t/c\max}=%.3f \; h/c_{\max}=%.3f \; x_{h/c\max}=%.3f$'%(maxt, xmaxt, maxc, xmaxc),fontsize=8)

curve, = plt.plot(xf, yf, lw=2, color='red')
plt.axis('equal')

axcolor = 'lightgoldenrodyellow'
axAl1 = plt.axes([0.25, 0.27, 0.65, 0.01], facecolor=axcolor)
axAl2 = plt.axes([0.25, 0.25, 0.65, 0.01], facecolor=axcolor)
axAl3 = plt.axes([0.25, 0.23, 0.65, 0.01], facecolor=axcolor)
axAl4 = plt.axes([0.25, 0.21, 0.65, 0.01], facecolor=axcolor)
axAu1 = plt.axes([0.25, 0.19, 0.65, 0.01], facecolor=axcolor)
axAu2 = plt.axes([0.25, 0.17, 0.65, 0.01], facecolor=axcolor)
axAu3 = plt.axes([0.25, 0.15, 0.65, 0.01], facecolor=axcolor)
axAu4 = plt.axes([0.25, 0.13, 0.65, 0.01], facecolor=axcolor)
axAoA = plt.axes([0.25, 0.11, 0.65, 0.01], facecolor=axcolor)

sAl1 = Slider(axAl1, 'Al1', -0.4, -0.05, valinit=Al10)
sAl2 = Slider(axAl2, 'Al2', -0.4, -0.05, valinit=Al20)
sAl3 = Slider(axAl3, 'Al3', -0.4, -0.05, valinit=Al30)
sAl4 = Slider(axAl4, 'Al4', -0.4, -0.05, valinit=Al40)
sAu1 = Slider(axAu1, 'Au1', 0.05, 0.4, valinit=Au10)
sAu2 = Slider(axAu2, 'Au2', 0.05, 0.4, valinit=Au20)
sAu3 = Slider(axAu3, 'Au3', 0.05, 0.4, valinit=Au30)
sAu4 = Slider(axAu4, 'Au4', 0.05, 0.4, valinit=Au40)
sAoA = Slider(axAoA, 'AoA', -5.0, 5.0, valinit=AoA0)


def update(val):
    Al1 = sAl1.val
    Al2 = sAl2.val
    Al3 = sAl3.val
    Al4 = sAl4.val
    Au1 = sAu1.val
    Au2 = sAu2.val
    Au3 = sAu3.val
    Au4 = sAu4.val
    AoA = sAoA.val*np.pi/180.0

    Al = [Al1, Al2, Al3, Al4]
    Au = [Au1, Au2, Au3, Au4]
    airfoil = am.cstfoil(Au, Al, xa, N1=0.5, N2=1.0, tte=0.0, plot=False)
    xf = airfoil['x_coord']
    yf = airfoil['y_coord']
    maxt = airfoil['max_thickness']
    xmaxt = airfoil['x_max_thickness']
    mint = airfoil['min_thickness']
    maxc = airfoil['max_camber']
    xmaxc = airfoil['x_max_camber']
    ax.set_title(r'$t/c_{\max}=%.3f \; x_{t/c\max}=%.3f \; h/c_{\max}=%.3f \; x_{h/c\max}=%.3f$'%(maxt, xmaxt, maxc, xmaxc),fontsize=8)

    curve.set_xdata(xf)
    curve.set_ydata(yf)
    fig.canvas.draw_idle()

for ss in [sAl1, sAl2, sAl3, sAl4, sAu1, sAu2, sAu3, sAu4, sAoA]:
    ss.on_changed(update)

#------------------------
# Run button

runax = plt.axes([0.65, 0.025, 0.1, 0.04])
button_run = Button(runax, 'Run', color=axcolor, hovercolor='0.975')

def run(event):
    Al1 = sAl1.val
    Al2 = sAl2.val
    Al3 = sAl3.val
    Al4 = sAl4.val
    Au1 = sAu1.val
    Au2 = sAu2.val
    Au3 = sAu3.val
    Au4 = sAu4.val
    AoA = sAoA.val*np.pi/180.0

    Al = [Al1, Al2, Al3, Al4]
    Au = [Au1, Au2, Au3, Au4]
    results = eb.run_cst(Al, Au, Nchord, AoA, mach,
                         gamma=1.4, order=2,
                         iter=iter, CFL=CFL, use_local_dt=1,
                         res_tol=res_tol, reinitialize=reinitialize, plot=True)

button_run.on_clicked(run)

#------------------------

# Reset button
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button_reset = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    for ss in [sAl1, sAl2, sAl3, sAl4, sAu1, sAu2, sAu3, sAu4, sAoA]:
        ss.reset()
button_reset.on_clicked(reset)

#------------------------

# Export button
exportax = plt.axes([0.5, 0.025, 0.1, 0.04])
button_export = Button(exportax, 'Export', color=axcolor, hovercolor='0.975')

def export(event):

    Al1 = sAl1.val
    Al2 = sAl2.val
    Al3 = sAl3.val
    Al4 = sAl4.val
    Au1 = sAu1.val
    Au2 = sAu2.val
    Au3 = sAu3.val
    Au4 = sAu4.val
    AoA = sAoA.val*np.pi/180.0

    Al = [Al1, Al2, Al3, Al4]
    Au = [Au1, Au2, Au3, Au4]

    xa = (1-np.cos(np.linspace(0, 1, 81)*np.pi))/2
    airfoil = am.cstfoil(Au, Al, xa, N1=0.5, N2=1.0, tte=0.0, plot=False)

    am.export_airfoil(airfoil)

    ctypes.windll.user32.MessageBoxW(0, "Coordinates exported to airfoil.dat", "Message", 0)

button_export.on_clicked(export)

#------------------------
# Mach
def submit_mach(expression):
    global mach
    mach = eval(expression)

axbox_mach = fig.add_axes([0.075, 0.7, 0.1, 0.075])
text_box_mach = TextBox(axbox_mach, "Mach")
text_box_mach.on_submit(submit_mach)
text_box_mach.set_val("%g"%mach)  # Trigger `submit` with the initial string.

#------------------------
# Iter
def submit_iter(expression):
    global iter
    iter = eval(expression)

axbox_iter = fig.add_axes([0.075, 0.6, 0.1, 0.075])
text_box_iter = TextBox(axbox_iter, "Iter")
text_box_iter.on_submit(submit_iter)
text_box_iter.set_val("%g"%iter)  # Trigger `submit` with the initial string.

#------------------------
# CFL
def submit_cfl(expression):
    global CFL
    CFL = eval(expression)

axbox_cfl = fig.add_axes([0.075, 0.5, 0.1, 0.075])
text_box_cfl = TextBox(axbox_cfl, "CFL")
text_box_cfl.on_submit(submit_cfl)
text_box_cfl.set_val("%g"%CFL)  # Trigger `submit` with the initial string.

#------------------------
# res_tol
def submit_res(expression):
    global res_tol
    res_tol = eval(expression)

axbox_res = fig.add_axes([0.075, 0.4, 0.1, 0.075])
text_box_res = TextBox(axbox_res, "res tol")
text_box_res.on_submit(submit_res)
text_box_res.set_val("%g"%res_tol)  # Trigger `submit` with the initial string.

#------------------------
# Reinitialize box

rax = plt.axes([0.025, 0.10, 0.15, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('yes', 'no'), active=1)

plt.title('Reinitialize?')
def colorfunc(status):
    global reinitialize
    if status == 'yes':
        reinitialize = 1
    elif status == 'no':
        reinitialize = 0
    print(reinitialize)

radio.on_clicked(colorfunc)

#------------------------
# Title
axtitle = plt.axes([0.15, 0.95, 0.75, 0.05], facecolor=None)
axtitle.text(0.1,0,'2D Euler Solver', fontsize=14)
axtitle.text(0.6,-0.25,'Cap. Ney Sêcco (ney@ita.br)\nInstituto Tecnológico de Aeronáutica - Brazil\nPRJ-23 Course / 2023', fontsize=6)
axtitle.spines['top'].set_visible(False)
axtitle.spines['right'].set_visible(False)
axtitle.spines['bottom'].set_visible(False)
axtitle.spines['left'].set_visible(False)
axtitle.get_xaxis().set_ticks([])
axtitle.get_yaxis().set_ticks([])


plt.show()
