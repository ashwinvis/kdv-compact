import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from ipywidgets import interact, widgets

import finite_diff


N = 256
xmax = 100.
tmax = 30.
x = np.linspace(0, xmax, N)
x0 = xmax / 2
dx = xmax / (N - 1)
dx3 = dx**3
dict_convection = dict(
    CD2=1, UD3=2, CD4=3, CD4_non_conservative=4, OUCS3=5)
dict_dispersion = dict(
    CD2=1, CD4=2, O3CS=3)
dict_convection[''] = 5
dict_dispersion[''] = 3


print('Available convection schemes are ', dict_convection.keys())
convection = dict_convection[input('Choose: ')]
print('Available dispersion schemes are ', dict_dispersion.keys())
dispersion = dict_dispersion[input('Choose: ')]



t_last_print = 0

def kdv(t, u):
    global t_last_print
    uux = finite_diff.convection(convection, u, dx)
    uxxx = finite_diff.dispersion(dispersion, u, dx3)
    rhs = -6 * uux - uxxx

    # Output params
    dt_print = 2
    if t > (t_last_print + dt_print):
        t_last_print = t
        print(f'Time t={t:3.3f}, u_max={u.max():3.3f}')

    return rhs


def initialize(x, k=1, nb_solitons=1):
    if nb_solitons == 1:
        amp = 2.
    elif nb_solitons == 2:
        amp = 6

    u0 = amp * (1 - np.tanh(k * (x - x0))**2)
    return u0

# Solve
u0 = initialize(x)
sol = solve_ivp(kdv, (0, tmax), u0)

# Plot
def surfplot():
    T, X = np.meshgrid(sol.t, x)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T, sol.y, rstride=10, cstride=10)
    plt.show()

ani_fig, ani_ax = plt.subplots()
line, = ani_ax.plot(x, u0)

def ani_update(i):
    t = sol.t[i]
    u = sol.y[:,i]
    line.set_ydata(u)
    ani_ax.set_title(f'Time = {t:3.3f}, u_max = {u.max():3.3f}')

def animate():
    interval = 0.2  # delay between frames in milliseconds
    FuncAnimation(
        ani_fig, ani_update, len(sol.t), blit=False, interval=interval, repeat=True)

        
slider = widgets.IntSlider(
    min=0, max=len(sol.t), step=1)
@interact(iteration=slider)
def animate_in_jupyter(iteration):
    ani_update(iteration)
    # ani_fig.canvas.draw()


plt.show()
