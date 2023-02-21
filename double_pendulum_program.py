import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from tkinter import ttk

# set parameters for simulation (masses and lengths)
m1 = 1.0
m2 = 1.0
l1 = 1.0
l2 = 1.0
g = -9.8  # like a true engineer
frame = (l1 + l2) * 1.125

# stepping information
t0 = 0  # initial time
dt = 0.025

# initial conditions with angles measured clockwise from vertical
theta1 = np.radians(45)
theta2 = np.radians(0)
omega1 = 0
omega2 = 0  # omegas are placeholders for derivs since we have second order eq's


# define functions for the double pendulum equations of motion
def f1(theta1, theta2, omega1, omega2):  # deriv of theta1
    global m1, m2, l1, l2, g
    return omega1


def f2(theta1, theta2, omega1, omega2):  # deriv of theta2
    global m1, m2, l1, l2, g
    return omega2


def f3(theta1, theta2, omega1, omega2):  # deriv of omega1
    global m1, m2, l1, l2, g

    numer = -g*(2*m1 + m2) * np.sin(theta1) - m2*g * np.sin(theta1 - 2*theta2) \
        - 2*np.sin(theta1 - theta2) * m2 * (omega2**2 * l2 + omega1**2 * l1 * np.cos(theta1 - theta2))
    denom = l1 * (2*m1 + m2 - m2*np.cos(2*(theta1 - theta2)))
    return numer / denom


def f4(theta1, theta2, omega1, omega2):  # deriv of omega2
    global m1, m2, l1, l2, g

    numer = 2*np.sin(theta1 - theta2) * (omega1**2 * l1 * (m1 + m2) + g*(m1 + m2)*np.cos(theta1) \
        + omega2**2 * l2*m2*np.cos(theta1 - theta2))
    denom = l2 * (2*m1 + m2 - m2*np.cos(2*(theta1 - theta2)))
    return numer / denom


# main loop using fourth order runge kutta
def pendulum_rk4(theta1, theta2, omega1, omega2, mode):
    global m1, m2
    m1 = mass1.get()
    m2 = mass2.get()
    plt.cla()
    t = t0
    p1 = [l1 * np.sin(theta1), l1 * np.cos(theta1)]
    p2 = [p1[0] + l2 * np.sin(theta2), p1[1] + l2 * np.cos(theta2)]

    while True:
        prev = [p1[0], p1[1], p2[0], p2[1]]

        # compute the various intermediate points for RK-4
        # first index corresponds to which function/equation were working on
        # second index corresponds to which intermediate point were working on calculation

        k11 = f1(theta1, theta2, omega1, omega2)
        k21 = f2(theta1, theta2, omega1, omega2)
        k31 = f3(theta1, theta2, omega1, omega2)
        k41 = f4(theta1, theta2, omega1, omega2)

        k12 = f1(theta1 + k11 * (dt/2), theta2 + k21 * (dt/2), omega1 + k31 * (dt/2), omega2 + k41 * (dt/2))
        k22 = f2(theta1 + k11 * (dt/2), theta2 + k21 * (dt/2), omega1 + k31 * (dt/2), omega2 + k41 * (dt/2))
        k32 = f3(theta1 + k11 * (dt/2), theta2 + k21 * (dt/2), omega1 + k31 * (dt/2), omega2 + k41 * (dt/2))
        k42 = f4(theta1 + k11 * (dt/2), theta2 + k21 * (dt/2), omega1 + k31 * (dt/2), omega2 + k41 * (dt/2))

        k13 = f1(theta1 + k12 * (dt/2), theta2 + k22 * (dt/2), omega1 + k32 * (dt/2), omega2 + k42 * (dt/2))
        k23 = f2(theta1 + k12 * (dt/2), theta2 + k22 * (dt/2), omega1 + k32 * (dt/2), omega2 + k42 * (dt/2))
        k33 = f3(theta1 + k12 * (dt/2), theta2 + k22 * (dt/2), omega1 + k32 * (dt/2), omega2 + k42 * (dt/2))
        k43 = f4(theta1 + k12 * (dt/2), theta2 + k22 * (dt/2), omega1 + k32 * (dt/2), omega2 + k42 * (dt/2))

        k14 = f1(theta1 + k13 * dt, theta2 + k23 * dt, omega1 + k33 * dt, omega2 + k43 * dt)
        k24 = f2(theta1 + k13 * dt, theta2 + k23 * dt, omega1 + k33 * dt, omega2 + k43 * dt)
        k34 = f3(theta1 + k13 * dt, theta2 + k23 * dt, omega1 + k33 * dt, omega2 + k43 * dt)
        k44 = f4(theta1 + k13 * dt, theta2 + k23 * dt, omega1 + k33 * dt, omega2 + k43 * dt)

        theta1 += (k11/6 + k12/3 + k13/3 + k14/6) * dt
        theta2 += (k21/6 + k22/3 + k23/3 + k24/6) * dt
        omega1 += (k31/6 + k32/3 + k33/3 + k34/6) * dt
        omega2 += (k41/6 + k42/3 + k43/3 + k44/6) * dt

        # update the graphics program
        if mode == 'pendulum':
            p1 = [l1 * np.sin(theta1), l1 * np.cos(theta1)]
            p2 = [p1[0] + l2 * np.sin(theta2), p1[1] + l2 * np.cos(theta2)]
            plt.xlim(-frame, frame)
            plt.ylim(-frame, frame)
            plt.plot([0, p1[0]], [0, p1[1]], color='green', linewidth=0.3)
            plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color='blue', linewidth=0.3)
            plt.pause(0.005)
            plt.cla()
        elif mode == 'path':
            p1 = [l1 * np.sin(theta1), l1 * np.cos(theta1)]
            p2 = [p1[0] + l2 * np.sin(theta2), p1[1] + l2 * np.cos(theta2)]
            plt.xlim(-frame, frame)
            plt.ylim(-frame, frame)
            plt.plot([prev[0], p1[0]], [prev[1], p1[1]], color='green', linewidth=0.3)
            plt.plot([prev[2], p2[0]], [prev[3], p2[1]], color='blue', linewidth=0.3)
            plt.pause(0.005)

        t += dt



# build the window
root = Tk()
title = Label(root, text='Double Pendulum!')

# set parameters to be determined by the user
angle1 = IntVar()
angle1.set(90)
angle2 = IntVar()
angle2.set(90)
mode = StringVar()
mode.set('pendulum')
mass1 = IntVar()
mass1.set(1)
mass2 = IntVar()
mass2.set(1)

label1 = ttk.Label(text='choose angle of inner pendulum')
label2 = ttk.Label(text='choose angle of outside pendulum')
label3 = ttk.Label(text='choose display mode')
label4 = ttk.Label(text='set the mass of the first pendulum')
label5 = ttk.Label(text='set the mass of the second pendulum')


inner_mass = ttk.Spinbox(root, textvariable=mass1, from_=1, to=20, increment=0.1)
outer_mass = ttk.Spinbox(root, textvariable=mass2, from_=1, to=20, increment=0.1)
inner = ttk.Spinbox(root, textvariable=angle1, from_=-180, to=180, increment=0.1)
outer = ttk.Spinbox(root, textvariable=angle2, from_=-180, to=180, increment=0.1)
setmode = ttk.Combobox(root, textvariable=mode)
setmode['values'] = ['pendulum', 'path']
run = ttk.Button(root, text='go go go!',
                 command=lambda: pendulum_rk4(np.radians(90-angle1.get()),
                                              np.radians(90-angle2.get()), 0, 0, mode.get()))
title.grid(row=0, column=0, columnspan=3)
label1.grid(row=1, column=0, padx=10)
label2.grid(row=1, column=1, padx=10)
label3.grid(row=1, column=2, padx=10)
inner.grid(row=2, column=0)
outer.grid(row=2, column=1)
setmode.grid(row=2, column=2)
label4.grid(row=3, column=0)
label5.grid(row=3, column=2)
inner_mass.grid(row=4, column=0)
outer_mass.grid(row=4, column=2)
run.grid(row=5, column=1, pady=25)



root.mainloop()
