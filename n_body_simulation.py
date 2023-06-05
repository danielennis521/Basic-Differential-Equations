import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from tkinter import ttk

dt = 0.01
G = 6.67e-2
#G = 6.67e-11


def on_closing():
    plt.close()
    root.destroy()


def newtons_law(p1, p2):
    r = np.sqrt((p2.pos[0] - p1.pos[0])**2 + (p2.pos[1] - p1.pos[1])**2)
    f = ((p2.pos[0] - p1.pos[0])/r, (p2.pos[1] - p1.pos[1])/r)
    F = G * p1.mass * p2.mass / r**2
    return f[0]*F, f[1]*F


class point:        # an object to be simulated 

    def __init__(self, pos, vel, mass):
        if len(pos)==len(vel):
            self.pos = pos
            self.vel = vel
            self.mass = mass
            self.dim = len(pos)



class nbody_sim:    # object that takes an array of points and uses them to simulate the systems gravity
    
    def __init__(self, points):
        self.points = points
        self.time = 0.0
        self.n = len(points)
        self.forces = np.zeros((self.n, 2)) 

        self.gravity()
        temp = []
        for i in range(self.n):
            x = self.points[i].pos[0] + dt*self.points[i].vel[0] + 0.5*self.forces[i][0]*dt**2/self.points[i].mass
            y = self.points[i].pos[1] + dt*self.points[i].vel[1] + 0.5*self.forces[i][1]*dt**2/self.points[i].mass
            vx = self.points[i].vel[0] + dt*self.forces[i][0]/self.points[i].mass
            vy = self.points[i].vel[1] + dt*self.forces[i][1]/self.points[i].mass

            p = point([x, y], [vx, vy], self.points[i].mass)
            temp.append(p)
        
        self.prev = self.points
        self.points = temp


    def gravity(self):   # determine the gravitational force acting on the particles
        for f in self.forces:
            f[0] = 0.0
            f[1] = 0.0

        for i in range(self.n//2):
            for j in range(i+1, self.n):
                new_f = newtons_law(self.points[i], self.points[j])
                self.forces[i][0] += new_f[0]
                self.forces[i][1] += new_f[1]
                self.forces[j][0] -= new_f[0]
                self.forces[j][1] -= new_f[1]



    def update_system(self):
        self.gravity()

        for i in range(self.n):
            x = 2.0*self.points[i].pos[0] - self.prev[i].pos[0] + self.forces[i][0]*dt**2/self.points[i].mass
            y = 2.0*self.points[i].pos[1] - self.prev[i].pos[1] + self.forces[i][1]*dt**2/self.points[i].mass
            self.prev[i].pos = self.points[i].pos
            self.points[i].pos = [x, y]


root = Tk()

# define variables set in window
mass = IntVar()
mass.set(1)
xpos = IntVar()
xpos.set(0)
ypos = IntVar()
ypos.set(0)
xvel = IntVar()
xvel.set(0)
yvel = IntVar()
yvel.set(0)

    
# build the window
title = Label(root, text='The N-body problem! (in 2-D)')
title.grid(row=0, column=0, columnspan=20)

label1 = ttk.Label(text='Create an object')
label1.grid(row=1, column=0, columnspan=12, padx=10, pady=10)

label2 = ttk.Label(text='set mass')
label2.grid(row=2, column=0, columnspan=4, padx=10)
set_mass = ttk.Spinbox(root, textvariable=mass, from_=1, to=20, increment=0.1)
set_mass.grid(row=3, column=1, columnspan=2, padx=10)

label3 = ttk.Label(text='set initial position')
label3.grid(row=2, column=4, columnspan=4, padx=10)

label4 = ttk.Label(text='set the initial velocity')
label4.grid(row=2, column=8, columnspan=4, padx=10)

label5 = ttk.Label(text='x pos')
label5.grid(row=3, column=4, padx=10)
set_xpos = ttk.Spinbox(root, textvariable=xpos, from_=-10, to=10, increment=0.01)
set_xpos.grid(row=3, column=5, padx=10)

label6 = ttk.Label(text='y pos')
label6.grid(row=3, column=6, padx=10)
set_ypos = ttk.Spinbox(root, textvariable=ypos, from_=-10, to=10, increment=0.01)
set_ypos.grid(row=3, column=7, padx=10)

label7 = ttk.Label(text='x vel')
label7.grid(row=3, column=8, padx=10)
set_xvel = ttk.Spinbox(root, textvariable=xvel, from_=-10, to=10, increment=0.01)
set_xvel.grid(row=3, column=9, padx=10)

label8 = ttk.Label(text='y vel')
label8.grid(row=3, column=10, padx=10)
set_yvel = ttk.Spinbox(root, textvariable=yvel, from_=-10, to=10, increment=0.01)
set_yvel.grid(row=3, column=11, padx=10)

add_point = ttk.Button(root, text='Add the object')
add_point.grid(row=4, column=6, columnspan=2, pady=20)

run_sim = ttk.Button(root, text='GO GO GO!!!')
run_sim.grid(row=5, column=6, columnspan=2, pady=20)

# frame = 12
# p1 = point([0.0, 5.0], [20.0, -5.0], 1)
# p2 = point([0.0, 0.0], [0.0, 0.0], 100000)
# sim = nbody_sim([p1, p2])


# for i in range(1000):
#     plt.xlim(-frame, frame)
#     plt.ylim(-frame, frame)
#     plt.plot(sim.points[0].pos[0], sim.points[0].pos[1], 'ro', color='blue', markersize=4)
#     plt.plot(sim.points[1].pos[0], sim.points[1].pos[1], 'ro', color='red', markersize=10)
#     plt.pause(0.001)
#     plt.cla()

#     sim.update_system()

root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()
