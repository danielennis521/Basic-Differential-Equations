import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from tkinter import ttk


def newtons_law(p1, p2):
    r = np.sqrt((p2.pos[0] - p1.pos[0])**2 + (p2.pos[1] - p1.pos[1])**2)
    f = ((p2.pos[0] - p1.pos[0])/r, (p2.pos[1] - p1.pos[1])/r)
    F = 6.67e-11 * p1.mass * p2.mass / r**2
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
        self.forces = np.zeros(self.n, 2) 

        self.gravity()
        temp = []
        for i in range(self.n):
            x = self.points[i].pos[0] + 0.001*self.points[i].vel[0] + 0.5*self.forces[i][0]*0.00001/self.points[i].mass
            y = self.points[i].pos[1] + 0.001*self.points[i].vel[1] + 0.5*self.forces[i][1]*0.00001/self.points[i].mass
            vx = self.points[i].vel[0] + 0.001*self.forces[i][0]/self.points[i].mass
            vy = self.points[i].vel[1] + 0.001*self.forces[i][1]/self.points[i].mass

            p = point([x, y], [vx, vy], self.points[i].mass)
            temp.append(p)
        
        self.prev = self.points
        self.points = temp


    def gravity(self):   # determine the gravitational force acting on the particles
        for f in self.forces:
            f[0] = 0.0
            f[1] = 0.0

        for i in range(self.n):
            for j in range(i+1, self.n):
                f = newtons_law(self.points[i], self.points[j])
                self.forces[i][0] += f[0]
                self.forces[i][1] += f[1]
                self.forces[j][0] -= f[0]
                self.forces[j][1] -= f[1]



    def update_system(self):
        self.gravity()

        for i in range(self.n):
            x = 2.0*self.points[i].pos[0] - self.prev[i].pos[0] + self.forces[i][0]*0.00001/self.points[i].mass
            y = 2.0*self.points[i].pos[1] - self.prev[i].pos[1] + self.forces[i][1]*0.00001/self.points[i].mass
            self.prev[i].pos = self.points[i].pos
            self.points[i].pos = [x, y]

    
    
# build the window
root = Tk()
title = Label(root, text='The N-body problem! (in 2-D)')

