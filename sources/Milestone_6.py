

from Library.ODEs import Cauchy, RKE
from Library.Mecanica_Orbital import bodies3, Lagrange, Stb_Lagrange

from numpy import array, linspace, zeros, around, shape
from random import random

import matplotlib.pyplot as plt

# Condiciones de la simulación
to = 0     # [s] tiempo de inicio
tf = 100   # [s] tiempo final
dt = 0.1   # [s] delta de t - paso de tiempo

n = int((tf-to)/dt)   # numero de pasos
t = linspace(to,tf,n) # vector tiempo

# Puntos de Lagrange
x0 = zeros([5,4])

x0[0,:] = array([0.8, 0.6, 0, 0])
x0[1,:] = array([0.8, -0.6, 0, 0])
x0[2,:] = array([-0.1, 0, 0, 0])
x0[3,:] = array([0.1, 0, 0, 0])
x0[4,:] = array([1.01, 0, 0, 0])

xp = Lagrange(x0, t)

# Estabilidad 

xp_stb = zeros(4)

for i in range(5):
    xp_stb[:2] = xp[i, :]
    stb = around(Stb_Lagrange(xp_stb,t), 5)

# Orbitas
x0_Lp = zeros([5,4])
aux0_Lp = zeros([4,1])

eps = 1e-3

for i in range(5):

    x0_Lp[i, :2] = xp[i, :] + eps
    x0_Lp[i, 2:] = eps

    aux0_Lp[:,0] = x0_Lp[i,:]
    print("1")
    x_Lp = Cauchy(RKE,bodies3,aux0_Lp,t,output_q=False)


# Plots

LP = ["L4", "L5", "L3", "L1", "L2"]

for j in range(5):
        
    fig, ax = plt.subplots( figsize = (10, 10) )

    for i in range( len(xp) ):
        ax.plot( xp[i, 0],xp [i, 1], 'o', label = LP[i] )
        
    ax.plot( x_Lp[0,:], x_Lp[1,:], color = 'b', label = "Orbit" )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Lagrange Points of Earth-Moon System and Orbit arround " + LP[j])
    plt.legend()
    ax.grid()
    print(j)
    if LP[j] == "L4" or LP[j] == "L5": 
        fig, ay = plt.subplots( figsize = (7,7) )
        ay.plot( xp[j, 0],xp [j, 1], 'o', label = LP[j] )
        ay.plot( x_Lp[0,:], x_Lp[1,:], color = 'b', label = "Orbit"  )
        ay.set_xlabel("x")
        ay.set_ylabel("y")
        ay.set_title("Orbit arround" + LP[j])
        plt.legend()
        ay.grid()


