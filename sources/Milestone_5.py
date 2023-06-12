from numpy import array, linspace, concatenate, shape, reshape, ceil, amax, absolute, transpose
from Library.Plots import Figure, Plot, animate_plot, ShowPlot, Legend, Fig3D, Xlabel, Ylabel
from Library.ODEs import Euler, Euler_implicito, CN, RK4, LF, Cauchy
from Library.Mecanica_Orbital import Nbodies, random_Ub, solarsystem
from matplotlib.animation import FuncAnimation


""" Milestone 5: Problema de los N cuerpos """

""" Estudio con distintos propagadores (Euler explícito e implícito, Crank Nicolson, RK, ...) """ 

""" Resolución con 4 cuerpos"""
nb = 4

r1 = array([[10],[0],[0]])
v1 = array([[0],[-0.3],[0]])

r2 = array([[-10],[0],[0]])
v2 = array([[0],[0.3],[0]])

r3 = array([[0],[10],[0]])
v3 = array([[0.3],[0],[0]])

r4 = array([[0],[-10],[0]])
v4 = array([[-0.3],[0],[0]])

U0 = concatenate((r1,v1,r2,v2,r3,v3,r4,v4),axis=0)

tf = 300
dt = 0.01

n = int(tf/dt)
t  = linspace(0,tf,n)

print("Comienzo de la propagación de " + str(nb) + " cuerpos mediante diferentes métodos numéricos")

U_E = Cauchy(Euler,Nbodies,U0,t,output_q=False,Nb=nb)  # [Nb*6 x n]
U_Ei = Cauchy(Euler_implicito,Nbodies,U0,t,output_q=False,Nb=nb)
U_CN = Cauchy(CN,Nbodies,U0,t,output_q=False,Nb=nb)
U_RK4 = Cauchy(RK4,Nbodies,U0,t,output_q=False,Nb=nb)
U_LF = Cauchy(LF,Nbodies,U0,t,output_q=False,Nb=nb)

print("Propagación de " + str(nb) + " cuerpos finalizada")

""" PLOTS """
# Euler
for i in range(0,nb):
    Plot(U_E[i*6,0],U_E[i*6+1,0],Show="no",Form="r+")
    Plot(U_E[i*6,:],U_E[i*6+1,:],Show="no")
Plot(0,0,Title="Euler",Form="k.")

# Euler Implicito
for i in range(0,nb):
    Plot(U_Ei[i*6,0],U_Ei[i*6+1,0],Show="no",Form="r+")
    Plot(U_Ei[i*6,:],U_Ei[i*6+1,:],Show="no")
Plot(0,0,Title="Euler Implicito",Form="k.")

# Crank Nicolson
for i in range(0,nb):
    Plot(U_CN[i*6,0],U_CN[i*6+1,0],Show="no",Form="r+")
    Plot(U_CN[i*6,:],U_CN[i*6+1,:],Show="no")
Plot(0,0,Title="Crank Nicolson",Form="k.")

# Runge Kutta 4
for i in range(0,nb):
    Plot(U_RK4[i*6,0],U_RK4[i*6+1,0],Show="no",Form="r+")
    Plot(U_RK4[i*6,:],U_RK4[i*6+1,:],Show="no")
Plot(0,0,Title="Runge Kutta",Form="k.")

# Leap Frog
for i in range(0,nb):
    Plot(U_LF[i*6,0],U_LF[i*6+1,0],Show="no",Form="r+")
    Plot(U_LF[i*6,:],U_LF[i*6+1,:],Show="no")
Plot(0,0,Title="Leap Frog",Form="k.")


""" Resolución con N cuerpos"""
nb = 10

U0 = random_Ub(nb,dim=3,kp=10,kv=0.3)

tf = 100
dt = 0.01

n = int(tf/dt)  
t  = linspace(0,tf,n)

print("Comienzo de la propagación de " + str(nb) + " cuerpos mediante Euler")

U = Cauchy(Euler,Nbodies,U0,t,output_q=False,Nb=nb)  # [Nb*6 x n]

print("Propagación de " + str(nb) + " cuerpos finalizada")

""" PLOT """
# Progresivo
U_s = reshape(U,(n,nb,2,3)) 
r = reshape(U_s[:,:,0,:],(n,nb,3))

animate_plot(r)

# Foto final
F = Fig3D()
lgd = []
for i in range(0,nb):
    F.plot(U[i*6,0], U[i*6+1,0],U[i*6+2,0],"r+")
    F.plot(U[i*6,:], U[i*6+1,:],U[i*6+2,:])
    lgd.append(""); lgd.append("Cuerpo " + str(i+1))
F.plot(0,0,0,"k.")
Legend(labels=lgd)
Xlabel("X")
Ylabel("Y")
ShowPlot()
