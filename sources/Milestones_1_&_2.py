from numpy import array, zeros, linspace
from Library.Plot import Plot
from Library.Kepler import Kepler
from Library.Esquemas_temporales import Euler, Euler_implicito, CN, RK4
from Library.Problema_Cauchy import Cauchy

""" Milestone 1 y 2: Propagadores y comparación """

""" Simulación de distintos propagadores (Euler explícito e implícito, 
    Crank Nicolson, RK, ...)                                                   """ 

# Posición y velocidad inicial (componentes en el plano)
xo = array([[0],[6978.14],[7],[0]])

# Condiciones de la simulación
to = 0           # [s] tiempo de inicio
tf = 10000       # [s] tiempo final
dt = 5           # [s] delta de t - paso de tiempo

n = int((tf-to)/dt) # numero de pasos

# Inicialización de matrices
r_E = zeros((4,n))
r_Ei = zeros((4,n))
r_CN = zeros((4,n))
r_RK4 = zeros((4,n))
r_E[:,0:1] = xo
r_Ei[:,0:1] = xo
r_CN[:,0:1] = xo
r_RK4[:,0:1] = xo


for i in  range(1,n):
    
    # Euler o Euler explicito
    r_E[:,i:i+1] = Euler(Kepler,r_E[:,i-1:i],dt,tf)
    
    # Euler implicito
    r_Ei[:,i:i+1] = Euler_implicito(Kepler,r_Ei[:,i-1:i],dt,tf)
    
    # Crank Nicolson
    r_CN[:,i:i+1] = CN(Kepler,r_CN[:,i-1:i],dt,tf)
    
    # Runge Kutta 4
    r_RK4[:,i:i+1] = RK4(Kepler,r_RK4[:,i-1:i],dt,tf)

# Inicialización del tiempo
t = linspace(to,tf,n+1)    # para introducir en una sola variable: t0, tf, dt y n

# Cauchy con Euler
r_CE = Cauchy(Kepler,xo,t,Euler)

# Cauchy con Euler implicito
r_CEi = Cauchy(Kepler,xo,t,Euler_implicito)

# Cauchy con Crank Nicolson
r_CCN = Cauchy(Kepler,xo,t,CN)

# Cauchy con Runge Kutta 4
r_CRK4 = Cauchy(Kepler,xo,t,RK4)

""" PLOTS """
Plot(r_E[0,:],r_E[1,:],Title="Euler")

Plot(r_Ei[0,:],r_Ei[1,:],Title="Euler implicito")

Plot(r_CN[0,:],r_CN[1,:],Title="Crank Nicolson")

Plot(r_RK4[0,:],r_RK4[1,:],Title="Runge Kutta 4")

Plot(r_CE[0,:],r_CE[1,:],Title="Cauchy con Euler")

Plot(r_CEi[0,:],r_CEi[1,:],Title="Cauchy con Euler implicito")

Plot(r_CCN[0,:],r_CCN[1,:],Title="Cauchy con Crank Nicolson")

Plot(r_CRK4[0,:],r_CRK4[1,:],Title="Cauchy con Runge Kutta 4")

