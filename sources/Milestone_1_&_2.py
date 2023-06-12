from numpy import array, zeros, linspace
from Library.Plots import Plot
from Library.Mecanica_Orbital import Kepler, v_orbital
from Library.ODEs import Euler, Euler_implicito, CN, RK4, Cauchy


""" Milestone 1 y 2: Propagadores y comparación """

""" Simulación de distintos propagadores: Euler explícito e implícito, Crank Nicolson, RK, ... """ 

# Posición y velocidad inicial (componentes en el plano)
r0 = array([0,6978.14])         # [km] posicion inicial
v0 = v_orbital(r0)*array([1,0]) # [km/s] velocidad inicial

xo = array([[r0[0]],[r0[1]],[v0[0]],[v0[1]]]) # vector de estado inicial

# Condiciones de la simulación
to = 0        # [s] tiempo de inicio
tf = 100000   # [s] tiempo final
dt = 0.1        # [s] delta de t - paso de tiempo

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


for i in range(1,n):
    
    # Euler o Euler explicito
    r_E[:,i:i+1], qE = Euler(Kepler,r_E[:,i-1:i],dt,tf)
    
    # Euler implicito
    r_Ei[:,i:i+1], qEi = Euler_implicito(Kepler,r_Ei[:,i-1:i],dt,tf)
    
    # Crank Nicolson
    r_CN[:,i:i+1], qCN = CN(Kepler,r_CN[:,i-1:i],dt,tf)
    
    # Runge Kutta 4
    r_RK4[:,i:i+1], qRK4 = RK4(Kepler,r_RK4[:,i-1:i],dt,tf)

# Inicialización del vector tiempo
t = linspace(to,tf,n+1)    # para introducir en una sola variable: t0, tf, dt y n

# Cauchy con Euler
r_CE = Cauchy(Euler,Kepler,xo,t,output_q=False)

# Cauchy con Euler implicito
r_CEi = Cauchy(Euler_implicito,Kepler,xo,t,output_q=False)

# Cauchy con Crank Nicolson
r_CCN = Cauchy(CN,Kepler,xo,t,output_q=False)

# Cauchy con Runge Kutta 4
r_CRK4 = Cauchy(RK4,Kepler,xo,t,output_q=False)


""" PLOTS """
Plot(r_E[0,:],r_E[1,:], Title = "Euler")

Plot(r_Ei[0,:],r_Ei[1,:], Title = "Euler implicito")

Plot(r_CN[0,:],r_CN[1,:], Title = "Crank Nicolson")

Plot(r_RK4[0,:],r_RK4[1,:], Title = "Runge Kutta 4")

Plot(r_CE[0,:],r_CE[1,:], Title = "Cauchy con Euler")

Plot(r_CEi[0,:],r_CEi[1,:], Title = "Cauchy con Euler implicito")

Plot(r_CCN[0,:],r_CCN[1,:], Title = "Cauchy con Crank Nicolson")

Plot(r_CRK4[0,:],r_CRK4[1,:], Title = "Cauchy con Runge Kutta 4")

