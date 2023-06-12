from numpy import array, linspace, transpose, zeros
from Library.Plots import Figure, Plot, ShowPlot, Legend, Grid, Contour
from Library.ODEs import Euler, Euler_implicito, CN, RK4, LF, Cauchy, Stability_region


""" Milestone 4: Integración de un oscilador lineal + Estudio de las regiones de estabilidad """

""" Estudio con distintos propagadores (Euler explícito e implícito, Crank Nicolson, RK, ...) """ 

# --------------------------------------------------------------------------------------------------- #

""" Integración de un oscilador lineal """

# Posición y velocidad inicial
x0 = array([[0],[1]])

# Condiciones de la simulación
to = 0     # [s] tiempo de inicio
tf = 100   # [s] tiempo final
dt = 0.1   # [s] delta de t - paso de tiempo

n = int((tf-to)/dt)   # numero de pasos
t = linspace(to,tf,n) # vector tiempo


# Función del oscilador lineal
def f_oscilator(x,t):  # d2(x)/dt2 + x = 0

    return array([x[1],-x[0]])


# Cauchy con Euler
x_E = Cauchy(Euler,f_oscilator,x0,t,output_q=False)

# Cauchy con Euler implicito
x_Ei = Cauchy(Euler_implicito,f_oscilator,x0,t,output_q=False)

# Cauchy con Crank Nicolson
x_CN = Cauchy(CN,f_oscilator,x0,t,output_q=False)

# Cauchy con Runge Kutta 4
x_RK4 = Cauchy(RK4,f_oscilator,x0,t,output_q=False)

# Cauchy con Leap Frog
x_LF = Cauchy(LF,f_oscilator,x0,t,output_q=False)


""" PLOTS """
Plot(x_E[0,:],x_E[1,:], Title = "Euler")

Plot(x_Ei[0,:],x_Ei[1,:], Title = "Euler implicito")

Plot(x_CN[0,:],x_CN[1,:], Title = "Crank Nicolson")

Plot(x_RK4[0,:],x_RK4[1,:], Title = "Runge Kutta 4")

Plot(x_LF[0,:],x_LF[1,:], Title = "Leap Frog")


# --------------------------------------------------------------------------------------------------- #

""" Estudio de las regiones de estabilidad """

SR_E = Stability_region(Euler)
SR_Ei = Stability_region(Euler_implicito)
SR_CN = Stability_region(CN)
SR_RK4 = Stability_region(RK4)
SR_LF = Stability_region(LF)

""" PLOTS """

x = linspace(-5,5,100)
y = linspace(-5,5,100)

# Euler
Figure()
px = Plot(x,zeros(100),Show="no",c='k')  # eje x
py = Plot(zeros(100),y,Show="no",c='k')  # eje y
Grid()
Contour(x, y, transpose(SR_E), ctype = "lines", Show="no", levels = [0, 1], colors = ['r'], linewidths = 2)
Contour(x, y, transpose(SR_E), ctype = "filled", levels = [0, 1], colors =['b'])     

# Euler implicito
Figure()
px = Plot(x,zeros(100),Show="no",c='k')  # eje x
py = Plot(zeros(100),y,Show="no",c='k')  # eje y
Grid()
Contour(x, y, transpose(SR_Ei), ctype = "lines", Show="no", levels = [0, 1], colors = ['r'], linewidths = 2)
Contour(x, y, transpose(SR_Ei), ctype = "filled", levels = [0, 1], colors =['b'])     

# Crank Nicolson
Figure()
px = Plot(x,zeros(100),Show="no",c='k')  # eje x
py = Plot(zeros(100),y,Show="no",c='k')  # eje y
Grid()
Contour(x, y, transpose(SR_CN), ctype = "lines", Show="no", levels = [0, 1], colors = ['r'], linewidths = 2)
Contour(x, y, transpose(SR_CN), ctype = "filled", levels = [0, 1], colors =['b'])     

# Runge Kutta
Figure()
px = Plot(x,zeros(100),Show="no",c='k')  # eje x
py = Plot(zeros(100),y,Show="no",c='k')  # eje y
Grid()
Contour(x, y, transpose(SR_RK4), ctype = "lines", Show="no", levels = [0, 1], colors = ['r'], linewidths = 2)
Contour(x, y, transpose(SR_RK4), ctype = "filled", levels = [0, 1], colors =['b'])     

# Leap Frog
Figure()
px = Plot(x,zeros(100),Show="no",c='k')  # eje x
py = Plot(zeros(100),y,Show="no",c='k')  # eje y
Grid()
Contour(x, y, transpose(SR_LF), ctype = "lines", Show="no", levels = [0, 1], colors = ['r'], linewidths = 2)
Contour(x, y, transpose(SR_LF), ctype = "filled", levels = [0, 1], colors =['b'])     

