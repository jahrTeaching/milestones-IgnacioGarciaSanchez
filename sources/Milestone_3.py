from numpy import array, linspace
from Library.Plots import Plot, Bar, Trendline, ShowPlot, Legend
from Library.Mecanica_Orbital import Kepler, v_orbital
from Library.ODEs import Euler, Euler_implicito, CN, RK4, Richardson, Convergence_rate, Cauchy


""" Milestone 3: Discursión del error y la velocidad de covergencia """

""" Estudio de distintos propagadores (Euler explícito e implícito, Crank Nicolson, RK, ...) """ 

# Posición y velocidad inicial (componentes en el plano)
r0 = array([0,6978.14])         # [km] posicion inicial
v0 = v_orbital(r0)*array([1,0]) # [km/s] velocidad inicial

xo = array([[r0[0]],[r0[1]],[v0[0]],[v0[1]]]) # vector de estado inicial

# Condiciones de la simulación
to = 0        # [s] tiempo de inicio
tf = 100000   # [s] tiempo final

dt1 = 1       # [s] delta de t - paso de tiempo
dt2 = 0.1     # [s] delta de t - paso de tiempo

t1 = linspace(to,tf,int((tf-to)/dt1))
t2 = linspace(to,tf,int((tf-to)/dt2))


# Desviaciones de Richardson de los esquemas temporales
E_Euler   = Richardson(Euler,Kepler,xo,t1,t2); print("\nError de Euler calculado.")
E_Euler_i = Richardson(Euler_implicito,Kepler,xo,t1,t2); print("Error de Euler Implicito calculado.")
E_CN      = Richardson(CN,Kepler,xo,t1,t2); print("Error de Crank-Nicolson calculado.")
E_RK4     = Richardson(RK4,Kepler,xo,t1,t2); print("Error de Runge Kutta calculado.")
print("\n ----- Errores de Richardson obtenidos ----- \n")

# Velocidad de convergencia de los esquemas temporales
log_NE, log_EE = Convergence_rate(Euler,Kepler,xo,t1,imax=10); print(" - Convergencia de Euler determinada.\n")
log_NEi, log_EEi = Convergence_rate(Euler_implicito,Kepler,xo,t1,imax=10); print(" - Convergencia de Euler Implicito determinada.\n")
log_NCN, log_ECN = Convergence_rate(CN,Kepler,xo,t1,imax=10); print(" - Convergencia de Crank-Nicolson determinada.\n")
log_NRK, log_ERK = Convergence_rate(RK4,Kepler,xo,t1,imax=10); print(" - Convergencia de Runge Kutta determinada.")
print("\n ----- Velocidades de converdencia realizadas ----- \n")

# Euler
r_E1 = Cauchy(Euler,Kepler,xo,t1,output_q=False)
r_E2 = Cauchy(Euler,Kepler,xo,t2,output_q=False)

# Euler implicito
r_Ei1 = Cauchy(Euler_implicito,Kepler,xo,t1,output_q=False)
r_Ei2 = Cauchy(Euler_implicito,Kepler,xo,t2,output_q=False)

# CN
r_CN1 = Cauchy(CN,Kepler,xo,t1,output_q=False)
r_CN2 = Cauchy(CN,Kepler,xo,t2,output_q=False)

# RK4
r_RK1 = Cauchy(RK4,Kepler,xo,t1,output_q=False)
r_RK2 = Cauchy(RK4,Kepler,xo,t2,output_q=False)


""" PLOTS """
E0 = Plot(xo[0],xo[1], Form="r+", Show="no", label="Initial point")
E1 = Plot(r_E1[0,:],r_E1[1,:], Show="no", label="t1")
E2 = Plot(r_E2[0,:],r_E2[1,:], Show="no", Title="Euler t1 vs t2", label="t2")
Legend(handles=[E0[0], E1[0], E2[0]])
ShowPlot()

Ei0 = Plot(xo[0],xo[1], Form="r+", Show="no", label="Initial point")
Ei1 = Plot(r_Ei1[0,:],r_Ei1[1,:], Show="no", label="t1")
Ei2 = Plot(r_Ei2[0,:],r_Ei2[1,:], Show="no", Title="Euler implicito t1 vs t2", label="t2")
Legend(handles=[Ei0[0], Ei1[0], Ei2[0]])
ShowPlot()

CN0 = Plot(xo[0],xo[1], Form="r+", Show="no", label="Initial point")
CN1 = Plot(r_CN1[0,:],r_CN1[1,:], Show="no", label="t1")
CN2 = Plot(r_CN2[0,:],r_CN2[1,:], Show="no", Title="CN t1 vs t2", label="t2")
Legend(handles=[CN0[0], CN1[0], CN2[0]])
ShowPlot()

RK0 = Plot(xo[0],xo[1], Form="r+", Show="no", label="Initial point")
RK1 = Plot(r_RK1[0,:],r_RK1[1,:], Show="no", label="t1")
RK2 = Plot(r_RK2[0,:],r_RK2[1,:], Show="no", Title="RK t1 vs t2", label="t2")
Legend(handles=[RK0[0], RK1[0], RK2[0]])
ShowPlot()

p1 = Plot(t1, (E_Euler[1,:]**2 + E_Euler[2,:]**2)**0.5,     Show="no", label="Euler")
p2 = Plot(t1, (E_Euler_i[1,:]**2 + E_Euler_i[2,:]**2)**0.5, Show="no", label="Euler Implicito")
p3 = Plot(t1 ,(E_CN[1,:]**2 + E_CN[2,:]**2)**0.5,           Show="no", label="CN")
p4 = Plot(t1, (E_RK4[1,:]**2 + E_RK4[2,:]**2)**0.5,         Show="no", Title="Evolución del error", label="RK")
Legend(handles=[p1[0], p2[0], p3[0], p4[0]])
ShowPlot()
       
Trendline(t1, (E_Euler[1,:]**2 + E_Euler[2,:]**2)**0.5,     Show="no", label="Euler")
Trendline(t1, (E_Euler_i[1,:]**2 + E_Euler_i[2,:]**2)**0.5, Show="no", label="Euler Implicito")
Trendline(t1, (E_CN[1,:]**2 + E_CN[2,:]**2)**0.5,           Show="no", label="CN")
Trendline(t1, (E_RK4[1,:]**2 + E_RK4[2,:]**2)**0.5,         Show="no", Title="Tendencia del error", label="RK")
Legend(handles=[p1[0], p2[0], p3[0], p4[0]])
ShowPlot()

Bar(["Euler","Euler Implicito","Crank-Nicolson","Runge-Kutta 4"],
    [(max((E_Euler[1,:]**2 + E_Euler[2,:]**2)**0.5)),
     (max((E_Euler_i[1,:]**2 + E_Euler_i[2,:]**2)**0.5)),
     (max((E_CN[1,:]**2 + E_CN[2,:]**2)**0.5)),
     (max((E_RK4[1,:]**2 + E_RK4[2,:]**2)**0.5))], Title="Error máx cometido en la posición")

Bar(["Euler","Euler Implicito","Crank-Nicolson","Runge-Kutta 4"],
    [(sum((E_Euler[1,:]**2 + E_Euler[2,:]**2)**0.5)),
     (sum((E_Euler_i[1,:]**2 + E_Euler_i[2,:]**2)**0.5)),
     (sum((E_CN[1,:]**2 + E_CN[2,:]**2)**0.5)),
     (sum((E_RK4[1,:]**2 + E_RK4[2,:]**2)**0.5))], Title="Error acumulado en la posición")

Plot(log_NE,log_EE, Title = "Velocidad de convergencia de Euler")
Plot(log_NEi,log_EEi, Title = "Velocidad de convergencia de Euler Implicito")
Plot(log_NCN,log_ECN, Title = "Velocidad de convergencia de Crank-Nicolson")
Plot(log_NRK,log_ERK, Title = "Velocidad de convergencia de Runge-Kutta 4")
