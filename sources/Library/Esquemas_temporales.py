from scipy.optimize import newton

""" Esquemas temporales """

""" Esquemas numericos para la resolucion de EDOs (Ecuaciones diferenciales) 
        Argumentos: f  = ecuacion diferencial -> dx/dt = f
                    xn = valor de x en el tiempo t
                    dt = paso de tiempo
                    t  = valor del tiempo en ese instante
                    
        Output:     valor de x en t + dt, orden del esquema                 """


""" Euler explicito """

def Euler(f,xn,dt,t):
    
    return [xn + dt*f(xn,t), 1]     # funcion dependiente de x_n para obtener x_n+1


""" Euler implicito """

def Euler_implicito(f,xn,dt,t): 

    def f_xn1(X):
          return X - xn - dt*f(X,t)      # funcion dependiente de x_n+1

    return [newton(f_xn1,xn), 1]


""" Crank - Nicolson """

def CN(f,xn,dt,t):
    
    def f_xn_xn1(X):
         return  X - (xn + (dt/2)*f(xn,t)) - (dt/2)*f(X,t+dt)    # funcion dependiente de x_n+1 y xn
    
    return [newton(f_xn_xn1,xn), 2]


""" Runge-Kutta (RK) """

def RK4(f,xn,dt,t):
    
    return [xn + (dt/6)*(f(xn,t) + 2*f(xn + f(xn,t)/2, t+dt/2) + 2*f(xn + f(xn + f(xn,t)/2, t+dt/2)/2, t+dt/2) + f(xn + f(xn + f(xn + f(xn,t)/2, t+dt/2)/2, t+dt/2), t+dt)), 4]
