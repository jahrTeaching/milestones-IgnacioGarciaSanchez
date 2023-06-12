from scipy.optimize import newton
from numpy import  zeros, linspace, log10, float64, reshape, matmul, shape
from numpy.linalg import norm

""" ODEs : Ecuaciones diferenciales ordinarias """

""" En esta librería se encuantran:
     - Esquemas temporales para la resolución de ecuaciones diferenciales ordinarias
     - El Problema de Cauchy o de valor inicial
     - Propiedades de las ODES: Error de Richardson, Ratio de convergencia, región de estabilidad """


# --------------------------------------------------------------------------------------------------- #


""" Esquemas temporales """

""" Esquemas numericos para la resolucion de EDOs (Ecuaciones diferenciales) 
        Argumentos: 
                    f  = ecuacion diferencial -> dx/dt = f
                    xn = valor de x en el tiempo t
                    dt = paso de tiempo
                    t  = valor del tiempo en ese instante
                    
        Output:     valor de x en t + dt, orden del esquema                 """

""" Euler explicito """

def Euler(f,xn,dt,t,**arg):
    Euler.__name__ = "Euler"
    
    return xn + dt*f(xn,t,**arg), 1     # funcion dependiente de x_n para obtener x_n+1


""" Euler implicito """

def Euler_implicito(f,xn,dt,t,**arg): 
    Euler_implicito.__name__ = "Euler Implicito"
    
    def f_xn1(X):
          return X - xn - dt*f(X,t,**arg)      # funcion dependiente de x_n+1

    return newton(f_xn1,xn,maxiter=200), 1


""" Crank - Nicolson """

def CN(f,xn,dt,t,**arg):
    CN.__name__ = "Crank-Nicolson"
    
    def f_xn_xn1(X):
         return  X - (xn + (dt/2)*f(xn,t,**arg)) - (dt/2)*f(X,t+dt,**arg)    # funcion dependiente de x_n+1 y xn
    
    return newton(f_xn_xn1,xn,maxiter=200), 2


""" Salto de rana (leap frog) """

def LF(f,xn,dt,t,**arg): 
    LF.__name__ = "Leap-Frog"

    return xn + dt*f(xn + (dt/2)*f(xn,t),t,**arg), 2


""" Runge-Kutta (RK) """

# Grado 4
def RK4(f,xn,dt,t,**arg):
    RK4.__name__ = "Runge-Kutta 4"
    
    return xn + (dt/6)*(f(xn,t,**arg) + 2*f(xn + dt*f(xn,t,**arg)/2, t+dt/2,**arg) + 2*f(xn + dt*f(xn + dt*f(xn,t,**arg)/2, t+dt/2,**arg)/2, t+dt/2,**arg) + f(xn + dt*f(xn + dt*f(xn + dt*f(xn,t,**arg)/2, t+dt/2,**arg)/2, t+dt/2,**arg), t+dt,**arg)), 4

# Grado 5
def RK5(f,xn,dt,t):
    RK5.__name__ = "Runge-Kutta 5"
    
    k1 = dt*f(xn,t)
    k2 = dt*f(xn + k1/4,t + dt/4)
    k3 = dt*f(xn + 3*k1/32 + 9*k2/32,t + 3*dt/8)
    k4 = dt*f(xn + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197,t + 12*dt/13)
    k5 = dt*f(xn + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104,t + dt)
    
    return xn + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5

# Grado 8 -En construcción
def RK8(f,xn,dt,t):
    RK8.__name__ = "Runge-Kutta 8"
    
    return xn + (dt/6)*(f(xn,t) + 2*f(xn + f(xn,t)/2, t+dt/2) + 2*f(xn + f(xn + f(xn,t)/2, t+dt/2)/2, t+dt/2) + f(xn + f(xn + f(xn + f(xn,t)/2, t+dt/2)/2, t+dt/2), t+dt)), 8

# Grado 9 -En construcción
def RK9(f,xn,dt,t):
    RK9.__name__ = "Runge-Kutta 9"
    
    return xn + (dt/6)*(f(xn,t) + 2*f(xn + f(xn,t)/2, t+dt/2) + 2*f(xn + f(xn + f(xn,t)/2, t+dt/2)/2, t+dt/2) + f(xn + f(xn + f(xn + f(xn,t)/2, t+dt/2)/2, t+dt/2), t+dt)), 9


""" Runge-Kutta Embebido """
def RKE(f,xn,dt,t):
    RKE.__name__ = "RK45"
    
    eps = 1e-10
    
    V1 = RK45("First", f,xn,dt,t)
    V2 = RK45("Second", f,xn,dt,t)

    (a, b, bs, c, q, Ne) = Butcher_array()
    
    h = min(dt, Step_size(V1 - V2, eps, min(q), dt))
    
    N = int( dt / h ) + 1
    h = dt / N
    
    V1 = xn; V2 = xn
    
    for i in range(N):
        time = t + i * dt / N
        V1 = V2
        
        V2 = RK45("First",f,V1,h,time)
    
    U2 = V2
    
    ierr = 0
    
    return U2, [4,5]

""" RK para embebido """
def RK45(tag,f,xn,dt,t):
    
    (a, b, bs, c, q, Ne) = Butcher_array()
    
    N = len(xn)
    k = zeros( [Ne, N] )
    
    aux = f( xn, t + c[0]*dt )
    
    k[0:1,:] = reshape(aux,(1,len(xn)))
    
    
    if tag == "First":
        
        for i in range(1,Ne):
            Up = xn
            
            for j in range(i):

                aux = dt * a[i,j]*k[j,:]
                Up[:,0:1] = Up[:,0:1] + reshape(aux,(len(xn),1))
                
            
            k[i:i+1,:] =  reshape(f( Up, t + c[i]*dt ),(1,len(xn)))
        
        U2 = xn + reshape(dt * matmul(b,k),(len(xn),1))
        
    elif tag == "Second":
        
        for i in range(1,Ne):
            Up = xn
            
            for j in range(i):
                aux = dt * a[i,j]*k[j,:]
                Up[:,0:1] = Up[:,0:1] + reshape(aux,(len(xn),1))
                
            k[i:i+1,:] =  reshape(f( Up, t + c[i]*dt ),(1,len(xn)))
            
        U2 = xn + reshape(dt * matmul(bs,k),(len(xn),1))
    
    return U2


# --------------------------------------------------------------------------------------------------- #


""" Problema de Cauchy o problema de valores iniciales """

""" Resolución del problema de Cauchy o de valor inicial de una ODE mediante un esquema temporal.
        Argumentos: 
                    f_et = funcion del esquema temporal.
                    f    = ecuacion diferencial -> dx/dt = f
                    xo   = valor inicial (to) de x.
                    t    = vector de tiempos (ti) -> len(t) = numero de pasos
                    output_q = variable de control para el return de q.
                               (True or False - valor por defecto = True)
                    
        Output:     
                    matriz de valores de x en el tiempo.
                           (filas = xi, columnas = ti)                             
                    q    = grado del esquema temporal, puede omitirse con 'output_q'. """


def Cauchy(f_et,f,xo,t,output_q=True,**arg): 

     x = zeros((len(xo),len(t))) 

     x[:,0:1] = xo[:]

     for i in range(1,len(t)):
        
        x[:,i:i+1], q = f_et(f,x[:,i-1:i],t[i]-t[i-1],t[i],**arg)

     if output_q: 
          return x, q
     else:
          return x


# --------------------------------------------------------------------------------------------------- #


""" Análisis de métodos numéricos : Convergencia, Desviación, Estabilidad """

"""     Diferentes funciones para caracterizar los métodos numéricos.     """


""" Extrapolación de Richardson - Error o desviación del método numérico 
        Argumentos: 
                    f_et  = funcion del esquema temporal 
                    f     = ecuacion diferencial -> dx/dt = f
                    xo    = valor inicial (to) de x
                    t1,t2 = vectores de tiempos (ti) -> len(ti) = numero de pasos
                            donde len(t2) > len(t1) para mismos t0 y tf
                    
        Output:     Error                                                   """


def Richardson(f_et,f,xo,t1,t2=[]):
    
    if not t2.any(): t2 = linspace(t1[0],t1[-1],2*len(t1))   # valor por defecto de t2
    
    C = (len(t2)/len(t1)) # C = dt1/dt2 , como ambos tienen mismo to y tf y aplicando dt = (tf-to)/n siendo n = numero de pasos -> C = n1/n2
    
    m = 1
    while m*C != int(m*C): m+=1
    
    U1, q = Cauchy(f_et,f,xo,t1)
    U2, q2 = Cauchy(f_et,f,xo,t2)
    
    return (U2[:,::int(m*C)] - U1[:,::m])/(1 - 1/C**q)


""" Tasa de convergencia del método numérico 
        Argumentos: 
                    f_et  = funcion del esquema temporal 
                    f     = ecuacion diferencial -> dx/dt = f
                    xo    = valor inicial (to) de x
                    t     = vector de tiempo -> len(t) = numero de pasos
                    imax  = numero de iteraciones, similar al parámetro maxiter de algunas funciones
                    
        Output:     
                    log_N = vector de logaritmos de los numero de pasos
                    log_E = vector de logaritmos de los errores                                  """

def Convergence_rate(f_et,f,xo,t,imax=[]):
    
    if imax:
        
        log_N = []
        log_E = []

        for i in range(1,imax):
            
            U1 = Cauchy(f_et,f,xo,linspace(t[0],t[-1],i*len(t)),output_q=False)
            U2 = Cauchy(f_et,f,xo,linspace(t[0],t[-1],(i+1)*len(t)),output_q=False)
            
            log_N.append(log10((i+1)*len(t)))
            log_E.append(log10(norm((U2[0:2,-1] - U1[0:2,-1]))))
            print("  Iteración " + str(i) + " de la convergencia.")
    else:
        
        log_N = [log10(2*len(t))]
        U1 = Cauchy(f_et,f,xo,linspace(t[0],t[-1],len(t)),output_q=False)
        U2 = Cauchy(f_et,f,xo,linspace(t[0],t[-1],2*len(t)),output_q=False)
        E2 = log10(norm((U2[0:2,-1] - U1[0:2,-1])))
        log_E = [E2]
        E1 = 0; i = 1

        while (E2 - E1) > 1e-2: # a cierta tolerancia
            
            i+=1
            U1 = Cauchy(f_et,f,xo,linspace(t[0],t[-1],i*len(t)),output_q=False)
            U2 = Cauchy(f_et,f,xo,linspace(t[0],t[-1],(i+1)*len(t)),output_q=False)
                        
            E1 = E2
            E2 = log10(norm((U2[0:2,-1] - U1[0:2,-1])))
            log_N.append(log10((i+1)*len(t)))
            log_E.append(E2)
            print("  Iteración " + str(i) + " de la convergencia.")
        
    return log_N, log_E


""" Region de estabilidad del método numérico 
        Argumento: f_et  = funcion del esquema temporal 
                    
        Output:                                       """

def Stability_region(fet): 

    x = linspace(-5,5,100)
    y = linspace(-5,5,100)
    SR = zeros((100, 100))

    for i in range(0,100): 
      for j in range(0,100):

        r, q = fet(lambda u, t: complex(x[i], y[j])*u, 1, 1, 0)
        SR[i, j] = abs(r) 

    return SR 


""" Bucher Array """

def Butcher_array():
    q = [5,4]
    Ne = 7 

    a = zeros( [Ne, Ne-1] )
    b = zeros( [Ne] )
    bs = zeros( [Ne] )
    c = zeros( [Ne] )
    
    c[:] = [ 0., 1./5, 3./10, 4./5, 8./9, 1., 1. ]

    a[0,:] = [          0.,           0.,           0.,         0.,           0.,     0. ]
    a[1,:] = [      1./5  ,           0.,           0.,         0.,           0.,     0. ]
    a[2,:]= [      3./40 ,        9./40,           0.,         0.,           0.,     0. ]
    a[3,:] = [     44./45 ,      -56./15,        32./9,         0.,           0.,     0. ]
    a[4,:] = [ 19372./6561, -25360./2187,  64448./6561,  -212./729,           0.,     0. ]
    a[5,:] = [  9017./3168,    -355./33 ,  46732./5247,    49./176, -5103./18656,     0. ]
    a[6,:]= [    35./384 ,           0.,    500./1113,   125./192, -2187./6784 , 11./84 ]

    b[:]  = [ 35./384   , 0.,   500./1113,  125./192,  -2187./6784  ,  11./84  ,     0.]
    bs[:] = [5179./57600, 0., 7571./16695,  393./640, -92097./339200, 187./2100, 1./40 ]
    
    return (a, b, bs, c, q, Ne)


""" Step size """
def Step_size(dU, eps, q, dt):
    
    if( norm(dU) > eps ):
        size =  dt *( eps / norm(dU) )**( 1 / ( q + 1) ) # Step_size
    else:
        size = dt # Step_size
        
    return size