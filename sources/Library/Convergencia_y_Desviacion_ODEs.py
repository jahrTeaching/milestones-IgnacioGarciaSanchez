from numpy import linspace, log10
from numpy.linalg import norm
from Library.Problema_Cauchy import Cauchy

""" Análisis de métodos numéricos : Convergencia y Desviación """

""" Diferentes funciones para caracterizar los métodos numéricos.
        Argumentos: f_et  = funcion del esquema temporal 
                    q     = orden del esquema numérico
                    f     = ecuacion diferencial -> dx/dt = f
                    xo    = valor inicial (to) de x
                    t1,t2 = vectores de tiempos (ti) -> len(t) = numero de pasos
                            donde len(t2) > len(t1)
                    
        Output:     Error                                                   """


""" Extrapolación de Richardson - Error o desviación del método numérico """

def Richardson(f_et,q,f,xo,t1,t2=[]):
    
    if not t2.any(): t2 = linspace(t1[0],t1[-1],2*len(t1))   # valor por defecto de t2
    
    C = (len(t2)/len(t1)) # C = dt1/dt2 , como ambos tienen mismo to y tf y aplicando dt = (tf-to)/n siendo n = numero de pasos -> C = n1/n2
    
    m = 1
    while m*C != int(m*C): m+=1
    
    return (Cauchy(f_et,f,xo,t2)[:,::int(m*C)] - Cauchy(f_et,f,xo,t1)[:,::m])/(1 - 1/C**q)



""" Tasa de convergencia del método numérico """

def Convergence_rate(f_et,q,f,xo,t,imax=[],show=[]):
    
    if imax:
        
        log_N = []
        log_E = []
        for i in range(1,imax):
            
            log_N.append(log10((i+1)*len(t)))
            log_E.append(log10(norm((Cauchy(f_et,f,xo,linspace(t[0],t[-1],(i+1)*len(t)))[0:2,-1] - Cauchy(f_et,f,xo,linspace(t[0],t[-1],i*len(t)))[0:2,-1]))))
            print("Iteración " + str(i) + " de la convergencia.\n")
    else:
        
        log_N = [log10(2*len(t))]
        E2 = log10(norm((Cauchy(f_et,f,xo,linspace(t[0],t[-1],2*len(t)))[0:2,-1] - Cauchy(f_et,f,xo,linspace(t[0],t[-1],len(t)))[0:2,-1])))
        log_E = [E2]
        E1 = 0; i = 1
        while (E2 - E1) > 1e-2:
            
            i+=1
            E1 = E2
            E2 = log10(norm((Cauchy(f_et,f,xo,linspace(t[0],t[-1],(i+1)*len(t)))[0:2,-1] - Cauchy(f_et,f,xo,linspace(t[0],t[-1],i*len(t)))[0:2,-1])))
            log_N.append(log10((i+1)*len(t)))
            log_E.append(E2)
            print("Iteración " + str(i) + " de la convergencia.\n")
        
    return log_N, log_E


""" Llamadas a Richarson con el orden definido """

def Richardson_D1(f_et,f,xo,t1,t2=[]):
    
    return Richardson(f_et,1,f,xo,t1,t2)

def Richardson_D2(f_et,f,xo,t1,t2):
    
    return Richardson(f_et,2,f,xo,t1,t2)

def Richardson_D4(f_et,f,xo,t1,t2):
    
    return Richardson(f_et,4,f,xo,t1,t2)


