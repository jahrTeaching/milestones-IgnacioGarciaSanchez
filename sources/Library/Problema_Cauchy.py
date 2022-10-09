from numpy import  zeros

""" Problema de Cauchy o problema de valores iniciales """

""" Resolución del problema de Cauchy o de valor inicial de una ODE mediante 
    un esquema temporal.
        Argumentos: f    = ecuacion diferencial -> dx/dt = f
                    xo   = valor inicial (to) de x
                    t    = vector de tiempos (ti) -> len(t) = numero de pasos
                    f_et = funcion del esquema temporal 
                    
        Output:     matriz de valores de x en el tiempo
                    (filas = xi, columnas = ti)                             """

def Cauchy(f,xo,t,f_et): 

     x = zeros((len(xo),len(t))) 

     x[:,0:1] = xo

     for i in range(1,len(t)):

        x[:,i:i+1] = f_et(f,x[:,i-1:i],t[i]-t[i-1],t[i]) 

     return x
 
    
