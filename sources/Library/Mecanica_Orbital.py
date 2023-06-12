from scipy.optimize import newton
from numpy import array, zeros, shape, reshape, size, cos, sin, pi, float64
from numpy.linalg import norm, eig
from random import uniform

""" FUNCIONES ORBITALES """

""" En esta librería se encuantran todas las funciones relacionadoas con la mecánica orbital: 
    - Funciones del movimiento kepleriano (problema de 2 cuerpos)
    - Funciones del problema de N cuerpos
    - Funciones de cálculo de características orbitales """


# Constantes
G = 6.67384E-11    # [m3/(kg s2)] constante gravitacional universal
mu = 3.986044418E5 # [km3/s2] constante gravitacional de la Tierra
RT = 6378.14       # [km] Radio terrestre
J2 = 1.0827E-3     # valor del parámetro de achatamiento terrestre J2


# --------------------------------------------------------------------------------------------------- #


""" Movimiento Kepleriano : Ecuación diferencial 

    Ecuación diferencial matricial del movimiento kepleriano: 
                    d2_r/d_t**2 = - mu*r/(nr**3)
                    
 en el plano:                                  en el espacio:  [  rx  ]     [       vx       ]
              [  rx  ]     [       vx       ]                  [  ry  ]     [       vy       ]
           d  [  ry  ]  =  [       vy       ]               d  [  rz  ]  =  [       vz       ]
           dt [  vx  ]     [ -mu*rx/(nr**3) ]               dt [  vx  ]     [ -mu*rx/(nr**3) ]
              [  vy  ]     [ -mu*ry/(nr**3) ]                  [  vy  ]     [ -mu*ry/(nr**3) ]
                                                               [  vz  ]     [ -mu*rz/(nr**3) ] 
"""

def Kepler(x,t,GM=mu):
    
    if len(x) == 4: # x = [rx,ry,vx,vy]
        
        return array([x[2], x[3], -GM*(x[0])/((x[0]**2 + x[1]**2)**1.5), -GM*(x[1])/((x[0]**2 + x[1]**2)**1.5)])
    
    elif len(x) == 6: # x = [rx,ry,rz,vx,vy,vz]
        
        return array([x[3], x[4], x[5], -GM*(x[0])/((x[0]**2 + x[1]**2 + x[2]**2)**1.5), -GM*(x[1])/((x[0]**2 + x[1]**2 + x[2]**2)**1.5), -GM*(x[2])/((x[0]**2 + x[1]**2 + x[2]**2)**1.5)])


# --------------------------------------------------------------------------------------------------- #


""" Problema de N cuerpos  
      Argumentos: 
                    Nb = Número de cuerpos (mínimo 2)
                    Ub = matriz [Nb x Nc] de estado (componentes de posicion y velocidad) ó
                          vector [Nb*Nc] de estado de los N cuerpos.
                            donde Nc = nº de componentes (en el plano: rx, ry, vx, vy ó en el espacio: rx, ry, rz, vx, vy, vz)
                    
        Output:     Ub1 = matriz [Nb x Nc] de estado (componentes de posicion y velocidad) ó
                           vector [Nb*Nc] de estado de los N cuerpos en tn+1                 
                            (dependiendo del formato de entrada)                             """

def Nbodies(Ub=[],t=[],Nb=2):

    if not Ub.any(): Ub = random_Ub(Nb)

    Nf, Nc = shape(Ub)

    if Nc == 1:  # vector de estado
        Nc = Nf/Nb
        formato = "Vector"
    else:
        formato = "Matriz"

    Us = reshape(Ub,(Nb,2,int(Nc/2)))        # [Nb x 2 x Nc/2] reorganizacion para separar componentes de velocidad de componentes de posicion
    r  = reshape(Us[:,0,:],(Nb,int(Nc/2)))   # [Nb x Nc/2] Posicion : [Nb x 2] si tenemos rx y ry ó [Nb x 3] si tenemos rx, ry y rz
    v  = reshape(Us[:,1,:],(Nb,int(Nc/2)))   # [Nb x Nc/2] Velocidad : [Nb x 2] si tenemos vx y vy ó [Nb x 3] si tenemos vx, vy y vz
    
    Ub1 = zeros((len(Ub),1))
    
    Us1 = reshape(Ub1,(Nb,2,int(Nc/2)))
    drdt = reshape(Us1[:,0,:],(Nb,int(Nc/2)))
    dvdt = reshape(Us1[:,1,:],(Nb,int(Nc/2)))

    for i in range(Nb):

        drdt[i,:] = v[i,:]

        for j in range(Nb):
            
            if j != i:  
                r_d = r[j,:] - r[i,:]
                dvdt[i,:] = dvdt[i,:] + r_d[:]/(norm(r_d)**3)

    if formato == "Vector":
        return Ub1
    elif formato == "Matriz":
        return reshape(Ub1,(Nb,Nc))


""" Inicializador random de N cuerpos 
    A partir de un número de cuerpos (Nb) se determinan posiciones y velocidades random
     dentro de un espacio, 2D o 3D, donde se puede controlar sus dimensiones con kp y
     la diferencia en el órden de magnitud entre posición y velocidad con kv. 
          Por defecto: 3D (dim=3), kp=1 y kv=0.01                                     """

def random_Ub(Nb,dim=3,kp=1,kv=0.01):

    # Inicializacion del vector de estado
    U = zeros((Nb*2*dim,1))
    Us = reshape(U,(Nb,2,dim))

    for i in range(0,Nb):
        phi = uniform(0,2*pi)
        theta = uniform(0,2*pi)
        pos = array([cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)])
        for j in range(dim):
            Us[i,0,j] = pos[j]*kp   # Espacio de -1*k a 1*k
            Us[i,1,j] = uniform(-1,1)*kp*kv # Velocidad log10(1/kv) órdenes de magnitud menor que la posición

    return U

def solarsystem():

    # SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO 
    M = array([1988500., 0.330, 4.87, 5.97, 0.073, 0.642, 1898., 568., 86.8, 102., 0.0130])*1e+24

    U0 = zeros((11*2*3,1))
    U_0 = reshape(U0, (11, 2, 3) ) # punteros
    r0 = reshape(U_0[:, 0, :], (11, 3))
    v0 = reshape(U_0[:, 1, :], (11, 3))

    with open('Initialposition_solarsystem.txt') as data:
        substrings = data.read().split()
        p = [float64(substring) for substring in substrings]

    p = reshape(p, (11, 3))
    r0[:,:] = p[:,:]
    
    with open('Initialvelocity_solarsystem.txt') as data:
        substrings = data.read().split()
        v = [float64(substring) for substring in substrings]

    v = reshape(v, (11, 3))
    v0[:,:] = v[:,:]

    return U0,M


# --------------------------------------------------------------------------------------------------- #

""" Problema restringido de 3 cuerpos """
def bodies3(U,t):
    r = U[0:2]
    drdt = U[2:4]

    r1 = ((r[0] + mu)**2 + r[1]**2)**0.5
    r2 = ((r[0] + mu - 1)**2 + r[1]**2)**0.5

    dvdt1 = -(1-mu)*(r[0] + mu)/(r1**3) - mu*(r[0] + mu - 1)/(r2**3)
    dvdt2 = -(1-mu)*r[1]/(r1**3) - mu*r[1]/(r2**3)

    return array([ drdt[0], drdt[1], r[0] + 2*drdt[1] + dvdt1, r[1] - 2*drdt[0] + dvdt2])

""" Puntos de Lagrange """
def Lagrange(U0,t):

    def f(r):
        x = zeros(4)
        x[0:2] = r
        U = bodies3(x,t)
        return U[2:4]
   
    Lp = zeros([5,2])

    for i in range(5):
        Lp[i,:] = newton(f, U0[i,0:2])

    return Lp

""" Estabilidad de los puntos de Lagrange """

def Stb_Lagrange(U0,t):

    def Jacobian (F, xp):
        
        N = size(xp)
        dx = 1e-10
        
        Jac = zeros([N,N])

        for j in range(N):
            x = zeros(N)
            x[j] = dx
            Jac[:,j] = ( F(xp + x, t) - F(xp - x, t) ) /(2*dx)

        return Jac

    A = Jacobian(bodies3, U0)
    values, vectors = eig(A)

    return values

# --------------------------------------------------------------------------------------------------- #


""" Características orbitales """


def v_circular(r,GM=mu):
    if len(r) > 1:
        nr = 0
        for i in range(0,len(r)):
            nr =+ r[i]**2
        nr = (nr)**0.5
    else:
        nr = r
        
    return (GM/nr)**0.5

def v_orbital(r,a=0,GM=mu):
    if len(r) > 1:
        nr = 0
        for i in range(0,len(r)):
            nr =+ r[i]**2
        nr = (nr)**0.5
    else:
        nr = r
    
    if a == nr or a == 0:
        return v_circular(r,GM)
    else:
        return (2*GM*(1/r + 1/(2*a)))**0.5

