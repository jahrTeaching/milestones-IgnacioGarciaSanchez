from numpy import array

""" Movimiento Kepleriano """

""" Ecuación diferencial matricial del movimiento kepleriano en el plano: 
                    d2_r/d_t**2 = - mu*r/(nr**3)
                    
                  [  rx  ]     [       vx       ]
                d [  ry  ]  =  [       vy       ]
               dt [  vx  ]     [ -mu*rx/(nr**3) ]
                  [  vy  ]     [ -mu*ry/(nr**3) ]                           """


def Kepler(x,t,mu=3.986044418E5):
    
    if len(x) == 4:
        
        return Kepler2D(x,t,mu)
    
    elif len(x) == 6:
        
        return Kepler3D(x,t,mu)


def Kepler2D(x,t,mu=3.986044418E5):
    
    return array([x[2], x[3], -mu*(x[0])/((x[0]**2 + x[1]**2)**1.5), -mu*(x[1])/((x[0]**2 + x[1]**2)**1.5)])
         

def Kepler3D(x,t,mu=3.986044418E5):
    
    return array([x[3], x[4], x[5], -mu*(x[0])/((x[0]**2 + x[1]**2)**1.5), -mu*(x[1])/((x[0]**2 + x[1]**2)**1.5), 0])
