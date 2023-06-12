from numpy import finfo, where, array
#from sympy import diff


# Método de Newton

def Newton(f,xo,err=1e-10,fp=[],maxiter=10):
    
    x = array(xo)
    e = 100 #Error
    i = 0
    if not fp:
        
        dx = finfo(float).eps**0.33  # obtención de un dierencial de x
        x1 = x*(1 + where(f(x) >= 0,dx,-dx))
        
        x = x - f(x)*((x1 - x)/(f(x1) - f(x)))
        
        while e > err or i == maxiter:
            i+=1
            
            x_aux = x
            x = x1 - f(x1)*((x - x1)/(f(x) - f(x1)))
            x1 = x_aux
            #e = max(f(x1)*((x - x1)/(f(x) - f(x1)))) if x is not float else f(x1)*((x - x1)/(f(x) - f(x1)))
    else:
        while e > err:
            i+=1
            
            x = x - f(x)/fp(x)
            e = max(f(x)/fp(x))
    
    return (x, e)

def func(x):
    return x - 4

A, err = Newton(func,1)



""" GUILLE """

from operator import matmul
from numpy import array, zeros, dot, size
from numpy.linalg import inv, norm

## OPERACIONES MATEMÁTICAS ##

 #  
 #
 #  Inputs: 
 #         U : Vector estado en tn
 #         dt: Paso de tiempo 
 #         F : Sistema del que se quiere obtener la matriz Jacobiana o que se quiere solucionar con Newton-Raphson
 #
 #  Return: 
 #
 #         J : Matriz jacobiana 
 #         newton : solución del sistema no linear por Newton-Raphson 

def jacobiano(F,U):

    dim = len(U)
    Dx = 1e-3
    jacobian = array(zeros((dim,dim)))

    for i in range(dim):

        xj = array(zeros(dim))
        xj[i] = Dx
        jacobian[:,i] = (F(U + xj) - F(U - xj))/(2 * Dx)
    
    return jacobian

# LU factorization

def factorization_LU(A):

	n = size(A,1)
	U = zeros([n,n])
	L = zeros([n,n])

	U[0,:] = A[0,:]
	for i in range(0,n):
		L[i,i] = 1

	L[1:n,0] = A[1:n,0]/U[0,0]


	for k in range(1,n):

		for j in range(k,n):
			U[k,j] = A[k,j] - dot(L[k,0:k], U[0:k,j])

		for i in range(k+1,n):
			L[i,k] =(A[i,k] - dot(U[0:k,k], L[i,0:k])) / (U[k,k])

	return [L@U, L, U]

def solve_LU(M,b):

	n=size(b)
	y=zeros(n)
	x=zeros(n)

	[A,L,U] = factorization_LU(M)
	y[0] = b[0]

	for i in range(0,n):
		y[i] = b[i] - dot(A[i,0:i], y[0:i])
		

	x[n-1] = y[n-1]/A[n-1,n-1]

	for i in range(n-2,-1,-1):
		x[i] = (y[i] - dot(A[i, i+1:n+1], x[i+1:n+1])) / A[i,i]
		
	return x

def Inverse(A):

	n = size(A,1)

	B = zeros([n,n])

	for i in range(0,n):
		one = zeros(n)
		one[i] = 1

		B[:,i] = solve_LU(A, one)

	return B


def newton(F, U0):

    dim = len(U0)
    Dx = zeros(dim)
    b  = zeros(dim)
    U1 = U0
    
    eps = 1
    iteration = 0
    itmax = 10000

    while eps > 1e-8 and iteration <= itmax:
   
        J = jacobiano(F,U1)
        b = F(U1)
        Dx = dot(Inverse(J),b)
        U  = U1 - Dx
        eps = norm(U - U1)
        U1  = U
        iteration = iteration + 1

        if iteration == itmax:
            print('Máximo número de iteraciones alcanzado')
    
    return U
