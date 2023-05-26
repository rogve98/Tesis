using Random
using LinearAlgebra
using Plots
using LaTeXStrings
using Graphs
using GraphPlot
using StatsBase

plotlyjs()

"""RK4

Runge-Kutta 4. Es un integrador para resolver sistemas de ecuaciones diferenciales aunque
probablemente también pueda resolver ecuaciones diferenciales normales.

Parámetros:

f := función de variables Real
x0 := condiciones iniciales del sistema dinámico
t0 := tiempo inicial
tf := tiempo final
h := paso de integración
"""

function RK4(f,x0,t0,tf,h)
    #=al igual que en la función de eulerND, definimos una matriz de dimensión 
    (número de iteraciones × dimensión del sistema dinámico) como conjunto solución=#
    t = range(t0, stop = tf, step = h)
    n = length(t)
    dim = length(x0)
    #lo hacemos en un arreglo de ceros
    xs = zeros(n,dim)
    #imponemos la condición inicial en el primer renglón
    xs[1,:] = x0
    #generamos un ciclo for con las iteraciones de runge-kutta de cuarto orden
    for i in  2:n
        k1 = f(xs[i-1,:])
        k2 = f(xs[i-1,:]+(h/2)*k1)
        k3 = f(xs[i-1,:]+(h/2)*k2)
        k4 = f(xs[i-1,:]+h*k3)
        
        xs[i,:] = xs[i-1,:] + (h/6)*(k1+2*k2+2*k3+k4)
    end
    #=regresamos el resultado en una tupla, con los tiempos en la primera entrada y 
    el conjunto solución en la segunda entrada=#
    return (t , xs)
end

"""
Método de integración de Euler.

f := función a integrar 
x0 := condición inicial
t0, tf := tiempo inicial y final
dt := paso de integración.
"""

function eulerND(f,x0,t0,tf,dt)                #x0 es un vector N-dimensional
    tiempos = range(t0, stop = tf, step = dt)  #Definimos una discretización del paso dt
    n = length(tiempos)                        #Número de iteraciones a realizar
    dim = length(x0)                           #La dimensión del sistema de EDO, en este caso N-dimensional
    
    xs = zeros(n,dim)                          #Arreglo solución del sistema. Es una matriz de 
                                               # n-iteraciones × dimensión del sistema.
    xs[1,:] = x0                               #En el primer renglón de nuestro conjunto solución, imponemos
                                               # las condiciones iniciales
    for i in 2:n #aplicamos las iteraciones
        xs[i,:] = xs[i-1,:] + dt*f(xs[i-1,:])
    end
    return (tiempos,xs)
end  

""" Modelo de competencia de 2 especies en competencia (Prueba). 
Nota importante, solo funciona con r = [2,3] y K = [2,3]"""

function pruebas(x0,t0,tf,dt,r,K)
    A = [r[1]/K[1] rand()*r[1]/K[1];rand()*r[2]/K[2] r[2]/K[2]]
    function sistema(X)
        #Hay que poner manualmente los coeficientes de la matriz aleatoria
        return [r[1]*X[1]*(1-X[1]/K[1]-7.84877*X[2]/K[1]),
        r[2]*X[2]*(1-X[2]/K[2]-3.52496*X[1]/K[2])]
    end
    
    return (RK4(sistema,x0,t0,tf,dt),A)
end

""" Modelo de red aleatoria, la primera función es un generador de enlaces enlacesAleatorios con
base en un número de enlaces N y una probabilidad p """

function enlacesAleatorios(N,p)
    Channel() do channel
        for i in 1:N
            for j in 1:i
                if rand() < p
                    if i == j
                        continue
                    else
                        put!(channel,(i,j))
                    end
                end
            end
        end
    end
end

""" Generamos la red aleatoria con el generador de enlaces aleatorios. """

function redAleatoria(N,p)
    g = Graph(N)
    enlaces = collect(enlacesAleatorios(N,p))
    for i in 1:length(enlaces)
        add_edge!(g,enlaces[i][1],enlaces[i][2])
    end
    return g
end

"""Esta función nos genera matrices aleatorias con base en el modelo de red aleatoria.
Lo hice de esta manera para que pudieramos tener una representación gráfica de las 
interacciones y al mimsmo tiempo obtener matrices aleatorias con valores aleatorios además."""

function randomMatrix(N,p)
    g = redAleatoria(N,p)
    M = adjacency_matrix(g)
    for i in 1:N 
        M[i,i] = 1
    end
    M = 10*M.*randn(N,N)
    for i in 1:N
        M[i,i] = 1
    end
    return (Matrix(M), g)
end

"""Voy a diseñar una función que sea compatible con la matiz que regresa randomMatrix(N,p), 
es decir que solo necesitará la matriz para trabajar sin tener en cuenta las capacidades de carga
ni las tasas de crecimiento. Vamos a ve que sale. La variable params contiene la dimensión, la 
probabilidad de contectancia, la tasa de crecimiento y la capacidad de carga

params := [N,p,r,K]

"""

function poblacionesLK(x0,t0,tf,dt,params)
    N = params[1]
    p = params[2]
    r = params[3]
    K = params[4]
    #A , g = randomMatrix(N,p)   
    A = [1.0        0.0      13.5989   -3.28364  0.0;
    0.0        1.0       0.0       6.10228  0.0;
    0.574493   0.0       1.0       2.74343  0.0;
    2.59557   -3.14685  -2.20031   1.0      0.0;
    0.0        0.0       0.0       0.0      1.0;
   ]
   g=2 
    function sistema(X::Vector)
        sis = zeros(N)
        xs = zeros(N)
        for i in 1:N
            for j in 1:N
                xs[i] += A[i,j]*X[j]
            end
            sis[i] = r[i]*X[i]*(1-xs[i]/K[i])
        end
        return sis
    end
    
    return (RK4(sistema,x0,t0,tf,dt),eulerND(sistema,x0,t0,tf,dt),A,g)
end

"""
Jacobiano n-dimensional

X := especies
P := parámetros, [r,K,N,A]
"""

function Jacobiano(X::Vector,P::Array)
    r = P[1]
    K = P[2]
    N = P[3]
    A = P[4]
    M = zeros(N,N)
    for i in 1:N
        for j in 1:N
            if i == j
                xs = zeros(N)
                for k in 1:N
                    xs[i] += A[i,k]*X[k]
                end
                M[i,i] = r[i]*(1-xs[i]/K[i])-r[i]*X[i]/K[i]
            else
                M[i,j] = -r[i]*X[i]*A[i,j]/K[i]
            end
        end
    end
    return M
end        

""" Sistema del LK

X := Especies
P := parámetros, [r,K,N,A]

"""

function sistema(X::Vector,P::Array)
    r = P[1]
    K = P[2]
    N = P[3]
    A = P[4]
    sis = zeros(N)
    xs = zeros(N)
    for i in 1:N
        for j in 1:N
            xs[i] += A[i,j]*X[j]
        end
        sis[i] = r[i]*X[i]*(1-xs[i]/K[i])
    end
    return sis
end

"""
Newton Rhapson multidimensional.

Jacobiano := función que nos regresa el jacobiano del sistema.
x0 := condición inicial.
P := parámetros, [r,K,N,A]
n := número de pasos para iterar Newton Rhapson
"""

function nrMulti(Jacobiano::Function,x0::Vector,P::Array,n::Int)
    h = 0.1
    N = P[3]
    sol = zeros(n,N)
    sol[1,:] = x0
    for i in 2:n
        sol[i,:] = sol[i-1,:] - inv(Jacobiano(sol[i-1,:],P))*sistema(sol[i-1,:],P)
    end
    return sol[end,:]
end