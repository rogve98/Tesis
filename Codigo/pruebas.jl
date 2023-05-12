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

""" Modelo de competencia de 5 especies en competencia. Funciona con tasa de crecimiento
y capacidaddes de carga por separado, es para analizar sistemas aleaotrios completamente.

Como están basadas en la función prueba, me di cuenta que este sistema está mal, por eso 
no llega a la capcidad de carga, las ri/ki no con iguales y por eso afecta en el resultado.
"""

function cincoEspecies(x0,t0,tf,dt,r,K)
    A = zeros(5,5)
    for i in 1:5
        A[i,i] = 1
        for j in 1:5
            if i != j
                A[i,j] = rand()*r[i]/K[i]
            end
        end
    end
    
    function sistema(X)
        sol = zeros(5)
        xs = zeros(5)
        for i in 1:5
            for j in 1:5
                xs[i] += A[i,j]*X[j]
            end
            sol[i] = r[i]*X[i]*(1-xs[i]/K[i])
        end
        return sol               
    end
    
    return RK4(sistema,x0,t0,tf,dt)
end

""" Modelo de competencia de 10 especies en competencia. Mismo caso en en el 
sistema de cinco especies."""

function diezEspecies(x0,t0,tf,dt,r,K)
    A = zeros(10,10)
    for i in 1:10
        A[i,i] = 1
        for j in 1:10
            if i != j
                A[i,j] = rand()*r[i]/K[i]
            end
        end
    end
    
    function sistema(X)
        sol = zeros(10)
        xs = zeros(10)
        for i in 1:10
            for j in 1:10
                xs[i] += A[i,j]*X[j]
            end
            sol[i] = r[i]*X[i]*(1-xs[i]/K[i])
        end
        return sol               
    end
    
    return RK4(sistema,x0,t0,tf,dt)
end

""" Modelo de competencia de 2 especies en competencia (Prueba). 
Nota importante, solo funciona con r = [2,3] y K = [2,3]"""

function pruebas(x0,t0,tf,dt,r,K)
    A = [r[1]/K[1] rand()*r[1]/K[1];rand()*r[2]/K[2] r[2]/K[2]]
    function sistema(X)
        sis = zeros(2)
        xs = zeros(2)
        for i in 1:2
            for j in 1:2
                xs[i] += A[i,j]*X[j]
            end
            sis[i] = r[i]*X[i]*(1-xs[i]/K[i])
        end
        return sis
        #return [r[1]*X[1]*(1-X[1]/K[1]-7.84877*X[2]/K[1]),r[2]*X[2]*(1-X[2]/K[2]-3.52496*X[1]/K[2])]
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
    M = M.*randn(N,N)
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
    #= A , g = randomMatrix(N,p)  
    A = 10 .*A
    for i in 1:N
        A[i,i] = 1
    end =#
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
    
    return (RK4(sistema,x0,t0,tf,dt),A,g)
end

"""Esta función lo que pretende es pintar el sistema resultante de lo que resuelve poblacionesLK.
Para luego encontrar puntos críticos y determinar Jacobiano de cada uno y sus estabilidades.
Lo que nos interesa con esta función es tratar de obtener La matriz de interacciones, sacar
sus eigenvalores y luego cotejarlo con el artículo."""

#Actualización(110523): escribirlo es muy dificil e innecesario. Mejor construimos el Jacobiano manualmente

function sistemaLK(A::Matrix,r::Vector,K::Vector)
    N = length(r)
    STR = string.(zeros(5,5))
    for i in 1:N
        for j in 1:N
            STR[i,j] = "$(A[i,j])*x$j"
            #println("x' = $(r[i])*x$i(1-)")
        end
    end
    return STR
end

"""Esta función solamente es un integrador de sistemas de ecuaciones diferenciales. La 
finalidad del mismo es observar las series de tiempo resultantes de las matrices de interacciones
generadas por la función randomMatrix(args)."""

#Actualización(110523): ESte sistema no funcuiona porque no analiza la serie de tiempo
#del sistema, lo hace de la matriza aleatoria y bueno eso no es precisamente 
#la dinámica del sistema.

function sistemas(x0,t0,tf,dt,A::Matrix)
    N = Int(sqrt(length(A)))
    function sistema(X::Vector)
        # sol = zeros(N)
        # for i in 1:N
        #     for j in 1:N
        #         sol[i] += A[i,j]*X[j]
        #     end
        # end
        return [A[1,1]*X[1]+A[1,2]*X[2],A[2,1]*X[1]+A[2,2]*X[2]]
    end
    RK4(sistema,x0,t0,tf,dt)
end