using Random
using LinearAlgebra
using Plots
using LaTeXStrings
using Graphs
using GraphPlot
using StatsBase
using Distributions

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
    return [tiempos,xs]
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
    d = Normal(0,0.1)
    g = redAleatoria(N,p)
    M = adjacency_matrix(g)
    M = M.*rand(d,N,N)
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
    A , g = randomMatrix(N,p)   
#     A = [1.0        0.0      13.5989   -3.28364  0.0;
#     0.0        1.0       0.0       6.10228  0.0;
#     0.574493   0.0       1.0       2.74343  0.0;
#     2.59557   -3.14685  -2.20031   1.0      0.0;
#     0.0        0.0       0.0       0.0      1.0;
#    ]
#    g=2 
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

"""Red Circular"""

function Circular(N)
    g = Graph(N)
    for i in 1:N
        add_edge!(g,i,i+1)
    end
    add_edge!(g,N,1)
    return g#gplot(g,nodelabel=1:N)
end

"""
n:=número de nodos totales
m:=número de enlaces nuevos por nodo; con este parámetro definimos la red inicial a partir de una red GNL

La idea del algorítmo es construir una red inicial aleatoria tipo GNL para asegurarnos que sea mayormente conexa
(la GNP no lo asegura al 100% a menos de que ocupemos la teoría de punto crítico pero eso lo veremos después). 
Se verifica si es conexa mediante un condicional y si si generamos el algorítmo de albert barabasi que consiste en
definit probabilidades de enlace en función del grado de cada nodo j dividido entre el grado total (suma de grados).
De aquí definimos un vector de m entradas que serán los nuevos nodos a los que irán conectados un último y nuevo 
nodo.

targets funciona de tal manera que hace una selección desde el nodo 1 hasta el nodo nv(G) y por medio de las
probabilidades se van eligiendo esos nodos. Por ejemplo para la primera entrada del vector será más probable que
se escoga el número entre 1 y nv(G) que tenga la probabilidad más alta en probabilities. Por ejemplo si el vector
de probabilities en su tercera entrada tiene un grado alto (probabilidad alta), es muy probable que el nodo que 
vaya a salir seleccionado entre 1 y nv(G) es el 3 (porque le corresponde la tercera entrada). De aquí actúa el 
replace y el 3 ya no puede ser escogido en futuras ocasiones pero el procedimiento se repite para la siguiente 
entrada del vector targets, y así hasta acabar con todos.

Posteriormente se agrega un nodo extra para poder conectar los targets a este nodo extra y la manera de conectarlos
es por medio del índice i del ciclo for que justamente coincide con el nodo agregado, de aquí nada más falta
conectar ese nodo extra con todos los nodos de target. Y bueno este procedimiento se realiza iteradamente hasta
terminar con los valores de n.
"""

function barabasi_albert(n, m)
    # Inicializar un grafo completo con m nodos
    G = Circular(m) #GNL(m,m)
    if is_connected(G)
        # Implementar el modelo de Barabási-Albert
        for i in m+1:n
            # Calcular las probabilidades de conexión para nodos existentes
            probabilities = [degree(G, j) / sum(degree(G, k) for k in vertices(G)) for j in vertices(G)]

            # Elegir m nodos existentes basados en las probabilidades
            targets = sample(1: nv(G),Weights(probabilities),m,replace=false)

            # Agregar un nuevo nodo
            add_vertex!(G)

            # Conectar el nuevo nodo a los nodos seleccionados
            for target in targets
                add_edge!(G,i,target)
            end
        end
    end
    M = Matrix(adjacency_matrix(G))
    M = 10*M.*randn(N,N)
    for i in 1:N
        M[i,i] = 1
    end

    return G,M
end


"""
En la semana del 4 de marzo, 2024 se hizo un experimento tratando de testear la calidad de poblacionesLK, se creó
este sistema para poder visualizar una comparativa entre ambas formas de sintetizar el algorítmo. Los resultados 
encontrados en su momento fue que que este sistema conserva desde las primeras iteraciones una precisión de 16 decimales,
mientras que poblacionesLK presenta una baja en precisión comenzando con 4 decimales y de ahí va desarrollando hasta
acoplarse a este sistema. Estuve revisando porque podría darse esa falta de precisión e intenté fijar los vectores Y
valores a BigFloat pero realmente no tuvo un gran efecto. Si es necesario reproducir nuevas pruebas dejemos esta
evidencia y construcción para retomar si se requiere.
"""

function prueba(x0,t0,tf,dt,params)
    N = params[1]
    p = params[2]
    r = params[3]
    K = params[4]
    function sistema(X::Vector) #[2X[1]*(1-(X[1]+X[2])/2),2X[2]*(1-(X[2]+2X[1])/3)]
        return [
            2X[1]*(1-(1.0X[1] +       0.0X[2] +      13.5989X[3]   -3.28364X[4] +  0.0X[5])/2),
            2X[2]*(1-(0.0X[1] +       1.0X[2] +      0.0X[3] +      6.10228X[4] +  0.0X[5])/3),
            2X[3]*(1-(0.574493X[1] +  0.0X[2] +      1.0X[3] +      2.74343X[4] +  0.0X[5])/2),
            2X[4]*(1-(2.59557X[1]    -3.14685X[2]   -2.20031X[3] +  1.0X[4] +      0.0X[5])/3),
            2X[5]*(1-(0.0X[1] +       0.0X[2] +      0.0X[3] +      0.0X[4] +      1.0X[5])/2)
        ]
    end
    return RK4(sistema,x0,t0,tf,dt)
end

function no_es_NaN(valor)
    !isnan(valor)
end

function esEstable(sol::Matrix)
    all(no_es_NaN,sol)
end

function contadorEstables(sol::Matrix)
    contador = []
    N = length(sol[1,:])
    for i in 1:N
        if in(1,isnan.(sol[:,i]))
            push!(contador,i)
        end
    end
    return length(contador)
end