"""Voy a diseñar una función que sea compatible con la matiz que regresa randomMatrix(N,p), 
es decir que solo necesitará la matriz para trabajar sin tener en cuenta las capacidades de carga
ni las tasas de crecimiento. Vamos a ve que sale. La variable params contiene la dimensión, la 
probabilidad de contectancia, la tasa de crecimiento y la capacidad de carga

params es el objeto que tiene todo lo necesario para esta función.

"""

function poblacionesLK(params::Parametros)
    N = params.N
    p = params.p
    r = params.r
    K = params.K
    σ = params.σ
    red = params.S
    A , g = randomMatrix(N,p,σ,red)   
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
    
    return Soluciones(RK4(sistema,params.x0,params.t0,params.tf,params.h)[1],
    RK4(sistema,params.x0,params.t0,params.tf,params.h)[2],
    eulerND(sistema,params.x0,params.t0,params.tf,params.h)[2],
    A,
    g
    )
end

"""
Función para determinar que el valor introducido es un real, o como tal no es un NaN (Not a Number)
"""

function no_es_NaN(valor)
    !isnan(valor)
end

"""
Esta función es capaz de mapear todos los valores de una matriz en donde determina si la solución del sistema
es estable o no. El mapeo determinar si existen NaN o no.
"""

function esEstable(sol::Matrix)
    all(no_es_NaN,sol)
end

"""
Esta función es igual a poblacionesLK solamente que el esfuerzo computacional se centra
en la generación de las matrices para poder reducir el tiempo de ejecución del conjunto de 
transición. Opera con el mismo objeto Parametros
"""

function integrador(params::Parametros)
    N = params.N
    p = params.p
    r = params.r
    K = params.K
    σ = params.σ
    red = params.S
    A , _ = randomMatrix(N,p,σ,red)
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
    return RK4(sistema,params.x0,params.t0,params.tf,params.h)[2]
end

"""
Esta función genera la gráfica de la transición. Opera sobre un conjunto de probabilidades
(en general de longitud 50)
"""

function transicion(params::Parametros,p)
    sol = []
    medidas = 200
    for i in p
        estables = []
        params.p = i
        for j in 1:medidas
            xs = integrador(params::Parametros)
            push!(estables, xs)
        end
        push!(sol,count(x -> x == 1,esEstable.(estables)))
    end
    return sol
end

"""

"""

function LKBarabasi(params::Parametros,m)
    N = params.N
    r = params.r
    K = params.K
    σ = params.σ
    A , g = BMatrix(N,m,σ)
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
    return Soluciones(RK4(sistema,params.x0,params.t0,params.tf,params.h)[1],
    RK4(sistema,params.x0,params.t0,params.tf,params.h)[2],
    eulerND(sistema,params.x0,params.t0,params.tf,params.h)[2],
    A,
    g
    )
end

"""

"""

function integradorB(params::Parametros,m)
    N = params.N
    r = params.r
    K = params.K
    σ = params.σ
    A , _ = BMatrix(N,m,σ)
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
    return RK4(sistema,params.x0,params.t0,params.tf,params.h)[2]
end

"""

"""

function transicionB(params::Parametros,m)
    sol = []
    medidas = 200
    for i in m
        estables = []
        for j in 1:medidas
            xs = integradorB(params::Parametros,i)
            push!(estables, xs)
        end
        push!(sol,count(x -> x == 1,esEstable.(estables)))
    end
    return sol
end