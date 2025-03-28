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
    red = params.Red
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
    r = params.r
    K = params.K
    A = incidencias(params)
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
    X = RK4(sistema,params.x0,params.t0,params.tf,params.h)[2][end,:]
    if any(isnan,X)
        return 0
    else
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
        evs = eigvals(M)
        return Int(all(x->x<0,real(evs)))
    end
end

"""
Esta función genera la gráfica de la transición. Opera sobre un conjunto de probabilidades
(en general de longitud 50 y 100)
"""

function transicion(params::Parametros,p)
    sol = []
    medidas = 10
    for i in p
        estables = []
        params.p = i
        for j in 1:medidas
            xs = integrador(params::Parametros)           
            push!(estables, xs)
        end
        push!(sol,count(x -> x == 1,estables))
    end
    return sol
end

"""
Esta función nos servirá para poder generar gráficos de transición para un conjunto de 100 sigmas.
El algoritmo es el mismo que en la función anterior solo que en lugar de evaluar sobre las probabilidades
se centra en las sigmas de la distribución normal.
"""

function transicionσ(params::Parametros,σ)
    sol = []
    medidas = 10
    for i in σ
        estables = []
        params.σ = i 
        for j in 1:medidas 
            xs = integrador(params::Parametros)
            push!(estables,xs)
        end
        push!(sol,count(x -> x  == 1, estables))
    end
    return sol    
end

"""

"""
function transicionN(params::Parametros,N)
    sol = []
    medidas = 1000
    for i in N
        estables = [] 
        params.N = Int(i) 
        params.x0 = ones(Int(i))
        params.r = 2*ones(Int(i))
        params.K = 5*ones(Int(i))
        for j in 1:medidas 
            xs = integrador(params::Parametros)
            push!(estables,xs)
        end
        push!(sol,count(x -> x == 1,esEstable.(estables)))
    end
    return sol    
end

"""
Sistema de competencia de especies (completo) para la red de barabasi albert.
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
Integrador para la red de barabasi albert.
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
Para realizar los gráficos de transición de fase para la red de albert barabasi.
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

"""
Segunda versión de la función intgeración para recuperar todos los Jacobianos que resultaron
estables de cada simulación exitosa.
"""

function integradorJacobs(params::Parametros)
    N = params.N
    r = params.r
    K = params.K
    A = incidencias(params)
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
    X = RK4(sistema,params.x0,params.t0,params.tf,params.h)[2][end,:]
    if any(isnan,X)
        return 0,nothing,A # A es la matriz de incidencias
    else
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
        evs = eigvals(M)
        if  all(x->x<0,real(evs))
            return 1,M,A
        else
            return 0,nothing,A
        end
    end
end