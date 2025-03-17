using LinearAlgebra
using Plots
using LaTeXStrings
using Graphs
using GraphPlot
using Distributions
using CSV
using DataFrames
using StatsBase
using DelimitedFiles

plotlyjs()


include("structs.jl") #Contiene los structs
include("Integradores.jl") #Contiene RK4 y eulerND
include("Redes.jl") #Contiene todo el material de Redes
include("Estabilidad.jl") #Contiene las funciones de estabilidad por punto critico
include("Sistema.jl") #Contiene las funciones del sistema y la transición.
include("pruebas.jl") #Contiene código reciclado y algunas pruebas importantes por hacer
include("Scripts.jl") #Contiene scripts para poder guardar los datos en archivos .CSV
include("paralel.jl") #Vamos a implementar cómputo en paralelo.
include("Estadistica.jl")

"""
Redes aleatorias que fungen como Jacobianos (matrices de interacciones). Es el estilo 
de May en donde se salta todo el sistema para solo trabajar con matrices aleatorias que 
dependen de una probabilidad y con la diagonal con -d, representando la capacidad de 
carga del sistema.
"""

function interacciones(N,p,σ,red)
    g = redAleatoria(N,p,red)
    M = adjacency_matrix(g)
    dist = Normal(0,σ)
    Id = -1 * Matrix(I,N,N)
    M = M.*rand(dist,N,N) + Id
    #M = 1/sqrt(N*p)*M
    return Matrix(M)#, g
end

"""
Interacciones de la red de Allesina, este es un prototipo para este momento. Hay que ver
como construir adecuadamente la distribución bivariada y/o verificar si es correcta la 
forma en que propuse la construcción de la matriz de interacciones.

En este caso σ contiene en la primer columna los valores de los centros (las μ's)
y la segunda columna tiene las desviaciones estandar σ's.
"""

function interaccionesEliptica(N,p,σ,red)
    g = redAleatoria(N,p,red)
    M = adjacency_matrix(g)
    dist1 = Normal(σ[1,1],σ[1,2])
    dist2 = Normal(σ[2,1],σ[2,2])
    Id = -3 * Matrix(I,N,N)
    M = M.*UpperTriangular(rand(dist1,N,N)) + M.*LowerTriangular(rand(dist2,N,N)) + Id
    return Matrix(M)
end

"""
Función para generar un círculo con centro en (h,k) y radio radio r.
"""

function circulo(h,k,r)
    θ = range(0, stop = 2π, length = 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

"""
Función para generar una eli´se con centro  (h,k) y con semiejes a y b
"""

function elipse(h, k, a, b)
    θ = range(0, stop = 2π, length = 500)
    h .+ a .* cos.(θ), k .+ b .* sin.(θ)
end

"""
Esta transición de fase lo que busca es ver como varía la estabilidad en función de p 
para poder comparar los resultados de May con los míos, dependiendo del tiempo que tarde 
la compilación veré para que Ns considero estos cálculos.
"""

function transicionMay(N,p,σ,red)
    sol = []
    medidas = 1000
    for i in p
        estables = []
        for _ in 1:medidas
            M = interacciones(N,i,σ,red)
            eigs = eigvals(M)
            push!(estables, eigs)
        end
        neg = all.(x -> real(x) < 0, estables)
        push!(sol,count(x -> x == 1,neg))
    end
    return sol
end

"""
Función para determinar la transición de May en función de sigma. Esta función me parece que es 
todavía más relevante que la anerior porque el parámetro crítico que propone may se ajusta mejor
a los resultados en función de sigma.
"""

function transicionσMay(N,p,σ,red)
    sol = []
    medidas = 1000
    for i in σ
        estables = []
        for _ in 1:medidas
            M = interacciones(N,p,i,red)
            eigs = eigvals(M)
            push!(estables, eigs)
        end
        neg = all.(x -> real(x) < 0, estables)
        push!(sol,count(x -> x == 1,neg))
    end
    return sol
end


"""
Las siguientes tres funciones servirán como alternativa para sacar los gráficos de transición
de may para redes dirigidias (completamente aleatorias). La motivación de estas funciones es la de 
optimizar tiempos de ejecución ya que mediante el método tradicional utilizando red = "dirigida"
se tarda alrededor de un día en compilar, y justo ahora necesitamos ahorrarnos el mayor tiempo
posible.
"""

function interaccionesMayDir(N,p,σ)
    M = zeros(N,N)
    dist = Normal(0,σ)
    for i in 1:N
        for j in 1:N
            if i == j 
                continue
            end
            if rand() < p
                M[i,j] = rand(dist)
            end
        end
    end
    Id = -1 * Matrix(I,N,N)
    return M+Id
end

function transicionMayDir(N,p,σ)
    sol = []
    medidas = 1000
    for i in p
        estables = []
        for _ in 1:medidas
            M = interaccionesMayDir(N,i,σ)
            eigs = eigvals(M)
            push!(estables, eigs)
        end
        neg = all.(x -> real(x) < 0, estables)
        push!(sol,count(x -> x == 1,neg))
    end
    return sol
end

function transicionσMayDir(N,p,σ)
    sol = []
    medidas = 1000
    for i in σ
        estables = []
        for _ in 1:medidas
            M = interaccionesMayDir(N,p,i)
            eigs = eigvals(M)
            push!(estables, eigs)
        end
        neg = all.(x -> real(x) < 0, estables)
        push!(sol,count(x -> x == 1,neg))
    end
    return sol
end