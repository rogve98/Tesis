using LinearAlgebra
using Plots
using LaTeXStrings
using Graphs
using GraphPlot
using Distributions
using CSV
using DataFrames

plotlyjs()



include("Integradores.jl") #Contiene RK4 y eulerND
include("Redes.jl") #Contiene todo el material de Redes
include("Estabilidad.jl") #Contiene las funciones de estabilidad por punto critico
include("structs.jl") #Contiene los structs
include("Sistema.jl") #Contiene las funciones del sistema y la transición.

"""
Redes aleatorias que fungen como Jacobianos (matrices de interacciones). Es el estilo 
de May en donde se salta todo el sistema para solo trabajar con matrices aleatorias que 
dependen de una probabilidad y con la diagonal con -d, representando la capacidad de 
carga del sistema.
"""

function interacciones(N,p,σ,d)
    g = redAleatoria(N,p)
    M = adjacency_matrix(g)
    dist = Normal(0,σ)
    M = M.*rand(dist,N,N)
    for i in 1:N
        M[i,i] = -d
    end
    M = 1/sqrt(N*p)*M
    return (Matrix(M), g)
end

"""
Función para generar un círculo con centro en (h,k) y radio radio r.
"""

function circulo(h,k,r)
    θ = range(0, stop = 2π, length = 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end