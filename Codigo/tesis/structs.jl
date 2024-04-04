"""
Definimos un objeto mutable que guarde todos los parámetros necesarios para poder realizar
las operaciones de las fuciones participantes en este código, en general definimos:
N := Número de especies
σ := Desviación estándar de la distribución normal
p := Probabilidad de la red aleatoria; estrictamente debe ser un parámetro en [0,1]
x0 := Vector de las condiciones iniciales del sistema
t0 := Tiempo inicial
tf := Tiempo final
h := Paso de integración
r := Vector N-dimensional que contiene las tasas de crecimiento de las especies
K := Vector N-dimensional que contiene las capacidades de carga de las poblaciones
"""

mutable struct Parametros
    N::Int          
    σ::Float64
    p::Float64
    x0::Vector
    t0::Int
    tf::Int
    h::Float64
    r::Vector
    K::Vector 
end

"""
Este objeto contendrá las soluciones requeridas del sistema, solo se enfoca en la integración
del sistema como tal. Quizás más adelante sea necesario definir más objetos para otras
cantidades/soluciones. De momento este objeto deine:
t := tiempos de integración, es un range()
rk4 := solución integrada por RK4
Euler := solución integrada por Euler
A := Refiere a la matriz de interacciones, los α_{ij}
g := Es la red del sistema
"""

mutable struct Soluciones
	t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    rk4::Matrix
    Euler::Matrix
    A::Matrix
    g::SimpleGraph{Int64}
end

"""
Objeto neceario para las funciones del módulo de Estabilidad, contiene toda la información e 
incluso contiene el vector atractor del sistema en donde se determina que las especies coexisten.
"""

mutable struct estabilidad
    params::Parametros
    M::Soluciones
    X::Vector
end