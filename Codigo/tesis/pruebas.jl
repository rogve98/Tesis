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

"""Sistema de 5 x 5 que podría ser útil"""
#     A = [1.0        0.0      13.5989   -3.28364  0.0;
#     0.0        1.0       0.0       6.10228  0.0;
#     0.574493   0.0       1.0       2.74343  0.0;
#     2.59557   -3.14685  -2.20031   1.0      0.0;
#     0.0        0.0       0.0       0.0      1.0;
#    ]
#    g=2 

function pruebaDiagCeros(params::Parametros)
    r = params.r
    K = params.K
    function sistema(X)
        return [r[1]X[1]-(r[1]*randn()*X[1]*X[2])/K[1],r[2]X[2]-(r[2]*randn()*X[1]*X[2])/K[2]]
    end
    return RK4(sistema,params.x0,params.t0,params.tf,params.h)
end

"""
Se guarda aquí la función del jacobiano porque es una versión anterior a la que se encuentra en Estabilidad.jl
"""

function jacobianoSimple(A::Matrix,X::Vector,params::Parametros)
    r = params.r
    K = params.K
    N = params.N
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

"""
Método de gauss newton para ajuste de curva sigmoidal, de momento no funciona xdxd
"""

function gauss_newton(x, y, β_init; tol=1e-6, max_iter=100)
    # Definir la función modelo (sigmoide)
    function modelo(x, β)
        return 1 ./ (1 .+ exp.(-β[1] .* (x .- β[2])))
    end

    # Definir la matriz Jacobiana
    function jacobiano(x, β)
        J = zeros(length(x), 2)
        f_x = modelo(x, β)
        J[:, 1] = (x .- β[2]) .* f_x .* (1 .- f_x)  # Derivada respecto a β1
        J[:, 2] = -β[1] .* f_x .* (1 .- f_x)        # Derivada respecto a β2
        return J
    end

    β = copy(β_init)
    for iter in 1:max_iter
        # Residuos
        r = y .- modelo(x, β)
        # Jacobiano
        J = jacobiano(x, β)
        # Actualización de parámetros
        Δβ = -inv(J' * J) * J' * r
        β += Δβ
        # Verificar la convergencia
        if norm(Δβ) < tol
            println("Convergencia alcanzada en la iteración $iter")
            break
        end
    end

    return β
end

# # Ejemplo de uso con datos simulados
# x_data = collect(-10:0.5:10)
# y_data = 1 ./ (1 .+ exp.(-1.5 .* (x_data .- 2))) + 0.1*randn(length(x_data))  # Agregar algo de ruido
# β_inicial = [1.0, 0.0]  # Parámetros iniciales para la sigmoide (β1, β2)

# # Llamar a la función de Gauss-Newton
# β_opt = gauss_newton(x_data, y_data, β_inicial)
# println("Parámetros ajustados: β1 = $(β_opt[1]), β2 = $(β_opt[2])")


"""Ejemplo para función multivariada para el caso de la ley elíptica de Allesina"""
# μ = [0, 0]
# Σ = [1  0.9;
#      0.9 10]
# p = MvNormal(μ, Σ)

# X = range(-8, 8, length=100)
# Y = range(-8, 8, length=100)
# Z = [pdf(p, [x,y]) for y in Y, x in X] # Note x-y "for" ordering
# contourf(X, Y, Z, color=:viridis)
