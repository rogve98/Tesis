# Guía de ejecución:
# La función Interacciones se ejecuta primero para poder obtener el objeto
# estabilidad necesario para la invocación del resto de funciones de este módulo
# Luego pasamos a calcular el Jacobiano con el objeto estabilidad y eso nos
# estaría generando la Matriz de interacciones ya evaluada en el atractor
# que es un vector n-dimensional. 

# Esto sería suficiente, el resto de funciones: Sistema y nrMulti ya no son 
# tan indispensables puesto que el objeto estabilidad contiene toda la información
# necesaria para determinar el Jacobiano: Matriz de interacciones.

# Lo último nada más sería determinar los eigenvalores del sistema y corroborar
# que tienen parte real negativa.

"""
Esta función solo reúne la información necesaria para poder utilizar el objeto 
estabilidad para podes usarlo como argumento de las funciones de este módulo.
El objeto estabilidad contiene la información completa de los parámetros y de las
soluciones del sistema y agrega el punto de equilibrio X para poder ser utilizado.
"""

function Interacciones(params::Parametros,sol::Soluciones)
    X = sol.rk4[end,:]
    return estabilidad(params,sol,X)
end

"""
Jacobiano n-dimensional. Se deteminó de manera analítica. Utiliza el objeto
estabilidad como argumento para sus cálculos.

estailidad(params::Parametros,sol::Soluciones,X::Vector)
"""

function Jacobiano(E::estabilidad)
    r = E.params.r
    K = E.params.K
    N = E.params.N
    A = E.sol.A
    X = E.X
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
Se implementa la siguiente función en donde se integra al sistema con RK4, se extrae el punto fijo y se
valida que el sistema es estable para poder sacar el Jacobiano del sistema. Luego se aplican las derivadas
parciales correspondientes para tener obtener finalmente la matriz de interacciones.
"""

function DiagonalJ(params::Parametros)
    N = params.N
    r = params.r
    K = params.K
    A = incidencias(params::Parametros)
    function Estables(A::Matrix)
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
        return LK = RK4(sistema,params.x0,params.t0,params.tf,params.h)[2]
    end
    LK = Estables(A)
    while !(esEstable(LK))
        A = incidencias(params)
        LK = Estables(A)        
    end
    X = LK[end,:]
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
    return diag(M)
end

"""
Esta función nos servirá para explorar la distribución de las diagonales de los
jacobianos, esto con el fin de poder hallar un parámetro de criticalidad en los
diagramas de transición de fase... Hay algunas ideas en el texto de mercedes que
encuadran con esto pero al día de hoy 11-Sep-24 no estoy seguro de la relación.
"""

function distDiagonal(params::Parametros,p)
    medidas = 5
    σ = params.σ
    for i in p
        diagonales = []
        #jacobianos = []
        params.p = i 
        for i in 1:medidas
            d = DiagonalJ(params::Parametros)
        #    push!(jacobianos,J)
            push!(diagonales,d)
        end
        writedlm("Diagonales_s$σ.p$i.csv",diagonales)
        #writedlm("Jacobianos p$i.csv",jacobianos)
    end
end

""" Sistema del LK. Únicamente es la reproducción del sistema n-dimensional. Utiliza la
información del objeto estabilidad para poder realizar sus cálculos.

estailidad(params::Parametros,sol::Soluciones,X::Vector)
"""

function sistema(E::estabilidad)
    r = E.params.r
    K = E.params.K
    N = E.params.N
    A = E.M.A
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
Newton Rhapson multidimensional. Se encarga de buscar raíces, sin embargo esta 
función ya no se estará ocupando ya que establidad.X tiene el vector solución
que se espera para poder ejecutar las funciones anterioes. Sin embargo lo dejaré
por aquí nada más pal registro/recuerdo jeje.

Jacobiano := función que nos regresa el jacobiano del sistema.
x0 := condición inicial; este fue remplazado por E.X
P := parámetros, [r,K,N,A]; este por el objeto E completo
n := número de pasos para iterar Newton Rhapson; con 100 o 1000 iteraciones esta chido.
"""

function nrMulti(Jacobiano::Function,E::estabilidad,n::Int)
    N = E.params.N
    sol = zeros(n,N)
    sol[1,:] = E.X
    for i in 2:n
        sol[i,:] = sol[i-1,:] - inv(Jacobiano(E))*sistema(E)
    end
    return sol[end,:]
end