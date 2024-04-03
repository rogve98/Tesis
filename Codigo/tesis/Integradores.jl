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
    t = range(t0, stop = tf, step = h)
    n = length(t)
    dim = length(x0)
    xs = zeros(n,dim)
    xs[1,:] = x0
    for i in  2:n
        k1 = f(xs[i-1,:])
        k2 = f(xs[i-1,:]+(h/2)*k1)
        k3 = f(xs[i-1,:]+(h/2)*k2)
        k4 = f(xs[i-1,:]+h*k3)
        
        xs[i,:] = xs[i-1,:] + (h/6)*(k1+2*k2+2*k3+k4)
    end
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