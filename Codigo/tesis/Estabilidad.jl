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