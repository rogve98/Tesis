"""
Esta función es la misma que se ubica en Sistema.jl solamente que se encuentra impementada con Threads.
Vamos a probar el cómputo en paralelo para observar si se optimizan los tiempos de ejecución.
"""

function transicionParalel(params::Parametros, p)
    sol = Vector{Int}(undef, length(p))  # Preasignamos el tamaño del vector para evitar accesos concurrentes
    medidas::Int = 1500
    Threads.@threads for idx in 1:length(p)
        i = p[idx]
        estables = []
        params_local = deepcopy(params)  # Asegúrate de copiar los parámetros para evitar conflictos entre hilos
        params_local.p = i
        for j in 1:medidas
            xs = integrador(params_local)
            push!(estables, xs)
        end

        sol[idx] = count(x -> x == 1, esEstable.(estables))
    end

    return sol
end


"""
Esta función es la misma que se ubica en Sistema.jl solamente que se encuentra implementada con Threads. 
Vamos a probar cómputo en paralelo para revisar si se pueden optimizar tiempos de ejecución.
"""

function transicionσParalel(params::Parametros,σ)
    sol = Vector{Int}(undef,length(σ))
    medidas::Int = 1500
    Threads.@threads for idx in 1:length(σ)
        i = σ[idx]
        estables = []
        params_local = deepcopy(params)
        params_local.σ = i
        for j in 1:medidas 
            xs = integrador(params_local)
            push!(estables,xs)
        end
        sol[idx] = count(x -> x == 1, esEstable.(estables))
    end
    return sol    
end

"""
Código alternativo para utilizar Threads con el fin de minimizar tiempos de ejecución usando varios
cores en la ejecución. La ejecución para σ = 0.5 tardó 9 días! Entonces para sigmas superiores 
esto se va a volver un relajo. Se busca minimzar tiempos de ejecución
"""

function distDiagonalParalel(params::Parametros,p)
    medidas = 10
    σ = params.σ
    for i in p
        #diagonales = []
        jacobianos = []
        params.p = i 

        # Crear el lock fuera del ciclo
        lck = ReentrantLock()

        # Paralelización en la generación de diagonales
        Threads.@threads for i in 1:medidas
            _,M = DiagonalJ(params::Parametros)

            # Bloque crítico con el lock
            lock(lck)  # Bloqueamos antes de modificar el arreglo compartido
            try
               #push!(diagonales, d)
                push!(jacobianos, M)
            finally
                unlock(lck)  # Desbloqueamos después de la operación
            end
        end

        # Escritura secuencial de los resultados a un archivo CSV
        #writedlm("Diagonales_s$σ.p$i.csv", diagonales)
        writedlm("Jacobianos_s$σ.p$i.csv", jacobianos)
    end
end
    
    
    
   #=  medidas = 100
    σ = params.σ
    for i in p
        diagonales = []
        params.p = i 

        Threads.@threads for i in 1:medidas
            d = DiagonalJ(params::Parametros)
            push!(diagonales,d)
        end

        writedlm("Diagonales_s$σ.p$i.csv",diagonales)
    end
end =#