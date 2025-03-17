"""
Esta función es la misma que se ubica en Sistema.jl solamente que se encuentra impementada con Threads.
Vamos a probar el cómputo en paralelo para observar si se optimizan los tiempos de ejecución.
"""

function transicionParalel(params::Parametros, p)
    sol = Vector{Int}(undef, length(p))  # Vector preasignado
    medidas = 1000

    Threads.@threads for idx in eachindex(p)
        params_local = deepcopy(params)  # Evita conflictos entre hilos
        params_local.p = p[idx]

        estables = Vector{Int}(undef, medidas)  # Prealoca para evitar push!

        for j in 1:medidas
            estables[j] = integrador(params_local)
        end

        sol[idx] = count(x -> x == 1, estables)
    end

    return sol
end


"""
Esta función es la misma que se ubica en Sistema.jl solamente que se encuentra implementada con Threads. 
Vamos a probar cómputo en paralelo para revisar si se pueden optimizar tiempos de ejecución.
"""

function transicionσParalel(params::Parametros,σ)
    sol = Vector{Int}(undef,length(σ))
    medidas::Int = 1000
    
    Threads.@threads for idx in 1:length(σ)
        params_local = deepcopy(params)
        params_local.σ = σ[idx]
        
        estables = Vector{Int}(undef,medidas)

        for j in 1:medidas 
            estables[j] = integrador(params_local)  
        end
        sol[idx] = count(x -> x == 1, estables)
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
        jacobianos = []
        params.p = i 

        # Crear el lock fuera del ciclo
        lck = ReentrantLock()

        # Paralelización en la generación de diagonales
        Threads.@threads for i in 1:medidas
            M = Jacobianos(params::Parametros)

            # Bloque crítico con el lock
            lock(lck)  # Bloqueamos antes de modificar el arreglo compartido
            try
                push!(jacobianos, M)
            finally
                unlock(lck)  # Desbloqueamos después de la operación
            end
        end

        # Escritura secuencial de los resultados a un archivo CSV
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