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

        sol[idx] = sum(estables)
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
        sol[idx] = sum(estables)
    end
    return sol    
end

"""
Código alternativo para utilizar Threads con el fin de minimizar tiempos de ejecución usando varios
cores en la ejecución. La ejecución para σ = 0.5 tardó 9 días! Entonces para sigmas superiores 
esto se va a volver un relajo. Se busca minimzar tiempos de ejecución
"""

function distDiagonalParalel(params::Parametros, p)
    medidas = 100
    σ = params.σ
    # Preasignar el arreglo de resultados (uno para cada valor de p)
    jacobianoss_total = Vector{Vector{Matrix{Float64}}}(undef, length(p))

    Threads.@threads for idx in 1:length(p)
        i = p[idx]
        jacobianoss = []  # Este vector será local a cada hilo
        params.p = i
        
        for _ in 1:medidas
            M = Jacobianos(params)
            push!(jacobianoss, M)
        end
        
        # Asignar los resultados del hilo a su respectivo índice
        jacobianoss_total[idx] = jacobianoss
    end

    # Escribir los resultados a archivos CSV después de la paralelización
    for (i, jacobianos) in enumerate(jacobianoss_total)
        writedlm("Jacobianos_s$σ.p$(p[i]).csv", jacobianos)
    end
end


"""
TRansición paralela con jacobianos y matrices de incidencais incluidas
"""

function transicionParalelCompleto(params::Parametros, p)
    sol = Vector{Int}(undef, length(p))  # Vector preasignado
    medidas = 1000

    Threads.@threads for idx in eachindex(p)
        params_local = deepcopy(params)  # Evita conflictos entre hilos
        params_local.p = p[idx]

        estables = Vector{Int}(undef, medidas)  # Prealoca para evitar push!
        Jbs = []
        IncEstables = []
        IncInestables = []

        for j in 1:medidas
            #solJ, solIn, Incidencias, Jacobiano
            e, J, Λ = integradorJacobs(params_local)
            estables[j] = e
            
            if e == 1
                push!(Jbs, J)
                push!(IncEstables,Λ)
            else
                push!(IncInestables,Λ)
            end
        end

        sol[idx] = sum(estables)

        # Guardar los resultados solo si hay datos
        if !isempty(Jbs)
            writedlm("Jacobianos_s$(params.σ)_p$(p[idx]).csv", Jbs)
        end
        if !isempty(IncEstables)
            writedlm("IncidenciasEstables_s$(params.σ)_p$(p[idx]).csv", IncEstables)
        end
        if !isempty(IncInestables)
            writedlm("IncidenciasInestables_s$(params.σ)_p$(p[idx]).csv", IncInestables)
        end
    end
    return sol
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