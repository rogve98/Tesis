mutable struct DatosCSV
    dominio::Union{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},Vector{Float64}}
    resultado::Vector{Any}
    ruta::String
    nombeArchivo::String
end


"""
Función para guardar las corridas de la función de transicion en el archivo 
correspondiente CSV.
"""

function transicionCSV(datos::DatosCSV)
    ruta = datos.ruta
    nombre = datos.nombeArchivo
    p = datos.dominio
    resultado = datos.resultado
    sN = DataFrame(CSV.File(ruta,header=false))   
    sN = sN[!,2]
    sNuevo = []
    for i in 1:length(p)
        push!(sNuevo,sN[i]+resultado[i])
    end
    writedlm(nombre,[p sNuevo])
end


"""
Función para explorar las diagonales de los Jacobianos según sus parámetros
"""

function vectorDiagonales(ruta::String,DF::DataFrame)
    if ruta != ""
        df = DataFrame(CSV.File(ruta,header=false))
        df = permutedims(df)
        return df,vcat(collect(eachcol(df))...)    
    else
        DF = permutedims(DF)
        return DF,vcat(collect(eachcol(DF))...)
    end
end

"""
Analizar parámetro críctico, en términos de las diagonales de las matrices jacobianas
"""

function propiedadesJacobianas(σ::Float64,propiedad::String)
    Jbs = []
    p = range(0,stop=1,length=100)
    for i in p
        ruta = "Jacobianos_s$(σ)_p$i.csv"
        if isfile(ruta)
            df = CSV.read(ruta,DataFrame,header=false)
            push!(Jbs,df)
        end
    end
    Jbs = Array.(Jbs)
    matrices = []
    for i in eachindex(Jbs)                                                                                                                                      conjuntos = []
        for j in eachindex(Jbs[i][:,1])
            ms = reshape(Jbs[i][j,:],100,100)
            push!(conjuntos,ms)
        end          
        push!(matrices,conjuntos)
    end
    diags = []
    for i in eachindex(matrices)
        push!(diags,diag.(matrices[i]))
    end
    if propiedad == "trazas"
        trazas = []
        for i in eachindex(diags)
            push!(trazas,sum(sum.(diags[i])))
        end
        return trazas
    elseif propiedad == "trazas-promedio"
        trazas = []
        for i in eachindex(diags)
            push!(trazas,mean(sum.(diags[i])))
        end
        return trazas 
    elseif propiedad == "promedio" 
        promedio = []
        for i in eachindex(diags)
            push!(promedio,mean(mean.(diags[i])))
        end
        return promedio 
    elseif propiedad == "energia"
        energia = []
        for i in eachindex(matrices)
            push!(energia,-mean(sum.(matrices[i])))
        end
        return energia
    elseif propiedad == "Z-particion"
        diagonales = []
        zParticion = []
        for i in eachindex(diags)
            push!(diagonales,mean(sum.(diags[i])))
            push!(zParticion,mean(sum.(diags[i])))
        end
        return diagonales,zParticion
    elseif  propiedad == "f-energiaLibre"
        diagonales = []
        fELibre = []
        for i in eachindex(diags)
            push!(diagonales,mean(sum.(diags[i])))
            push!(fELibre,log(abs(mean(sum.(diags[i])))))
        end
        return diagonales,fELibre
    elseif propiedad == "derivada-f"
        diagonales = []
        derivada = []
        for i in eachindex(diags)
            push!(diagonales,mean(sum.(diags[i])))
            push!(derivada,mean(1/sum.(diags[i])))
        end
        return diagonales,derivada
    else
        return "Opción no válida"
    end
end