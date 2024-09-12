mutable struct DatosCSV
    dominio::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
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

function vectorDiagonales(ruta::String)
    df = DataFrame(CSV.File(ruta,header=false))
    df = permutedims(df)
    return df,vcat(collect(eachcol(df))...)    
end