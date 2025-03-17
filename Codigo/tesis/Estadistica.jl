"""
Función para determinar las modas de una distribucion con base en un histograma, esta función
se ha utilizado como alternativa para determinar las modas en lugar de la función mode()
"""

function Modas(distribucion::Array,bins::Int)
    hist = fit(Histogram, distribucion, nbins=bins)  
    idx = argmax(hist.weights)
    mode_bin = hist.edges[1][idx]  # Limite izquierdo del bin más frecuente
    return mode_bin
end