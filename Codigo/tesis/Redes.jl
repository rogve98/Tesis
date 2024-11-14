""" Modelo de red aleatoria, la primera función es un generador de enlaces enlacesAleatorios con
base en un número de enlaces N y una probabilidad p """

function enlacesAleatorios(N,p)
    Channel() do channel
        for i in 1:N
            for j in 1:i
                if rand() < p
                    if i == j
                        continue
                    else
                        put!(channel,(i,j))
                    end
                end
            end
        end
    end
end

""" Generamos la red aleatoria con el generador de enlaces aleatorios. Ya sea dirigida o 
no dirigida."""

function redAleatoria(N,p,red::String)
    if red == "dirigida"
        g = DiGraph(N)
        enlaces_i = collect(enlacesAleatorios(N,p))
        enlaces_j = collect(enlacesAleatorios(N,p))
        while length(enlaces_j)!=length(enlaces_i)
            enlaces_j = collect(enlacesAleatorios(N,p))
        end   
        for i in 1:length(enlaces_i)
            add_edge!(g,enlaces_i[i][1],enlaces_i[i][2])
            add_edge!(g,enlaces_j[i][2],enlaces_j[i][1])
        end
        return g
    elseif red == "no dirigida"
        g = Graph(N)
        enlaces = collect(enlacesAleatorios(N,p))
        for i in 1:length(enlaces)
            add_edge!(g,enlaces[i][1],enlaces[i][2])
        end
        return g
    end
end

"""Esta función nos genera matrices aleatorias con base en el modelo de red aleatoria.
Lo hice de esta manera para que pudieramos tener una representación gráfica de las 
interacciones y al mimsmo tiempo obtener matrices aleatorias con valores aleatorios además."""

function randomMatrix(N,p,σ,red::String)
    d = Normal(0,σ)
    g = redAleatoria(N,p,red)
    M = adjacency_matrix(g)
    Id = 1* Matrix(I, N, N)
    M = M.*rand(d,N,N)+I
    # for i in 1:N
    #     M[i,i] = 1
    # end
    return (Matrix(M), g)
end

"""
Combina las tres funciones anteriores para poder optimizar los tiempos de compilación. Utiliza
el struct Parametros como argumento.
"""

function incidencias(parametros::Parametros)
    N = parametros.N
    p = parametros.p
    σ = parametros.σ
    red = parametros.Red
    d = Normal(0,σ)

    function enlacesAleatorios(N,p)
        Channel() do channel
            for i in 1:N
                for j in 1:i
                    if rand() < p
                        if i == j
                            continue
                        else
                            put!(channel,(i,j))
                        end
                    end
                end
            end
        end
    end
    # Se puede optimizar más dejando solo el caso "no dirigida".
    if red == "dirigida"
        g = DiGraph(N)
        enlaces_i = collect(enlacesAleatorios(N,p))
        enlaces_j = collect(enlacesAleatorios(N,p))
        while length(enlaces_j)!=length(enlaces_i)
            enlaces_j = collect(enlacesAleatorios(N,p))
        end   
        for i in 1:length(enlaces_i)
            add_edge!(g,enlaces_i[i][1],enlaces_i[i][2])
            add_edge!(g,enlaces_j[i][2],enlaces_j[i][1])
        end
    elseif red == "no dirigida"
        g = Graph(N)
        enlaces = collect(enlacesAleatorios(N,p))
        for i in 1:length(enlaces)
            add_edge!(g,enlaces[i][1],enlaces[i][2])
        end
    end

    M = adjacency_matrix(g)
    M = M.*rand(d,N,N)
    for i in 1:N
        M[i,i] = 1
    end
    return Matrix(M)

end


"""Red Circular"""

function Circular(N)
    g = Graph(N)
    for i in 1:N
        add_edge!(g,i,i+1)
    end
    add_edge!(g,N,1)
    return g#gplot(g,nodelabel=1:N)
end

"""
n:=número de nodos totales
m:=número de enlaces nuevos por nodo; con este parámetro definimos la red inicial a partir de una red GNL

La idea del algorítmo es construir una red inicial aleatoria tipo GNL para asegurarnos que sea mayormente conexa
(la GNP no lo asegura al 100% a menos de que ocupemos la teoría de punto crítico pero eso lo veremos después). 
Se verifica si es conexa mediante un condicional y si si generamos el algorítmo de albert barabasi que consiste en
definit probabilidades de enlace en función del grado de cada nodo j dividido entre el grado total (suma de grados).
De aquí definimos un vector de m entradas que serán los nuevos nodos a los que irán conectados un último y nuevo 
nodo.

targets funciona de tal manera que hace una selección desde el nodo 1 hasta el nodo nv(G) y por medio de las
probabilidades se van eligiendo esos nodos. Por ejemplo para la primera entrada del vector será más probable que
se escoga el número entre 1 y nv(G) que tenga la probabilidad más alta en probabilities. Por ejemplo si el vector
de probabilities en su tercera entrada tiene un grado alto (probabilidad alta), es muy probable que el nodo que 
vaya a salir seleccionado entre 1 y nv(G) es el 3 (porque le corresponde la tercera entrada). De aquí actúa el 
replace y el 3 ya no puede ser escogido en futuras ocasiones pero el procedimiento se repite para la siguiente 
entrada del vector targets, y así hasta acabar con todos.

Posteriormente se agrega un nodo extra para poder conectar los targets a este nodo extra y la manera de conectarlos
es por medio del índice i del ciclo for que justamente coincide con el nodo agregado, de aquí nada más falta
conectar ese nodo extra con todos los nodos de target. Y bueno este procedimiento se realiza iteradamente hasta
terminar con los valores de n.
"""

function barabasi_albert(n, m)
    # Inicializar un grafo completo con m nodos
    G = Circular(m) #GNL(m,m)
    if is_connected(G)
        # Implementar el modelo de Barabási-Albert
        for i in m+1:n
            # Calcular las probabilidades de conexión para nodos existentes
            probabilities = [degree(G, j) / sum(degree(G, k) for k in vertices(G)) for j in vertices(G)]

            # Elegir m nodos existentes basados en las probabilidades
            targets = sample(1: nv(G),Weights(probabilities),m,replace=false)

            # Agregar un nuevo nodo
            add_vertex!(G)

            # Conectar el nuevo nodo a los nodos seleccionados
            for target in targets
                add_edge!(G,i,target)
            end
        end
    end
    M = Matrix(adjacency_matrix(G))
    M = 10*M.*randn(N,N)
    for i in 1:N
        M[i,i] = 1
    end

    return G,M
end

"""

"""

function BMatrix(N,m,σ)
    d = Normal(0,σ)
    g,_ = barabasi_albert(N,m)
    M = adjacency_matrix(g)
    M = M.*rand(d,N,N)
    for i in 1:N
        M[i,i] = 1
    end
    return (Matrix(M), g)
end