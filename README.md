# Tesis
Este repositorio contiene todos los archivos trabajados durante mi proceso de tesis. El tema es transiciones de fase en sistemas de N especies en competencia.

La estructura del reposotorio es la siguiente

* En la carpeta [Codigo](/Codigo) se encuentra los archivos .jl del proyecto en donde existe uno de pruebas y otro definitivo llamado main.jl. Recientemente se agregó uno de avances.jl que contiene el mismo código que main.jl solo que usa PyPlot como graficador. Esto para poder pintar espacios fase y otras cosas.
* En la carpeta de [Notebooks](/Notebooks) se encuentran cuatro notebooks, dos de apoyo para pintar espacios fase con Makie.jl y con matplotlib.py que respectivamente son Elementos básicos y Espacios fase. El notebook de pruebas son donde se hallan las pruebas y errores del código y en avances se encuentran los resultados importantes hasta el momento.
* En la carpeta [Imagenes](/Imagenes) se encuentran todas imágenes generadas y guardadas de los experimentos y avances que se hagan.

## Avances

### Elementos introductorios

Se dice que en nuestro universo no existe nada que premanezca en un reposo constante o absoluto, las entidades que componen el universo están en constante movimiento; y el movimiento ha sido ampliamente estudiado desde los tiempos de los griegos. Con el paso del tiempo, la evolución del conocimiento nos llevó a cuantificar el movimiento con base en las leyes de Newton, más precisamente con la segunda que nos dice

$$F=ma$$

De aquí surgirán modelos para poder estudiar el movimiento, tales como péndulos, planos inclinados etc. Pero para movimientos más complejos esta perspectiva puede llegar a ser muy laboriosa; como respuesta a ello surgieron dos interpretaciones de la mecánica las cuales son la formulación Lagrangiana y Hamiltoniana, que lo que hacen es ordenar el problema en un conjunto de ecuaciones diferenciales (la cantidad y el orden dependen del caso) y con ello pasar a resolver de forma más simple.

Con estas formulaciones los sistemas que se van estudiando son cada vez más difíciles pero capaces de resolverse, la única condición para eso es que las ecuaciones que compongan la dinámica sean lineales, es decir, que las ecuaciones tengan puros términos lineales. Que un sistema sea lineal implica que cumplen el principio de superposición, el cual nos dice que para cualquier conjunto de soluciones, la suma de las mismas también es una solución.

Sin embargo en la naturaleza es difícil observar sistemas con comportamiento lineal, situaciones en donde la suma de las partes sea igual al total. Podemos proponer sistemas bien comportados como los péndulos o cualquier tipo de oscilación, siempre cumplirán el principio de superposición. Por otra parte, la observación nos ha permitido constatar que el total de un sistema no es igual a la suma de sus partes, y el primer ejemplo evidente son los seres vivos, que compuestos de millares de células somos más complejos en conjunto que por separado.

De aquí surge un problema, porque las ecuaciones de las que hablábamos anteriormente no son capaces de resolver este tipo de dinámicas, y de aquí surgen las ciencias de la complejidad sostenidas por la dinámica no lineal. Cuando las ecuaciones que gobiernan algún fenómeno llegan a ser no lineales ocurren una serie de cosas: primero no pueden ser resueltas de manera analítica, existen otras técnicas (computacionales evidentemente); producen caos, es decir, comportamientos no periódicos en el tiempo y difíciles de predecir, sin embargo las ecuaciones siguen siendo deterministas. 

La no linealidad viene a partir de las interacciones entre las variables en cuestión, ya sean las células, individuos, activos etc. Esta comunicación afecta en el comportamiento de cada una de las variables en cuestión de forma que sea complicado cuantificar que esta sucediendo en tiempo real. Los modelos computacionales son una forma de estimar los comportamientos de este tipo de sistemas, sin embargo deben ser considerados como una aproximación, la vida real es más compleja que ello.

Para poder resolver este tipo de problemas es necesario integrar numéricamente las ecuaciones no-lineales que componen el sistema, para ello usamos integradores como Euler o Runge-Kutta de orden cuatro. Por tanto hay que pelearnos con el lenguaje para hacer que la computadora nos entienda de manera adecuada y nos resuelva lo que queremos en cuestión.


### Avances en el código

**050523**

Ya que tenemos nuestro integrador de RK4, definimos un sistema $N-$ dimensional con sus interacciones pertinentes. Para ello podemos agarrarnos el modelo de red aleatoria de Erdös-Rènyi (modelo GNP) que es un modelo que atribuye un enlace (interacción) en función de una probabilidad. Extraemos su matriz de adyacencia y la manipulamos para poder obtener una matriz aleatoria con entradas que cumplan la distribución normal. 

El problema no lineal que queremos resolver es el de un sistema de $N$ especies en competencia, por tanto nos agarramos las ecuaciones de Lotka-Volterra en dicha versión. Armando el sistema notamos que funciona para dos especies muy bien pero para para cinco comienza a hacer cosas extrañas como sobrepasar la capacidad de carga, descubrimos de forma empírica que se debía a la escala en la que se encuentran los números de la distribución y la diagonal que representan los términos logísticos de las ecuaciones.

Elegimos un factor de escala de 10 y dicho problema se ha resuelto pero hace falta demostrar porque funciona.

**110523**

He determinado el jacobiano del sistema en su forma general, deben cumplir dos reglas determinadas para la identidad y para la parte triangular superior e inferior. He generado código al respecto y hace falta verificar que funcione. Ahora hace falta determinar por medio de un algorítmo los puntos críticos del sistema y posteriormente su estabilidad asociada; luego vendrá lo de la ley circular y el resto de la chamba.

**170523**

Una vez generado el jacobiano del sistema hay que determinar los puntos críticos del mismo. Como tenemos un sistema de ecuaciones no lineales, y queremos encontrar raíces: lo ideal sería generar un algorítmo de Newton Rhapson n-dimensional para generar las raíces. Para ello se debe generar la siguiente regla.

$$
v_{i+1}=v_i-(J(v_i))^{-1}f(v_i)
$$

Donde las $v's$ es el vector solución en donde tiene la raíz de cada entrada, y $J(v_i)$ es el Jacobiano del sistema el cual hemos determinado con anterioridad. El avance de esta semana es que se generado dicho algorítmo y se ha probado para un sistema de $2\times 2$ que ya se había estudiado con cuidado anteriormente y las soluciones de sus puntos críticos coinciden con los cálculos analíticos. Quizás valga la pena analizar un sistema de $3\times 3$ para observar como se comportan los puntos críticos en dimensiones superiores.

Lo anterior con el fin de ver como vienen configurados los puntos críticos, si se quedan en un plano o si andan combinados por el hiper espacio etc etc.

**180523**

Se ha comprobado que el código funciona para sistemas de $2\times 2$. Se quizo extender hacia sistemas de $3\times 3$, pero no supe como hacer los espacios fase en 3D. Agregué todos los avances de ahora al notebook de avances donde viene todo un resumen. Lo siguiente será ir experimentando con las probabilidades (conectancia) y la magnitud de las interacciones. De ahí en adelante el camino es desconocido.

**220523**

Ahora viene la siguiente fase, necesitamos comprender como se dan las transiciones de fase en este tipo de sistemas. Recordemos que uno de los artículos propone la emergencia de varios puntos críticos atractores a partir de ciertos parámetros, mientras que el otro artículo reporta velocidades de mutación para cierto punto crítico. Hay que checar con cuidado ambos artículos para ver que podemos recoger de ellos. 

**290523**

El avance de ahora fue generar los primeros inputs sobre la ley circular y graficar los eigenvalores.

**030623**

Avances en la ley circular, ya tenemos forma de alterar la distribución normal con diferentes desviaciones estándar. Se hicieron pruebas con desviaciones superiores a 1 pero no se encontró lo esperado sino que solo se halló un desplazamiento del círculo más no una muestra de inestabilidad.

# Retomando...

04/01/2023

Tuvimos un largo periodo de inactividad, 6 meses al menos, la última entrada fue el 3 de julio cuando andábamos saliendo de "vacaciones" y Sergio andaría fuera de descanso. Ese descanso se prolongó, comencé en la música, en el trabajo etc etc.

Vamos a retomar pero con un nuevo enfoque, el enfoque que teníamos estaba bien, ya están construidas varias cosas pero la dirección que se estaba teniendo no me estaba convenciendo del todo; más que nada porque estábamos dejando fuera el tema de las transiciones de fase que es mi punto central de la investigación; nos estábamos enfocando en la estabilidad del sistema que si bien es interesante, veía un poco conflictiva la forma en la que se vería ante el sínodo si puede que "no exista" mucha física al respecto. 

Por ello desde noviembre del año pasado estuve pensando en un nuevo enfoque que me regresara las transiciones de fase y justo llegué a un clímax. Vamos a mezclar la teoría de redes enfocada en sistemas complejos con el sistema de dinámica no lineal y con eso trabajaremos para investigar una transición de fase pero debe ser teórica.

Los elementos principales ahora son:

1. Construir una red de Albert-Barabasi, es decir, una red con independencia de escala, modelar su ley de potencias y de ahí ligar los sistemas (matrices de adyacencia) con los sistemas de dinámica no lineal.

Esta es la principal misión, por eso quiero generar una meta para Enero: Construir el modelo para finales de Enero y avanzar a la siguiente fase que sería el ligamento con la parte del sistema no lineal.
