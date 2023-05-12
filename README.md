# Tesis
Este repositorio contiene todos los archivos trabajados durante mi proceso de tesis. El tema es transiciones de fase en sistemas de N especies en competencia.

La estructura del reposotorio es la siguiente

* En la carpeta [Codigo](/master/Codigo) se encuentra los archivos .jl del proyecto en donde existe uno de pruebas y otro definitivo llamado main.jl
* En la carpeta de [Notebooks](/master/Notebooks) se encuentran cuatro notebooks, dos de apoyo para pintar espacios fase con Makie.jl y con matplotlib.py que respectivamente son Elementos básicos y Espacios fase. El notebook de pruebas son donde se hallan las pruebas y errores del código y en avances se encuentran los resultados importantes hasta el momento.
* En la carpeta [Imagenes](/master/Imagenes) se encuentran todas imágenes generadas y guardadas de los experimentos y avances que se hagan.

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