# Resumen
Modelar, simular y representar redes metabólicas aún hoy en día representa un reto. Se han propuesto diversas aproximaciones para abordar este problema, los cuales se pueden resumir en tres principales modelos: los modelos continuos, discretos e híbridos. En este trabajo se propuso modelar las redes metabólicas con redes metabólicas con redes de Petri. Con este objetivo se implementó Metapenta, una herramienta que es capaz de realizar diferentes consultas sobre las redes metabólicas como encontrar el camino mínimo entre dos metabolitos, identificar metabolitos claves como los sumideros y las fuentes,  visualizar las reacciones de un metabolito dato por parámetro, entre otras funcionalidades.​

## Contexto biologico
Una red metabólica son el conjunto de todas las reacciones que ocurran en una célula. Está compuesta de tres componentes principales: metabolitos, reacciones y enzimas. Entender, visualizar y analizar estas redes metabólicas tiene una gran aplicación en la industria ya que ayuda a entender mecanismos de enfermedades como la diabetes y a mejorar producción de ciertos alimentos tales como pan y cerveza. En este proyecto se planteó una solución para analizar redes metabólicos basados en el formalismo de redes de Petri.

## Metodologia
Se implementó el acercamiento discreto al modelo de redes metabólicas. Para esto se montó la red biológica sobre una red de Petri, y sobre esta se hicieron las consultas. Todo esto constituyó la lógica del modelo. Toda la visualización de la red, se hizo usando Processing 2.0. Sin embargo, el modo como este modulo debía recibir una forma de la red que fuera apta pare ser pintada, por lo que hay una clase especializada en hacer esa traducción. Por otro lado, toda la interfaz gráfica se implementó utilizando JavaFX y este mismo fue embebido en Processing para de manera que se combinarán ambas tecnologías. 

## Diagrama de clases
![MetaPenta(2)](https://user-images.githubusercontent.com/32238112/120860438-4c677600-c54b-11eb-97ed-ac45ea69ec99.png)

## Dependencias
1. Java 11
2. [Processing](https://processing.org/): Una librería para realizar gráficos de manera gráfica, se utilizó para realizar la representación de lad metabólica de la red de Petri

## Interfaz de usuario

### Panel 1
Todo este panel está hecho en Processing 2.0 y tiene una lógica completamente independiente del resto de la aplicación, aquí sólo se maneja visualización. El usuario puede mover y arrastrar los metabolitos y las reacciones que están en este panel para lograr la conformación de la red que más le facilite el análisis. Adicionalmente, según los requerimientos, los nodos pueden mostrarse de distintos colores. Por ejemplo, en el requerimiento de ruta mínima, los metabolitos iniciales se muestran en naranja, el metabolito objetivo se muestra en verde y los metabolitos intermedios se muestran en morado.

### Panel 2
En este panel es donde ocurre la mayor parte de la interacción con el usuario. En la parte superior del panel, el usuario puede introducir los IDs de los metabolitos para hacer la búsqueda de la ruta mínima. Si bajamos un poco, llegamos al campo de texto que nos permite filtrar y mostrar la subnet asociada a un metabolito. Para realizar este filtro, el usuario debe colocar el ID del metabolito y seleccionar si quiere que se muestren las reacciones donde el metabolito es substrato, producto o ambas en el checkbok que se encuentra justo debajo. En la última parte de este panel podemos seleccionar, también por medio de un checkbok, los documentos que queremos descargar y el programa generará los archivos en la ruta con el prefijo especificado cuando el usuario oprima el botón download.

### Panel 3
Este panel es el que tiene menos interacción directa con el usuario debido a que, sólo se encarga de mostrar la información del nodo seleccionado. Si es un metabolito, se muestra su id, el compartimento, su fórmula química y las reacciones donde participa. Por otro lado, si es una reacción se muestra su nombre y los IDs de los metabolitos que en ella participan.

### Panel 4
Este panel tiene dos botones. El primero, el de \textit{"load"}, permite cargar el archivo en formato XML de la red metabólica, mientras que el botón de \textit{"clean"} elimina la red que se encuentre cargada para una nueva red pueda ser cargada.

![ui](https://user-images.githubusercontent.com/32238112/120860331-217d2200-c54b-11eb-922e-02896dc9a20e.png)

## Running
Si estás utilizando un IDE (por ejemplo eclipse), puedes correrla la aplicación corriendo la clase principal de la aplicación:

https://github.com/FLAGlab/MetaPeNTA/blob/master/src/metapenta/gui/Main.java

Si lo quieres correr por línea de comandos necesitas los archivos `jar` de la libería de `Processing` y los binarios de `Java FX` y correr el siguiente comando

`C:\Program Files\Java\jdk-11.0.7\bin\javaw.exe -Dfile.encoding=Cp1252 -classpath “<List of jar files>” metapenta.gui.Main <metabolic_file_path> <out_file>`

Donde <List of jar files> son los archivos antes mencionados. Los `jar` de processing, pueden encontrarse en la carpeta `lib` 

## License

MIT License

Copyright (c) 2020 vparrac

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
