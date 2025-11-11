# NombreApellido_Trabajo2.R
# Pablo Badolato López-Tormos
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento
#install.packages("vioplot")
#install.packages("ggplot2")

# 1. Cargar librerías (si necesarias) y datos del archivo "datos_biomed.csv". (0.5 pts)
# Cargo la librería para poder leer el archivo
library(readr)
#Cargo la librería para poder hacer el violin plot
library(vioplot)
#Cargo la librería del ggplot2, porque he visto que queda mucho mejor que usando funciones base de R
library(ggplot2)
# Abro el archivo
archivo <- read_csv("datos_biomed.csv")

# 2. Exploración inicial con las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? (0.5 pts)
#Hay 5 variables y 3 tratamientos
head(archivo) #Me da las 6 primeras filas 
summary (archivo) #Me da un resumen del archivo con información estadística y acerca del tamaño 
dim (archivo)#Me informa del número de filas y columnas (en este orden)
str(archivo) #Me da información acerca del tipo de variable con ejemplos. Es decir, me muestra la estructura.

# 3. Una gráfica que incluya todos los boxplots por tratamiento. (1 pt)
# Selecciono las variables que quiero comparar (Glucosa, presión y colesterol)
variables <- archivo[, c("Glucosa", "Presion", "Colesterol")]
#Para convertir "Tratamiento" en un factor y así asegurarme de que R no confunde el grupo de tratamientos como un número.Además, permite que los gráficos separen los datos por categoría
archivo$Tratamiento <- factor(archivo$Tratamiento)
# Para que los 3 boxplot queden juntos en una misma imagen pero diferenciados uno de otro (como si estuvieran en la misma fila pero en distintas columnas)
par(mfrow = c(1, 3))   # 1 fila, 3 columnas
# Creo los boxplot uno a uno usando un bucle for
for (variable in names(variables)){
	boxplot(variables[[variable]]~archivo$Tratamiento, #"boxplot": indica que la función que voy a realizar es generar un gráfico de esta forma. "~archivo$Tratamiento": Esto dice que coja las variables (que he establecido antes), en función de tratamiento 
	main = variable, #Establece el título de cada boxplot (el título será el nombre de la variable dependiendo de en cual esté el bucle)
	xlab = "Tratamiento", #El título del eje x
	ylab = variable, #El título del eje y
	col = c("red", "blue", "green")) #Le pongo un color diferente a cada boxplot
}
# 4. Realiza un violin plot (investiga qué es). (1 pt)
#Igual que antes, para ver los 3 gráficos juntos
par(mfrow = c(1, 3))   # 1 fila, 3 columnas
#Igual que antes, para asegurarme
archivo$Tratamiento <- factor(archivo$Tratamiento)
#Voy a hacer un violin plot que tenga fusionados los 3 como he hecho con el boxplot, mediante un bucle for
for (variable in names(variables)){
	vioplot(variables [[variable]]~archivo$Tratamiento, # "vioplot": es la función que indica que voy a hacer un violin plot. El resto es como en el gráfico anterior
	col = c("red", "lightblue", "green"), #establezco los colores
	main = paste("Distribución de ", variable), # Título, la función paste es parecida al print(f"") en python, que me permite incluir una variable en un str
	xlab = "Tratamiento", #Nombro al eje x
	ylab = variable) #Nombro al eje y
}

# 5. Realiza un gráfico de dispersión "Glucosa vs Presión". Emplea legend() para incluir una leyenda en la parte inferior derecha. (1 pt)
#Igual que antes, para ver los 3 gráficos juntos
archivo$Tratamiento <- factor(archivo$Tratamiento)
# Definir una paleta de colores, uno por nivel de Tratamiento
colores <- c("red", "blue", "green")
# Asignar un color a cada punto según su grupo
colores_trat <- colores[archivo$Tratamiento]
# Creo el gráfico de dispersión
plot(archivo$Glucosa, archivo$Presion,
	col = colores_trat, # Le asigno un color a cada grupo
	pch = 19, # Para ver puntos rellenos que se vean mejor
	main = "Relación entre Glucosa y Presión", #El título de la gráfica
	xlab = "Glucosa", #El título del eje x
	ylab = "Presión" #El título del eje y
)
# Añado la leyenda en la parte inferior derecha
legend("bottomright", #Pongo la leyenda abajo a la derecha
	legend = unique(archivo$Tratamiento), #Digo que cada punto corresponde a un paciente/sujeto que ha sido tratado con un tipo de tratamiento
	col = colores,
	pch = 19)


# 6. Realiza un facet Grid (investiga qué es): Colesterol vs Presión por tratamiento. (1 pt)
#Igual que antes, para ver los 3 gráficos juntos
archivo$Tratamiento <- factor(archivo$Tratamiento)
ggplot(archivo, aes(x = Colesterol, y = Presion)) + #Crea la base del gráfico al definir los ejes y de donde se sacan los datos. Además, aes (de "aesthetic mappings") lo que hace es asignar mis variables a elementos visuales (dice qué se dibuja)
	geom_point(alpha = 0.8) + #geom_point: agrega los puntos de dispersión, mientras que "alpha = 0.8" hace que sean un poco transparentes, en caso de que se superpongan
	geom_smooth(method = "lm", se = FALSE) + #Con geo_smooth añado una línea de tendencia de modelo lineal (lm)en cada panel. Con se = FALSE lo que consigo es que no se sombree el intervalo de confianza en la línea
	facet_grid(. ~ Tratamiento) + #facet_grid() divide el gráfico en una rejilla, y ". ~ Tratamiento" es para que se vea en horizontal en vez de en vertical, que sería "Tratamiento ~ ." 
	labs(title = "Colesterol vs Presión por tratamiento", # labs() define qué se etiqueta en cada eje (nombra los ejes, no es como aes que define lo que se dibuja)
	x = "Colesterol", y = "Presión") +
	theme_minimal()#Aplica un estilo visual limpio, sin bordes ni fondo gris, ampliamente utilizado para presentaciones o artículos científicos

# 7. Realiza un histogramas para cada variable. (0.5 pts)
# Configuro la ventana gráfica: 1 fila, 3 columnas (igual que en los ejercicios anteriores)
par(mfrow = c(1, 3))
#Creo el histograma a través de un bucle for para generar los 3 histogramas (uno por cada variable de variables, definidas en el ejercicio 3)
for(variable in names(variables)){
	hist(variables[[variable]], #Esto define la función, que ahora será un histograma (hist)
	main = paste("Histograma de ", variable), #Define el título de cada histograma, usando la función paste() explicada anteriormente
	xlab = variable, #Etiqueto los ejes 
	ylab = "Frecuencia")
}

# 8. Crea un factor a partir del tratamiento. Investifa factor(). (1 pt)
#Con la función "factor()" lo que consigo es que la columna "Tratamiento", compuesta por "FarmacoA", "FarmacoB" y "Placebo", en vez de ser leída como texto normal (character) que sea leída como grupos. De esta forma trata las 3 categorías de la columna de Tratamiento como niveles. Con este código lo que hago es sustituir la versión "texto" de la columna por la versión "factor"
archivo$Tratamiento <- factor(archivo$Tratamiento)
#Para ver la estructura del archivo, donde puedo observar "Tratamiento: Factor w/ 3 levels "FarmacoA","FarmacoB"...", lo que me indica que se han creado bien los niveles
str(archivo)
#De esta forma puedo ver los niveles que he creado
levels(archivo$Tratamiento)

# 9. Obtén la media y desviación estándar de los niveles de glucosa por tratamiento. Emplea aggregate() o apply(). (0.5 pts)
#Voy a usar la función "aggregate()"  en vez de "apply()", ya que esta última trabaja con matrices o data frames numéricos, no con factores directamente. Por lo tanto, si quisiera usar esta función tendrías que dividir los datos antes. Además, tener Tratamiento como factor es lo ideal para este tipo de cálculos
#La función aggregate() tiene la siguiente estructura básica: aggregate(variable ~ grupo, data = archivo, FUN = función). Esta función permite agrupar los datos de "variable" por "grupo" (Glucosa)la función que le indique (media y desviación estándar
#Uso "agreggate()" para la media de glucosa por tratamiento
aggregate(Glucosa ~ Tratamiento, data = archivo, FUN = mean)
#Ahora la desviación estándar
aggregate(Glucosa ~ Tratamiento, data = archivo, FUN = sd)

# 10. Extrae los datos para cada tratamiento y almacenalos en una variable. Ejemplo todos los datos de Placebo en una variable llamada placebo. (1 pt)
#Voy a usar la función "split()" para divir los datos. La forma general de esta función es: split(datos, grupo), donde datos corresponde al data frame y grupo corresponde a cómo quiero dividir el data frame. Con esta función se crea una lista por cada subconjunto de datos que consigas al dividir los datos
lista_tratamientos <- split(archivo, archivo$Tratamiento)
#Para ver la lista, que ya no es un data frame, sino una lista con varios “mini data frames” dentro
str(lista_tratamientos)
#Ahora para ver los 6 primeros de cada grupo, uso un bucle for para ir yendo grupo a grupo. Además, en el bucle incluyo la función "cat()" para que se ponga el nombre del grupo como encabezado y así mejorar la visualización de los datos
for (grupo in names(lista_tratamientos)) {
  cat("\n-----", grupo, "-----\n")   # Muestra el nombre del grupo con un salto (\n) y entre - para verlo más fácil y bonito
  print(head(lista_tratamientos[[grupo]]))  # Muestra las primeras 6 filas
}

# 11. Evalúa si los datos siguen una distribución normal y realiza una comparativa de medias acorde. (1 pt)
#Lo primero que hay que hacer es evaluar la normalidad de los datos. Yo voy a usar el test de Shapiro-Wilk, que es el test más usado (o eso tengo entendido) y se escribe tal que así: shapiro.test(valores)
#Para ir variable a variable (Glucosa, Presión y Colesterol) y medir la normalidad de todos los datos, creo un bule for para ir variable a variable en la lista vars. Y uso la función "cat()" para ir visualizando mejor qué estoy calculando
vars <- c("Glucosa", "Presion", "Colesterol")
for (variable in vars) {
  cat("\n### Variable:", variable, "###\n")
	# Bucle interno para recorrer cada grupo de tratamiento. Como yo ya tengo creada la lista separada (lista_tratamientos), puedo hacer un bucle for para comprobar la normalidad de los datos de cada lista
  	for (grupo in names(lista_tratamientos)) {
    	cat("\n---", grupo, "---\n")
    	print(shapiro.test(lista_tratamientos[[grupo]][[variable]]))
  	}
}
#Como todos los valores tienen un p-value superior a 0.05, no se rechaza la normalidad. Es decir, los datos siguen una distribución normal
#Sabiendo que los datos tienen una distribución normal, puedo hacer la prueba de comparación de medias con una ANOVA (si no tuvieran esta distribución habría que usar un test no paramétrico, como Kruskal–Wallis)
#Vuelvo a crear el mismo bucle de antes
for (variable in vars) {
	cat("\n### Variable:", variable, "###\n")
	# ANOVA clásico, uso este porque hay que comparar más de 2 grupos, si fueran sólo 2 usaría una t de Student
  	modelo <- aov(as.formula(paste(variable, "~ Tratamiento")), data = archivo) #El ANOVA
  	print(summary(modelo))     # resultado global, me da información general de los resultados del ANOVA
	#En caso de que haya diferencias, hago un Tukey post-hoc para ver entre qué grupos están las diferencias
	cat("\n--- Post-hoc Tukey ---\n")
 	print(TukeyHSD(modelo))
}

#Con el Tukey, podemos observar como en Glucosa no hay ninguna diferencia significativa entre los grupos, pero en Presión sí, por ejemplo, entre Placebo y Fármaco B

# 12. Realiza un ANOVA sobre la glucosa para cada tratamiento. (1 pt)
# ANOVA para comparar los valores de glucosa entre los distintos tratamientos
anova_glucosa <- aov(Glucosa ~ Tratamiento, data = archivo)
# Mostrar el resumen del ANOVA
summary(anova_glucosa)
# Si el ANOVA es significativo (p < 0.05), hacemos un test post-hoc de Tukey, que es igual que en el ejercicio anterior, para poder ver cuales son los valores que son estadísticamente distintos
TukeyHSD(anova_glucosa)
#Como se puede observar, no hay ningún tratamiento que cause un cambio significativo en los niveles de glucosa, porque los p-value son > 0.05

