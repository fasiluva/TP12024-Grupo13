setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(readxl)
library(ggplot2)
library(magrittr)
library(friends)
library(dplyr)
library(scales)
library(forcats)
library(stringr)

source("tp_functions.R")

Datos_LP <-
  read_excel("Datos_LP.xlsx", 
             skip = 2)

fuente <- "Fuente: La Poderosa. Encuesta realizada a barrios populares de Argentina."
title <- "Barrios populares de Argentina, año 2020."

### Para variable categórica nominal

# Tomo la columna el lugar que habitan actualmente es
# y hago una tabla con las columnas 
# `El lugar que habitan actualmente es`
# `frecuencia` (Se refiere a cuantas perosonas viven de tal forma,
#  alquilada, prestada, etc)
# `perc` Porcentaje de personas que viven en una casa, prestada, etc
lugar <- Datos_LP %>%
  group_by(`El lugar que habitan actualmente es:`) %>%
  summarize(frecuencia = n()) %>% # Ordena por frecuencia
  arrange(frecuencia) %>%
  mutate(`El lugar que habitan actualmente es:` = 
          fct_inorder(`El lugar que habitan actualmente es:`),
          frecRelativa = frecuencia/sum(frecuencia))

# Gráfico de barras horizontales de lugar, el cual indica
# cantidad de personas y porcentaje por situación de vivienda
# en la que viven las personas encuestadas
lugar %>%
  ggplot(aes(frecRelativa,
            `El lugar que habitan actualmente es:`))+
  geom_col(fill = "#458B74") +
  scale_y_discrete(labels = label_wrap(20)) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = paste("Situación de vivienda. ",title),
       y = "El lugar que habitan actualmente es",
       x = "Porcentaje de viviendas",
       caption = fuente)

moda(Datos_LP$`El lugar que habitan actualmente es:`)

### Para variable categórica ordinal

# Tomo la columna `¿Cómo es la presión del agua?`
# y hago una tabla con las columnas 
# `¿Cómo es la presión del agua?`
# `frecuencia` (Se refiere a cuantas personas tienen 
#  tal presión de agua)
# `perc` Porcentaje de personas que tienen tal presión
# de agua
presAgua <- Datos_LP %>%
  group_by(`¿Cómo es la presión del agua?`) %>%
  summarize(frecuencia = n()) %>%
  arrange(factor(`¿Cómo es la presión del agua?`,
                 levels = c("Buena", "Débil", "Muy débil"))) %>%
  mutate(frecRelativa = frecuencia/sum(frecuencia))

# Gráfico de sectores circulares para la variable presión de agua
presAgua %>%
  ggplot(aes(x = "", y = frecRelativa*100, fill = `¿Cómo es la presión del agua?`)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y") +
  geom_text(aes(label = paste0(round(frecRelativa*100), "%")),
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = paste("Presión de agua. ",title),
       caption = fuente) +
  guides(fill = guide_legend(reverse = FALSE)) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "black"))

# Gráfico de barras para la variable presión de agua
presAgua %>%
  ggplot(aes(`¿Cómo es la presión del agua?`,
             frecRelativa))+
  geom_col(width = 0.7,fill = "#458B74") +
  scale_x_discrete(limits = c("Muy débil", "Débil", "Buena")) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() + 
  labs(title = paste("Presión de agua. ",title),
       x = "Nivel de presión de agua",
       y = "Porcentaje de viviendas",
       caption = fuente)

### Para variable categórica con selección múltiple

# Filtro los datos de Datos_LP para tener en una variable sólo
# las columnas 
# `¿Cuál es la principal fuente de energía que utiliza para calefaccionar la vivienda?`,
# `...44`, `...45`, `...46`
princCalef <-
  Datos_LP[, c("¿Cuál es la principal fuente de energía que utiliza para calefaccionar la vivienda?",
               "...44", "...45", "...46", "...47", "...48")]

# Cambio los nombres de las columnas por nat, env, elec, 
# lencar, no tengo y noneces, dónde nat representa si la 
# persona utliza gas natural, env si utiliza gas envasado,
# lencar si utiliza leña o carbon si utiliza carbón,
# no tengo si no tiene para calefaccionar la vivienda y
# noneces si no necesita calefaccionar la vivienda
colnames(princCalef)  <- c("nat","env","elec","lencar", 
                           "notengo", "noneces")

# Obtengo una tabla con las frecuencias de cada columna de
# princCalef
princCal_df <- data.frame(
  type = c("Gas natural", "Gas envasado",
      "Electricidad","Leña/Carbon", "No tengo para calefaccionar",
      "No necesito calefaccionar"),
  freq = c(sum(grepl("Gas natural", princCalef$nat)),
           sum(grepl("Gas envasado", princCalef$env)),
           sum(grepl("Electricidad", princCalef$elec)),
           sum(grepl("/", princCalef$lencar)),
           sum(grepl("No tengo para calefaccionar", princCalef$notengo)),
           sum(grepl("No necesito calefaccionar", princCalef$noneces)))
  ) %>% arrange(freq) %>%
  mutate(type = fct_inorder(type))

# Creo una gráfico de barras horizontales, ordenadas 
# desde la mayor, hasta la menor frecuencia con los
# datos que están en pricCal_df
princCal_df %>%
  ggplot(aes(freq, type)) +
  geom_col(fill = "#458B74") +
  theme_bw() +
  scale_x_continuous(breaks = c(100,200,300,400,500)) +
  scale_y_discrete(labels = label_wrap(20)) +
  labs(y = "Fuente de calefacción",
       x = "Cantidad de viviendas",
       title = paste("Principales fuentes de calefacción en barrios populares.\n", 
                     title, sep = ""),
       caption = fuente)



### Para variable cuantitativa discreta

# Creo una tabla con los datos de la columna
# `¿Cuántos abonos/prepagos de datos móviles sostienen por mes en la vivienda?`
# el cual indica la frecuencia y el porcentaje de las viviendas
# que tienen la cantidad de abonos dada
abonos <- Datos_LP %>%
  group_by(`¿Cuántos abonos/prepagos de datos móviles sostienen por mes en la vivienda?`) %>%
  summarize(frecuencia = n()) %>%
  mutate(frecAcum = cumsum(frecuencia), frecRel = frecuencia/sum(frecuencia),
         cumsum(frecuencia/sum(frecuencia)))
colnames(abonos) <- c("abonosxcasa", "frecuencia", "frecAcum", "frecRelativa","frecRelativaAcum")

# Creo un gráfico de bastones para la información en abonos
abonos %>%
  ggplot(aes(abonosxcasa,
             frecRelativa))+
  geom_col(fill = "#458B74", width = 0.15) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(x = "abonos/prepagos de datos móviles",
       y = "Porcentaje de viviendas",
       title = paste("Distribución de abonos/prepagos de datos móviles que se sostiene\npor mes en las viviendas\n",
                     title, sep=""),
       caption = fuente)

abonoslist <- as.integer(Datos_LP$`¿Cuántos abonos/prepagos de datos móviles sostienen por mes en la vivienda?`)
IQR(abonoslist)
median(abonoslist)

### Para variable cuantitativa continua

# Obtenga de Datos_LP sólo la columna
# `Porcentaje de aumento en el último año`
# y omito los campos que no están completos
# (los compos con NA)
percAumentos <- Datos_LP[, "Porcentaje de aumento en el último año"] %>%
  na.omit()
colnames(percAumentos) <- c("aumentos_perc")

# Creo un histograma para percAumentos
percAumentos %>%
    ggplot(aes(x=aumentos_perc)) +
    geom_histogram(binwidth = 0.5,fill = "#458B74", bins = 80,
                   alpha = 0.8, boundary = 0) +
    theme_bw() +
    scale_x_continuous(labels = scales :: percent) +
    labs(x = "Porcentaje de aumentos",
         y = "Cantidad de viviendas",
         title = paste("Porcentaje de aumento de alquiler en el último año\n",
                       title, sep=""),
         caption = fuente)

# Mediana
median(percAumentos$aumentos_perc)

#Rango intercurtilico
IQR(percAumentos$aumentos_perc)

### Relación entre una variable categórica y una variable
### cuantitativa

# Obtengo sólo la columna de tiempo de residencia de Datos_LP
tiempoReside <- as.data.frame(sapply(
  Datos_LP[, "Tiempo de residencia en la vivienda actual (en años)"],
  as.numeric))
# Cambio el nombre de la columna a tiempo
colnames(tiempoReside) <- c("tiempo")

# Obtengo sólo la columna de el lugar que habitan actualmente de
# Datos_LP
estadoLugar <- Datos_LP[, "El lugar que habitan actualmente es:"]
# Cambio el nombre de la columna a estado
colnames(estadoLugar) <- c("estado")

# Boxplot comparativo
# Posesión del lugar que habitan con respecto a los
# años viviendo en ese lugar
Datos_LP %>%
  ggplot(aes(estadoLugar$estado, tiempoReside$tiempo))+
  geom_boxplot(fill = "#EEC591") +
  scale_x_discrete(labels = label_wrap(20)) +
  scale_y_continuous(breaks = seq(0,150,10)) +
  theme_bw() +
  labs(x = "El lugar que habitan actualmente es",
       y = "Cantidad de años en la vivienda",
       title = paste("Tipo posesión de vivienda con respecto a los años que viven en la misma\n",
                     title, sep = ""),
       caption = fuente)


formaObtAgua <- Datos_LP$`¿De qué forma obtiene el agua dentro de su vivienda?`
formaObtAgua <- str_replace_all(formaObtAgua, 
  "A través de una conexión sin medidor, es decir “informalmente”, sea a través de una conexión directa a la red pública o a través de una conexión indirecta a través de un vecinx “informalmente”",
  "A través de una conexión sin medidor, es decir “informalmente”")

### Relación entre dos variables categóricas

# Graáfico de barras subdivididas para las columnas de
# Datos_LP
# `¿De qué forma obtiene el agua dentro de su vivienda?`
# `¿Cómo es la presión del agua?`
Datos_LP %>%
  ggplot(aes(formaObtAgua,
    fill = `¿Cómo es la presión del agua?`))+
  geom_bar(position = "dodge", alpha = 0.5) +
  scale_x_discrete(labels = label_wrap(20)) +
  theme_bw() +
  labs(x = "Forma de obtener el agua",
       y = "Cantidad de viviendas",
       title = paste("Presión del agua respecto a la forma en que se obtiene\n",
                     title, sep=""),
       caption = fuente)


### Relación entre dos variables cuantitativas

# Obtengo sólo la columna de `Edad jefe/a del hogar`
# de Datos_LP
edadJefeH <- as.data.frame(sapply(
  Datos_LP[, "Edad jefe/a del hogar"],
  as.numeric))
# Cambio el nombre de la columna a c1
colnames(edadJefeH) <- c("c1")

# Diagrama de dispersión para variables
# edad jefe hogar y tiempo de residencia
# La variable es edad jefe hogar y la de
# respuesta tiempo de residencia
Datos_LP %>%
  ggplot(aes(x = edadJefeH$c1,
             y = tiempoReside$tiempo)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0,100,10)) +
  scale_y_continuous(breaks = seq(0,110,10)) +
  labs(x = "Edad jefe hogar",
       y = "Tiempo de residencia",
       title = paste("Edad jefe hogar con respecto al tiempo de residencia\n",
                     title, sep = ""))

