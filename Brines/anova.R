



file2="DataPerMass.csv"

X2=read.csv(file2, row.names=1, check.names = FALSE, sep=";")


#Only Masses 57 and 74
X6=X2[c(39:76,115:152),]



#Without Mass 45
X1=X2[39:152,]

#Only Mass 45
X4=X2[1:38,]


model4 <- aov(Intensity ~ Temp*Mass*Day, data=X4)
summary(model4)
coefficients(model4)





model6 <- aov(Intensity ~ Temp*Mass*Day, data=X6)
summary(model6)

# Instala el paquete ggplot2 si aún no lo tienes instalado
# install.packages("ggplot2")

# Carga el paquete ggplot2
library(ggplot2)

# Verificar la instalación y cargar paquetes
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)


library(emmeans)
library(ggplot2)

# Obtener medias ajustadas y intervalos de confianza
means <- emmeans(model6, ~ Temp * Mass * Day)

# Graficar los resultados
plot(means, by = c("Temp", "Mass"))

library(ggplot2)

# Asumiendo que 'datos' es tu conjunto de datos
# Ajustar el modelo mixto
model6 <- aov(Intensity ~ Temp * Mass * Day, data = X6)

# Obtener las medias ajustadas
means <- emmeans(model6, ~ Temp * Day)

# Convertir a formato de datos
means_data <- as.data.frame(means)

# Crear el gráfico
ggplot(means_data, aes(x = Day, y = emmean, color = as.factor(Temp))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = as.factor(Temp)), alpha = 0.2) +
  labs(title = "Influencia de la Temperatura a través de los Días",
       x = "Día",
       y = "Intensidad") +
  scale_fill_manual(values = c("blue", "red", "green")) +
  scale_color_manual(values = c("blue", "red", "green")) +
  facet_wrap(~Temp, scales = "free_y") +
  theme_minimal()


