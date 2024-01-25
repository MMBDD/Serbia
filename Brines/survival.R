rm(list=ls())


#install.packages("Matrix")
#install.packages("lme4")

library(Matrix)
library(lme4)

#install.packages("Matrix", dependencies=TRUE)


#install.packages("tidyverse")
library(tidyverse)

library(survival)
file="DataBrines.csv"

X=read.csv(file, row.names=1, check.names = FALSE, sep=";")

row_names_df_to_remove<-c("2023-12-19_30C (3)", "2023-12-19_15C (1)", "11-12-2023_15C (4)", "11-12-2023_15C (3)", "20-22-2023_15C (3)")


X=(X[!(row.names(X) %in% row_names_df_to_remove),])
LargeData <- X %>%
  pivot_longer(cols = -c(Name, Date, Day, Temp), names_to = "Mass", values_to = "Intensity")





# Suponiendo que tus datos se llaman 'tus_datos'
LargeData$Event <- ifelse(LargeData$Day >= 1, 1, 0)

# Muestra los primeros registros para verificar la adición de la columna 'Evento'
head(LargeData)


# Ejemplo de modelo de Kaplan-Meier

modelo_supervivencia <- survfit(Surv(Day, Event) ~ Temp, data = LargeData)
summary(modelo_supervivencia)




# Suponiendo que tienes una columna llamada 'Temp' en tus datos que indica la temperatura
# Puedes tener diferentes colores para diferentes niveles de 'Temp'

# Crear un objeto de supervivencia para cada nivel de temperatura
surv_object_temp_15 <- Surv(time = LargeData$Day[LargeData$Temp == 15], event = LargeData$Event[LargeData$Temp == 15])
surv_object_temp_25 <- Surv(time = LargeData$Day[LargeData$Temp == 25], event = LargeData$Event[LargeData$Temp == 25])
surv_object_temp_30 <- Surv(time = LargeData$Day[LargeData$Temp == 30], event = LargeData$Event[LargeData$Temp == 30])

# Ajustar modelos de Kaplan-Meier para cada nivel de temperatura
modelo_supervivencia_temp_15 <- survfit(surv_object_temp_15 ~ 1, conf.type = "log-log")
modelo_supervivencia_temp_25 <- survfit(surv_object_temp_25 ~ 1, conf.type = "log-log")
modelo_supervivencia_temp_30 <- survfit(surv_object_temp_30 ~ 1, conf.type = "log-log")

# Dibujar la curva de supervivencia para la temperatura 15
plot(modelo_supervivencia_temp_15, col = "blue", lwd = 2, main = "Survival Curve", xlab = "Day", ylab = "Survival probability")

# Añadir la curva de supervivencia para la temperatura 25
lines(modelo_supervivencia_temp_25, col = "red", lwd = 2)

# Añadir la curva de supervivencia para la temperatura 30
lines(modelo_supervivencia_temp_30, col = "green", lwd = 2)

# Añadir leyenda
legend("topright", legend = c("Temp 15", "Temp 25", "Temp 30"), col = c("blue", "red", "green"), lwd = 2)


