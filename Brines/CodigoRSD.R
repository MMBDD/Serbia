install.packages("ggplot2")
library(ggplot2)
library(dplyr)
library (rchemo)

file="DataBrines.csv"

X=read.csv(file, row.names=1, check.names = FALSE, sep=";")

row_names_df_to_remove<-c("2023-12-19_30C (3)", "2023-12-19_15C (1)", "11-12-2023_15C (4)", "11-12-2023_15C (3)", "20-22-2023_15C (3)")


X=(X[!(row.names(X) %in% row_names_df_to_remove),])


# Prepare the file for PCA analysis


xPCA=X[,5:1199]

# Perform PCA to remove outliers and identify important masses


pca=prcomp(xPCA, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL, rank.=20)
summary(pca)

##Calculate scores and loadings
pcascores <- pca$x
pcaloadings <- pca$rotation

##Plot PCA scores
# Visualize the plots to detect extreme outliers
# "sample 2646" in dataset XB can be seen as an oulier in plot PC7 vs PC8
# Please, go back to line 81 and add this sample to row_names_df_to_remove1

PCAcolors <- c("#66c2a5","#fc8d62","#8da0cb", "#E5C494", "#B3B3B3")
plot(pcascores[,1:2], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,1:2],labels=rownames(pcascores))

plot(pcascores[,2:3], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,2:3],labels=rownames(pcascores))

plot(pcascores[,c(1,3)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(1,3)],labels=rownames(pcascores))

plot(pcascores[,4:3], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,4:3],labels=rownames(pcascores))

plot(pcascores[,c(4,5)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(4,5)],labels=rownames(pcascores))

plot(pcascores[,c(6,5)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(6,5)],labels=rownames(pcascores))

plot(pcascores[,c(7,6)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(7,6)],labels=rownames(pcascores))

plot(pcascores[,c(7,8)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(7,8)],labels=rownames(pcascores))


PCAcolors <- c("#66c2a5","#fc8d62","#8da0cb", "#E5C494", "#B3B3B3")
plot(pca$rotation[,2], pch=18, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pca$rotation[,2],labels=rownames(pca$rotation))

# You can plot raw of pretreated spectra

par(mar = c(6, 6, 6, 6))
col = "red"
plotsp(pca$rotation[,2], ylab = "PC1 Loadings",xlab = "Masses",type = "l", col = col, zeroes = FALSE, labels = rownames(pca$rotation), add = FALSE)
title(main = "PC1 Loadings",line = NA, outer = FALSE)

# Detect important masses


# Change colour of bar plot
c.pc1 <- ifelse(pca$rotation[,1] > 0, yes="green2", no="red2")
c.pc2 <- ifelse(pca$rotation[,2] > 0, "green2", "red2")
c.pc3 <- ifelse(pca$rotation[,3] > 0, "green2", "red2")


n.pc1 <- ifelse(pca$rotation[,1] > 0, yes=-0.05, no=pca$rotation[,1]-0.01)
n.pc2 <- ifelse(pca$rotation[,2] > 0, yes=-0.0, no=pca$rotation[,2]-0.01)
n.pc3 <- ifelse(pca$rotation[,3] > 0, yes=-0.0, no=pca$rotation[,3]-0.01)



par(mar=c(5,5,5,2)) # Set margins
par(mfrow = c(2, 1))
barplot(pca$rotation[,1], main="PC 1",cex.main=4, cex.axis=3, ylim = c(-0.03, 0.03), las=2, col=c.pc1,cex.sub=2, cex.lab=2, cex.names=2)
barplot(pca$rotation[,2], main="PC 2",cex.main=4, cex.axis=3, ylim = c(-0.03, 0.03), las=2, col=c.pc1,cex.sub=2, cex.lab=2, cex.names=2)
barplot(pca$rotation[,3], main="PC 3",cex.main=4, cex.axis=3, ylim = c(-0.03, 0.03), las=2, col=c.pc1,cex.sub=2, cex.lab=2, cex.names=2)




#Without Mass 45
X1=X[39:152,]

#Only Mass 45
X2=X[1:38,]

#Only Masses 57 and 74
X6=X[c(39:76,115:152),]


model6 <- aov(Intensity ~ Temp*Mass*Day, data=X6)
summary(model6)


coef(model6)

ggplot(X6, aes(x=Day, y=Intensity, color=factor(Mass), linetype=factor(Temp))) +
  geom_line(size = 1.5) +
  labs(title="Day vs Intensity",
       x="Day",
       y="Intensity") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))


model <- aov(Intensity ~ Temp*Mass*Day, data=X)
summary(model)


coef(model)

ggplot(X, aes(x=Day, y=Intensity, color=factor(Mass), linetype=factor(Temp))) +
  geom_line(size = 1.5) +
  labs(title="Day vs Intensity",
       x="Day",
       y="Intensity") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))


boxplot(Intensity ~ Temp, data = X)
labs(x="Temperature °C",y="Intensity")




model1 <- aov(Intensity ~ Temp*Mass*Day, data=X1)
summary(model1)


coef(model1)

ggplot(X1, aes(x=Day, y=Intensity, color=factor(Mass), linetype=factor(Temp))) +
  geom_line(size = 1.5) +
  labs(title="Day vs Intensity",
       x="Day",
       y="Intensity")

boxplot(Intensity ~ Temp, data = X1)
labs(x="Temperature °C",y="Intensity")




model2 <- aov(Intensity ~ Temp*Mass*Day, data=X2)
summary(model2)


coef(model2)

ggplot(X2, aes(x=Day, y=Intensity, color=factor(Temp), linetype=factor(Temp))) +
  geom_line(size = 1.5) +
  labs(title="Mass 45: Day vs Intensity",
       x="Day",
       y="Intensity")

boxplot(Intensity ~ Temp, data = X2)
labs(x="Temperature °C",y="Intensity")

X18=X[X$Day == 18, ]
modelDay18 <- aov(Intensity ~ Temp*Mass, data=X18)
summary(modelDay18)
coef(modelDay18)


X15=X[X$Day == 1,]
modelDay15 <- aov(Intensity ~ Temp, data=X15)
summary(modelDay15)
coef(modelDay15)

######################## Important masses were not found, so another method was applied
rm(list=ls())


install.packages("Matrix")
install.packages("lme4")

library(Matrix)
library(lme4)

install.packages("Matrix", dependencies=TRUE)


install.packages("tidyverse")
library(tidyverse)
file="DataBrines.csv"

X=read.csv(file, row.names=1, check.names = FALSE, sep=";")

row_names_df_to_remove<-c("2023-12-19_30C (3)", "2023-12-19_15C (1)", "11-12-2023_15C (4)", "11-12-2023_15C (3)", "20-22-2023_15C (3)")


X=(X[!(row.names(X) %in% row_names_df_to_remove),])
LargeData <- X %>%
  pivot_longer(cols = -c(Name, Date, Day, Temp), names_to = "Mass", values_to = "Intensity")

#Fit a mixed lineal model
#mixedModel <- lmer(Intensity ~ Temp * Day + (1|Mass), data = LargeData)
# Instala el paquete nlme si aún no lo tienes instalado
install.packages("nlme")

# Carga el paquete nlme
library(nlme)

# Ajusta un modelo lineal mixto
model_mix_nlme <- lme(Intensity ~ Temp * Day, random = ~1 | Mass, data = LargeData)

# Muestra un resumen del modelo
summary(model_mix_nlme)


# Gráfico de Intensidad a lo largo del Tiempo
library(ggplot2)

ggplot(LargeData, aes(x = Day, y = Intensity, color = as.factor(Temp))) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Intensity changes with Time",
       x = "Day",
       y = "Intensity",
       color = "Temp")+
  theme(
    text = element_text(size = 12),  
    title = element_text(size = 16),  
    axis.title = element_text(size = 20),  
    legend.text = element_text(size = 20)  )

# Boxplot por Día
ggplot(LargeData, aes(x = factor(Day), y = Intensity, fill = as.factor(Temp))) +
  geom_boxplot() +
  labs(title = "Intensity variation per day",
       x = "Day",
       y = "Intensity",
       fill = "Tempe")

library(ggplot2)

ggplot(LargeData, aes(x = factor(Day), y = Intensity, fill = as.factor(Temp))) +
  geom_boxplot(outlier.shape = NA) +  # Elimina los outliers para que no afecten a la escala
  geom_jitter(position = position_jitter(width = 0.3, height = 0), size = 3, alpha = 0.7, color = "black") +  # Puntos agrandados
  labs(title = "Intensity variation per day",
       x = "Day",
       y = "Intensity",
       fill = "Temp") +
  theme(
    text = element_text(size = 14),  # Tamaño del texto general
    title = element_text(size = 18),  # Tamaño del título
    axis.title = element_text(size = 16),  # Tamaño de las etiquetas de ejes
    legend.text = element_text(size = 12)  # Tamaño del texto de la leyenda
  ) +
  scale_fill_manual(values = c("blue", "red", "green"))  # Colores personalizados para las temperaturas

library(ggplot2)

ggplot(LargeData, aes(x = factor(Day), y = Intensity, fill = as.factor(Temp))) +
  geom_boxplot(outlier.shape = NA, fill = c("blue", "red", "green")) +  # Relleno de los boxplots
  geom_jitter(position = position_jitter(width = 0.3, height = 0), size = 3, alpha = 0.7, color = "black") +  # Puntos agrandados
  labs(title = "Intensity variation per day",
       x = "Day",
       y = "Intensity",
       fill = "Temp") +
  theme(
    text = element_text(size = 14),  # Tamaño del texto general
    title = element_text(size = 18),  # Tamaño del título
    axis.title = element_text(size = 16),  # Tamaño de las etiquetas de ejes
    legend.text = element_text(size = 12)  # Tamaño del texto de la leyenda
  )
modelo_no_lineal <- nls(Intensity ~ b0 + b1 * Day + b2 * Day^2, data = LargeData,
                        start = list(b0 = 0, b1 = 0, b2 = 0))
summary(modelo_no_lineal)

coef(modelo_no_lineal)

# Coeficcients of the quadratic model
b0 <- 6.098e-04
b1 <- 8.894e-06
b2 <- -2.776e-07

# Calculation of the vertex of the quadratic parabola.
vertex_day <- -b1 / (2 * b2)
vertex_intensity <- b0 + b1 * vertex_day + b2 * vertex_day^2


cat("Day at which the intensity reaches its minimum or maximum:", round(vertex_day, 2), "\n")


cat("Intensity on that day:", round(vertex_intensity, 6), "\n")




# Fit a mixed-effects model with fixed effects for Temperature and Day, and a random effect for Mass
model_mixed_lme <- lme(Intensity ~ Temp * Day, random = ~1 | Mass, data = LargeData)

# Model summary
summary(model_mixed_lme)

