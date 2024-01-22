install.packages("ggplot2")
library(ggplot2)
library(dplyr)

file="DataBrines.csv"

X=read.csv(file, row.names=1, check.names = FALSE, sep=";")

X=X[,5:1199]

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


##Explore by PCA to remove outliers at sepal level

pca=prcomp(X, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL, rank.=20)
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
