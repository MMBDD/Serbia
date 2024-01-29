
install.packages(c("pls", "pROC", "caret"))

library(rchemo)
library(mdatools)

#Import Data
file="GTdata.csv"
data=read.csv(file, row.names=1, check.names = FALSE, sep=";")


#Remove outliers

row_names_df_to_remove<-c("2023-12-19_30C (3)", "2023-12-19_15C (1)", "11-12-2023_15C (4)", "11-12-2023_15C (3)", "20-22-2023_15C (3)")


data=(data[!(row.names(data) %in% row_names_df_to_remove),])

# Split data set per Period



# Specify the periods manually
data$Period <- ifelse(data$Day <= 8, "Period 1",
                      ifelse(data$Day <= 15, "Period 2",
                             ifelse(data$Day <= 20, "Period 3",
                                    ifelse(data$Day <= 27, "Period 4", "Period 5"))))

# Split the data frame by the "Period" column
split_data <- split(data, data$Period)

period_dataframes <- list()

# Convert each element of the list to a data frame
for (i in 1:length(split_data)) {
  period_dataframes[[i]] <- as.data.frame(split_data[[i]])
}


# Create separate data frames with unique names
for (i in seq_along(period_dataframes)) {
  assign(paste0("df_Period", i), period_dataframes[[i]])
}



y1=as.data.frame(df_Period1$Fermentation)
y2=as.data.frame(df_Period2$Fermentation)
y3=df_Period3$Fermentation
y4=df_Period4$Fermentation
y5=df_Period5$Fermentation

x1=df_Period1[,5:1200]
x2=df_Period2[,5:1200]
x3=df_Period3[,5:1200]
x4=df_Period4[,5:1200]
x5=df_Period5[,5:1200]

x3<- t(apply(x3, 1, function(row) row - mean(row)))
x4<- t(apply(x4, 1, function(row) row - mean(row)))
x5<- t(apply(x5, 1, function(row) row - mean(row)))

# Mostrar el resultado


# PLSDA model with CV
y3 <- as.factor(y3)
# Supongamos que 'y' es tu vector original
# Cambiar todos los 0 a 1 y todos los 1 a 2
y3 <- ifelse(y3 == 0, 1, ifelse(y3 == 1, 2, y3))

# Instala el paquete rchemo si no lo has hecho
install.packages("rchemo")

# Carga la biblioteca
library(rchemo)
#Optimal number of latent variables; (2 for x3, 3 for x4, 2 for x5)
fmD <- plsrda(x5,  y5, nlv = 2)

pred <- predict(fmD, x5)$pred


mx=as.matrix(y5)


###### Now we can see the Results of the Predictions

class(pred)
tab1=table(pred, mx)

tab1
###### Results in Validation


re1=confusionMatrix(tab1)

re1

###### Confusion matrix

as.table(re1)

as.matrix(re1)

###### Classification parameters

t=as.matrix(re1, what="overall")
t
e=as.matrix(re1, what="classes")
e

###### Summary of classification parameters

ComparativeTable=as.data.frame(c(t[1], e[1], e[2], e[5], e[11]))

rownames(ComparativeTable) <- c("Accuracy", "Sensitivity", "Specificity", "Precision", "Balanced accuracy") 

ComparativeTable



######################### ANOTHER PACKAGE TO COMPARE


y3 = factor(y3, labels = c("Normal", "Unnormal"))
y4 = factor(y4, labels = c("Normal", "Unnormal"))
y5 = factor(y5, labels = c("Normal", "Unnormal"))

# PLSDA
m = plsda(x5, y5, ncomp = 20, cv = 1)
summary(m)
plot(m)

# Confusion Matrix
getConfusionMatrix(m$calres)  

par(mfrow = c(1, 2))
plotPredictions(m)


par(mfrow = c(3, 2))
plotMisclassified(m, nc = 2)

plotSensitivity(m, nc = 2)

plotSpecificity(m, nc = 2)


par(mfrow = c(1, 2))
plotRegcoeffs(m, ncomp = 3)



#Show results in a report
#Understand what happens in period 5








