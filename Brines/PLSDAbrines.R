
library(mdatools)
library(rchemo)
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



y1=df_Period1$Fermentation
y2=df_Period2$Fermentation
y3=df_Period3$Fermentation
y4=df_Period4$Fermentation
y5=df_Period5$Fermentation

x1=df_Period1[,5:1200]
x2=df_Period2[,5:1200]
x3=df_Period3[,5:1200]
x4=df_Period4[,5:1200]
x5=df_Period5[,5:1200]


# Divide each peak by the total sum of intensities (Normalization 1)
normalize_by_sum <- function(data) {
  normalized_data <- data / rowSums(data)
  return(normalized_data)
}

# Apply normalization 1
Normalized.x2.1 <- normalize_by_sum(x2)
Normalized.x3.1 <- normalize_by_sum(x3)
Normalized.x4.1 <- normalize_by_sum(x4)
Normalized.x5.1 <- normalize_by_sum(x5)

#Substract the average and divide by standard deviation (Normalization 2)

#Substract the average-Normalizarion 2

x2.a<- t(apply(x2, 1, function(row) row - mean(row)))
x3.a<- t(apply(x3, 1, function(row) row - mean(row)))
x4.a<- t(apply(x4, 1, function(row) row - mean(row)))
x5.a<- t(apply(x5, 1, function(row) row - mean(row)))

# Divide by Standard deviation- NOrmalization 2


scale_by_row_sd <- function(data) {
  scaled_data <- t(apply(t(data), 1, function(row) row / sd(row)))
  return(scaled_data)
}

# Apply Normalization 2
Normalized.x2.2 <- as.data.frame(t(scale_by_row_sd(x2.a)))
Normalized.x3.2 <- as.data.frame(t(scale_by_row_sd(x3.a)))
Normalized.x4.2 <- as.data.frame(t(scale_by_row_sd(x4.a)))
Normalized.x5.2 <- as.data.frame(t(scale_by_row_sd(x5.a)))
############################################################################################################################


#######################################################
#Cross validation

y2 <- unlist(y2)

# Max number of components to consider in the model
max_nlv <- 9

# Cross Validation results
cv_results <- matrix(NA, nrow = max_nlv, ncol = 2, 
                     dimnames = list(NULL, c("Error Rate", "Accuracy")))

# Implementación de LOOCV
for (i in 1:max_nlv) {
  # Entrenar modelo PLS-DA con i componentes
  fmD <- plsrda(Normalized.x5.1, y5, nlv = i)
  
  # Make predictions and evaluate performances
  predictions <- predict(fmD, Normalized.x5.1)
  confusion_matrix <- table(predictions$pred, y5)
  
  # Calculate error and prediction
  error_rate <- 1 - sum(diag(confusion_matrix)) / sum(confusion_matrix)
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  # Save results in the matrix
  cv_results[i, ] <- c(error_rate, accuracy)
}

# Show results
print("Cross validation results:")
print(cv_results)
print(confusion_matrix)




###############################################
# PLSDA model with CV

fmD <- plsrda(Normalized.x2.1,  y2, nlv = 7)

pred <- predict(fmD, Normalized.x2.1)$pred


mx=as.matrix(y2)

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

#####################################################################################################################

