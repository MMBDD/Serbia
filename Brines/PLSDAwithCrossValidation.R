library(mdatools)
library(rchemo)
#Import Data
file="GTdata.csv"
data=read.csv(file, row.names=1, check.names = FALSE, sep=";")


#Remove outliers

row_names_df_to_remove<-c("2023-12-19_30C (3)", "2023-12-19_15C (1)", "11-12-2023_15C (4)", "11-12-2023_15C (3)", "20-22-2023_15C (3)")


data=(data[!(row.names(data) %in% row_names_df_to_remove),])

X=data[,5:1200]
y=data[,4]



###############################
######## Again from the beginning


# Install and load necessary packages
install.packages(c("pls", "caret"))
library(pls)
library(caret)

# Assume your data sets are called 'X' for predictors and 'y' for responses

# Set up PLSDA model
plsdamodel <- pls::plsr(y ~ ., data = cbind(X, y), ncomp = 2)  # Adjust 'ncomp' as needed

# Combine 'X' and 'y' into one data frame
my_data <- cbind(X, y)
my_data$y <- as.factor(my_data$y)

# Create training control for cross-validation
ctrl <- trainControl(method = "cv", number = 5)  # 5-fold cross-validation, adjust as needed

# Train the PLSDA model with cross-validation
plsdacv <- train(y ~ ., data = my_data, method = "pls", trControl = ctrl, tuneLength = 10)  # Adjust 'tuneLength' as needed

# Make predictions on the test set
predictions <- as.data.frame(predict(plsdacv, newdata = X))
 g=as.factor(predictions$`predict(plsdacv, newdata = X)`)
y=as.factor(y)
# Confusion matrix
conf_matrix <- confusionMatrix(g, y)

# Extract performance metrics
accuracy <- conf_matrix$overall["Accuracy"]
balanced_accuracy <- conf_matrix$byClass["Balanced Accuracy"]
sensitivity <- conf_matrix$byClass["Sensitivity"]
specificity <- conf_matrix$byClass["Specificity"]
precision <- conf_matrix$byClass["Positive Predictive Value"]

# Print the results
cat("Accuracy:", accuracy, "\n")
cat("Balanced Accuracy:", balanced_accuracy, "\n")
cat("Sensitivity:", sensitivity, "\n")
cat("Specificity:", specificity, "\n")
cat("Precision:", precision, "\n")

# Display the confusion matrix
print(conf_matrix)
