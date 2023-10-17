library(prospectr)
library(caret)
library(rchemo)

#The aim of this code is to create optimized PLSDA models
# To that end, different pretreatments to raw spectra will be tested
# Fuerthermore, models will be created with different number of important variables as input for PLSDA
# You can use this code to create Intravariety, Intervariety and Global Models
#Global Model, combining two or three varieties together.
#Set working directory

setwd("C:/temp/ICNIRS/ICNIRS/'GlobalModel")
# Import datasets, each data set should have only one spectrum per sepal
# First use the Code  "OutlierDetection.R", to build averaged data sets for each varieties

file1= "BRIOSOready.csv"
file2="Labelsbrioso.csv"

file3 = "Imagen1y2.csv"
file4="YlabelsCappriccia2.csv"

file5 = "XProvine.csv"
file6="LabelsProvine.csv"

# Import Table
XB=read.csv(file1, row.names=1, check.names=FALSE)

# Select only Spectra Columns 
XB=XB[,4:115]

# Go down to line 117 and perform Exploration by PCA
# Detect outliers and add them to row_names_df_to_remove1
row_names_df_to_remove1<-c("2599", "4038", "983", "3222", "3955", "3651", "440", "5160")

#Import Labels
YB=read.csv(file2, row.names=1)

#Select column related to Label Scenario 2
YB=as.data.frame(YB[,12]) 

# Before removing outliers, join X and Y in order to remove the same samples from both datasets
joinB=cbind(XB,YB)

#Remove outliers at sepal level
joinB=(joinB[!(row.names(joinB) %in% row_names_df_to_remove1),])

#Create cleaned X and Y matrices
XB=joinB[,1:112]
YB=as.data.frame(joinB[,113])

#Repeat the same procedure with the variety Cappricia 
XC=read.csv(file3, row.names=1, check.names=FALSE)
row_names_df_to_remove<-c("i2T11S5", "i2T5S3", "i1T1S3", "i1T1S2")
row_names_df_to_remover<-c("i1T1S3", "i1T1S2")
YC=read.csv(file4, row.names=1)
YC=as.data.frame(YC[,12]) 
joinC=cbind(XC,YC)
joinC=as.data.frame(joinC[!(row.names(joinC) %in% row_names_df_to_remove),])
joinC=as.data.frame(joinC[!(row.names(joinC) %in% row_names_df_to_remover),])
XC=joinC[,1:112]
YC=as.data.frame(joinC[,113])

#Repeat the same procedure with the variety Provine
XP=read.csv(file5, row.names=1, check.names=FALSE)
XP=XP[,4:115]
YP=read.csv(file6, row.names=1)
YP=as.data.frame(YP[,11]) 
joinP=cbind(XP,YP)
row_names_df_to_remove2<-c("7790", "9639", "0", "7232","10340", "1002", "4117")
joinP=joinP[!(row.names(joinP) %in% row_names_df_to_remove2),]
XP=joinP[,1:112]
YP=as.data.frame(joinP[,113])



##Explore by PCA to remove outliers at sepal level

pca=prcomp(XB, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL, rank.=20)
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


#Standarize Colnames

colnames(XC)=colnames(XP)
colnames(YC)=colnames(YP)
colnames(XC)=colnames(XB)
colnames(XP)= colnames(XB)
colnames(YB)="Class"
colnames(YC)="Class"
colnames(YP)="Class"

# Build one data set
# You can change the data sets to try a different global model

X=rbind(XC, XP)
Y=rbind(YC,YP)

#Apply pretreatment to raw spectra
# Standard Normal Variate

snvC=as.data.frame(snv(X))

#First or Second Derivative, to SNV or to raw data
SGC=as.data.frame(savgol(snvC,1,17,2))

#Detrend
wav <- as.numeric(colnames(X))
DT=as.data.frame(detrend(X, wav, 2)) 

## Plot spectra

# You can plot raw of pretreated spectra

par(mar = c(6, 6, 6, 6))
col = "red"
plotsp(snvC, ylab = "Measured Reflectance",xlab = "Wavelength (nm)",type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE)
title(main = " SNV NIR spectra",line = NA, outer = FALSE)

#Join X and Y before splitting dataset
# This is the moment where you can choose raw or pretreated spectra
# After you run the code until the end, it is advisable to return to this line and test another pretreatment method

join=cbind(Y,SGC)
row_names_df_to_remove<-c("2599")
join=join[!(row.names(join) %in% row_names_df_to_remove),]
nv=112

#Split dataset in a representative way for each class 70% for calibration, and 30% for validation

JoinSepals<- split(join, join$Class)
xtrain = data.frame()
ytrain =as.integer()
xtest = data.frame()
ytest=as.integer()

for (class in JoinSepals)
{
  dt = sort(sample(nrow(class), nrow(class)*.7))
  train<-class[dt,]
  test<-class[-dt,]
  tmpxtrain = train[,2:(nv+1)]
  xtrain= rbind(xtrain,tmpxtrain)
  tmpytrain =  train[,1]
  ytrain=append (ytrain, tmpytrain)
  tmpxtest = test[,2:(nv+1)]
  xtest =  rbind(xtest,tmpxtest)
  tmpytest =  test[,1]
  ytest = append(ytest,tmpytest)
}

# Check the dimensionality of the results

dim(xtrain)
dim(as.data.frame(ytrain))
dim(xtest)
dim(as.data.frame(ytest))

# Select different number of variables as input for PLSDA
# Select the best model according to classification errors in the validation set
# We are going to study from 5 to 39 variables selected by CovSel

## COVSEL 
i =5
to = 40
results = data.frame()
results = data.frame(matrix(nrow = 5, ncol = 0)) 

while( i < to)
{
  nvar = i
  #Covsel algorithm 
  vs=(covsel(xtrain, ytrain, nvar, scaly = TRUE, weights = NULL))
  variables=vs$sel[1]
  variables[1]
  
  ## columns subset
  v <- unlist(variables)
  v = sort(v)
  
  # Select the important variables from the original matrices
  CovSelTrain =xtrain[,c(v)]
  CovSelTest =xtest[,c(v)]
  
  ## Split TRAINING SET into VALIDATION AND TUNING in order to train the PLSDA model
  # This is important in order to select the optimal number of Latent Variables (LVs)
  OneTable= cbind(CovSelTrain,ytrain)
  noSplit = TRUE 
  while (noSplit)
  { 
    indices= sample(223, 40)
    validation =  OneTable[indices,]
    tuning  =  OneTable[-indices,]   
    Yvalidation=validation[,(nvar+1)]                    
    Xvalidation=validation[,-(nvar+1)] 
    Ytuning=tuning[,(nvar+1)]                    
    Xtuning=tuning[,-(nvar+1)] 
    ## Train a PLSDA model 
    # Calibrate a model in Validation set, Tune the model with the Tuning Set
    #m <- 50 ; p <- 8
    
    nlv <- nvar-1
    res <- gridscorelv(Xvalidation, Yvalidation, Xtuning, Ytuning, score = err,fun = plsrda,nlv = 1:nlv, verb = TRUE)
    
    # Plot to observe the latent variables with lower classification error (Accuracy)
    plotscore(res$nlv, res$y1,main = "ERR", xlab = "Nb. LVs", ylab = "Value")
    u<- res[res$y1 == min(res$y1), ][1, , drop = FALSE]
    ytrain=as.data.frame(ytrain)
    YvalInteger = as.integer(ytrain[,1])
  
    # Calibrate a PLSDA model in the Training Set
    fmD <- plsrda(CovSelTrain,  YvalInteger, nlv = u$nlv)
    
    # Make predictions in the Test Set
    pred <- predict(fmD, CovSelTest)$pred
    noSplit = (var(pred) == 0)
  }
  
  mx=as.matrix(ytest)
  ###### Results       
  class(pred)
  tab1=table(pred, mx)
  ######
  re1=confusionMatrix(tab1)
  ######
  as.table(re1)
  as.matrix(re1)
  ######
  t=as.matrix(re1, what="overall")
  e=as.matrix(re1, what="classes")
  ######
  comparativeTable=as.data.frame(c(t[1], e[1], e[2], e[5], e[11]))
  rownames(comparativeTable) <- c("Accuracy", "Sensitivity", "Specificity", "Precision", "Balanced accuracy") 
  colnames(comparativeTable)<- (nvar)
  results = cbind(results,comparativeTable)
  results
  i = i+1
}

# Compare results with different number of imporant variables as input for PLSDA
results

# Export the results for the selected datasets and pretreatments
# Please, add the pretreatment used in the name of the file

write.csv(results, "c:/temp/OPT/resultsNN.csv")


# Change pretreatments and run the code again

# Change data sets and run the code again


























































