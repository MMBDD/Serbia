library(prospectr)
library(caret)
library(rchemo)


#Global Model with all varieties together


file1= "BRIOSOa.csv"
file2="BRIOSOb.csv"

file3 = "CAPPRICIAa.csv"
file4="CAPPRICIAb.csv"

file5 = "PROVINEa.csv"
file6="PROVINEb.csv"


# Import Brioso

XBI=read.csv(file1, row.names=1, check.names = FALSE)

XB=XBI[,4:115]
row_names_df_to_remove1<-c("2599", "4038", "983", "3222", "3955", "3651", "440", "5160")
YB=read.csv(file2, row.names=1)
YB=as.data.frame(YB[,11]) 
joinB=cbind(XB,YB)
joinB=(joinB[!(row.names(joinB) %in% row_names_df_to_remove1),])
XB=joinB[,1:112]
YB=as.data.frame(joinB[,113])

#Import Cappricia

XC=read.csv(file3, row.names=1, check.names = FALSE)
row_names_df_to_remove<-c("i2T11S5", "i2T5S3", "i1T1S3", "i1T1S2")
row_names_df_to_remover<-c("i1T1S3", "i1T1S2")
YC=read.csv(file4, row.names=1)
YC=as.data.frame(YC[,11]) 
joinC=cbind(XC,YC)
joinC=as.data.frame(joinC[!(row.names(joinC) %in% row_names_df_to_remove),])
joinC=as.data.frame(joinC[!(row.names(joinC) %in% row_names_df_to_remover),])
XC=joinC[,1:112]
YC=as.data.frame(joinC[,113])

#Import Provine

XP=read.csv(file5, row.names=1, check.names = FALSE)
XP=XP[,4:115]
YP=read.csv(file6, row.names=1)
YP=as.data.frame(YP[,11]) 
joinP=cbind(XP,YP)
row_names_df_to_remove2<-c("7790", "9639", "0", "7232","10340", "1002", "4117")
joinP=joinP[!(row.names(joinP) %in% row_names_df_to_remove2),]
XP=joinP[,1:112]
YP=as.data.frame(joinP[,113])

#Write the same column names for every matrix


colnames(XC)=colnames(XB)
colnames(XP)= colnames(XB)
colnames(XC)=colnames(XP)
colnames(YC)=colnames(YP)

colnames(YB)="Class"
colnames(YC)="Class"
colnames(YP)="Class"

# Build one data set by combining two different varieties together


X=rbind(XC, XP)

Y=rbind(YC,YP)
row.names(Y)=row.names(X)
snv=as.data.frame(snv(X))

#First or Second Derivative, to SNV or to raw data
SGC=as.data.frame(savgol(snvC,2,17,2))


# Apply pretreatments

# Standar Normal Variate (SNV)

snv=as.data.frame(snv(X))

# First or Second Derivative

SGC=as.data.frame(savgol(snv,2,17,2))

# Detrend 


wav <- as.numeric(colnames(X))
DT=detrend(X, wav, p = 2)

# Plot Spectra

# You can plot Raw Spectra, or pretreated spectra


## Plots

par(mar = c(6, 6, 6, 6))

col  = "red"

plotsp(snv, ylab = "Measured Reflectance",xlab = "Wavelength (nm)",type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE)

title(main = " SNV NIR spectra",line = NA, outer = FALSE)


#Join  X and Y to build only one table 

join=cbind(Y,SGC)
row_names_df_to_remove<-c("2599")
join=join[!(row.names(join) %in% row_names_df_to_remove),]
nv=112


#Split dataset 70/30 in a representative way for each class

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



## According to the optimal model found previously, we select the most important variables by CovSel

  i =6
  nvar = i
  vs=(covsel(xtrain, ytrain, nvar, scaly = TRUE, weights = NULL))
  variables=vs$sel[1]
  variables[1]
  
  
## We build the matrix with the important variables found by CovSel
  v <- unlist(variables)
  v = sort(v)
  CovSelTrain =xtrain[,c(v)]
  CovSelTest =xtest[,c(v)]
  
  
  
#PLSDA modeling  
## Split TRAINING SET into VALIDATION AND TUNING, in order to choose the number of latent variables
#Please, add the optimal number of samples belonging to the Training set in  indices
  
  OneTable= cbind(CovSelTrain,ytrain)
  
  noseparo = TRUE 
  #while (noseparo)
  { 
    
    indices= sample(nrow(CovSelTrain), nrow(CovSelTrain) *0.30)
    
    validation =  OneTable[indices,]
    tuning  =  OneTable[-indices,]   
    
    Yvalidation=validation[,(nvar+1)]                    
    Xvalidation=validation[,-(nvar+1)]                
    
    
    Ytuning=tuning[,(nvar+1)]                    
    Xtuning=tuning[,-(nvar+1)] 
    #
    Ytuning=as.data.frame(Ytuning)
    #row.names(Ytuning)=row.names(Xtuning)
    #Yvalidation=as.data.frame(Yvalidation)
    #row.names(Yvalidation)=row.names(Xvalidation)
    ## Train a PLSDA model 
    
    m <- 50 ; p <- 8
    
    nlv <- nvar-1
    res <- gridscorelv(Xvalidation, Yvalidation, Xtuning, Ytuning, score = err,fun = plsrda,nlv = 1:nlv, verb = TRUE)
    
    
    plotscore(res$nlv, res$y1,main = "ERR", xlab = "Nb. LVs", ylab = "Value")
    u<- res[res$y1 == min(res$y1), ][1, , drop = FALSE]
    ytrain=as.data.frame(ytrain)
    YTrainInteger = as.integer(ytrain[,1])
    
    
    fmD <- plsrda(CovSelTrain,  YTrainInteger, nlv = u$nlv)
    
    pred <- predict(fmD, CovSelTest)$pred
    noseparo = (var(pred) == 0)
    
    
  }
  
  d=fmD$fm
  plot(d$P)
  plot(d$T)
  
  
  #The optimal number of latent variables was chosen according to the lower error in prediction
  #With 33 initial variables, the optimal value of latent variables is 19
  
  NumberofLatentVariables=u$nlv
  NumberofLatentVariables
  
  mx=as.matrix(ytest)
  
  
  ###### Now we can see the Results of the Predictions
  
  mx=as.matrix(ytest)
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


  ###### Export results to a project folder

  write.csv(e, "c:/temp/OPT/Workshop.csv")




  
  