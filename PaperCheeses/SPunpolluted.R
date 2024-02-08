
#Install and open libraries

install.packages("dplyr")


library(dplyr)
library(rchemo)
library(EMSC)
library(mdatools)


#This is the complete file

file = "totalFX17df.csv"

#This is the summarized file

file2="resultadosampleNumber.csv"
colnames()

ObjectsperSample <- read.csv(file2, row.names = NULL, check.names = FALSE)
Ch <- read.csv(file, row.names = NULL, check.names = FALSE)

colnames(Ch)



# Divide the complete file per sample and object, then compute the average
ObjectsperSample <- Ch %>%
  group_by(Samplenumber, Object) %>%
  summarize(across(everything(), mean))

# Divide the complete file per File name and object, then compute the average

ObjectsperFile <- Ch %>%
  group_by(Filenumber, Object) %>%
  summarize(across(everything(), mean))



#Save the summarized files

write.csv(resultados, "c:/temp/resultadosFileNumber.csv")

write.csv(ObjectsperSample, "c:/temp/resultadosampleNumber.csv")


#Build X and Y matrices

X=as.data.frame(ObjectsperSample[,3:114])
Yp=as.data.frame(ObjectsperSample[,122])

Yf=as.data.frame(ObjectsperSample[,119])

Yz=as.data.frame(ObjectsperSample[,123])

#Plot spectra

par(mar = c(6, 6, 8, 6))
par(mfrow = c(1, 1)) 
col = "red"
plotsp(log(1/X), ylab = "Absorbance",xlab = "Wavelength (nm)",main = "Raw NIR Spectra", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=2, cex.lab=2, cex.axis=1.5)



par(mar = c(9, 6, 9, 6))
par(mfrow = c(1, 1)) 
col = "red"
#plotsp(emsc$corrected, ylab = "Absorbance",xlab = "Wavelength (nm)",main = "After EMSC d=6", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=2, cex.lab=2, cex.axis=1.5)
plotsp(savgol, ylab = "Absorbance",xlab = "Wavelength (nm)",main = "After SNV and SD (2,17,2)", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=2, cex.lab=2, cex.axis=1.5)




#Apply Extended Multiplicative Scatter correction
savgol=savgol(X,2,17,2)

#Apply Standard Normal Variate
snv=snv(log(1/X))

#Apply Extended Multiplicative Scatter correction
emsc=EMSC(log(1/X), degree = 6)

#Apply Extended Multiplicative Scatter correction
savgol=savgol(snv,2,17,2)

#Plot Pretreated spectra

par(mar = c(5, 6, 4, 10) + 0.1)
par(mfrow = c(3, 1)) 
col = "red"
plotsp(snv, ylab = "Absorbance",xlab = "Wavelength (nm)",main = "After SNV", lwd = 2, type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=3, cex.lab=2, cex.axis=1.5)
plotsp(emsc$corrected, ylab = "Absorbance",xlab = "Wavelength (nm)",main = "After EMSC, degree=6", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=3, cex.lab=2, cex.axis=1.5)
plotsp(savgol, ylab = "Absorbance",xlab = "Wavelength (nm)",main = "After Second Derivative", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=3, cex.lab=2, cex.axis=1.5)


##Explore by PCA
pca=prcomp(savgol, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL, rank.=15)
summary(pca)

##Calculate scores and loadings
pcascores <- pca$x

pcaloadings <- pca$rotation

##Plot PCA scores

PCAcolors <- c("#66c2a5","#fc8d62","#8da0cb", "#E5C494", "#B3B3B3")
plot(pcascores[,1:2], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,1:2],labels=rownames(pcascores))


plot(pcascores[,2:3], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,2:3],labels=rownames(pcascores))

plot(pcascores[,c(1,3)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(1,3)],labels=rownames(pcascores))

plot(pcascores[,c(4,3)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(4,3)],labels=rownames(pcascores))

plot(pcascores[,c(4,5)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(4,5)],labels=rownames(pcascores))

plot(pcascores[,c(2,4)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(2,4)],labels=rownames(pcascores))

plot(pcascores[,c(6,5)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(6,5)],labels=rownames(pcascores))

plot(pcascores[,c(6,7)], pch=21, col=PCAcolors, bg=PCAcolors, cex=1.5)
text(pcascores[,c(6,7)],labels=rownames(pcascores))


##REMOVE OUTLIERS  #No outliers to remove this time

#row_names_df_to_remove<-c("51", "193")
#X=X[!(row.names(X) %in% row_names_df_to_remove),]
#Yp=Yp[!(row.names(Yp) %in% row_names_df_to_remove),]

#Change data to class matrix

#Feature selection

IPLS=ipls(savgol,Yf,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)



#Glucide

# Utilizar na.omit() para eliminar filas con NA en la columna A

# Utilizar na.omit() para eliminar filas con NA en la columna A

Objectsper <- ObjectsperSample [complete.cases(ObjectsperSample$GLUCIDE), ]


Xg=as.data.frame(Objectsper[,3:114])
Yg=as.data.frame(Objectsper[,124])

snvg=snv(Xg)
savgolg=savgol(snvg,2,17,2)

IPLS=ipls(savgolg,Yg,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)


IPLS=ipls(snvg,Yg,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(Xg,Yg,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

savgolgg=savgol(Xg,2,17,2)
IPLS=ipls(savgolgg,Yg,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

savgolggg=savgol(Xg,1,17,2)
IPLS=ipls(savgolggg,Yg,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

emscg=EMSC(Xg, degree = 6)
IPLS=ipls(emscg$corrected,Yg,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

# ZOUT

ObjectsperSample <- Ch %>%
  group_by(Samplenumber, Object) %>%
  summarize(across(everything(), mean))



Xz=ObjectsperSample[,3:114]
Yz=ObjectsperSample[,123]


snv=snv(Xz)
savgol2=savgol(snv,2,17,2)
savgol1=savgol(snv,1,17,2)
emsc=EMSC(Xz, degree = 6)
savgol22=savgol(Xz,2,17,2)


IPLS=ipls(snv,Yz,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(savgol2,Yz,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(savgol1,Yz,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(emsc$corrected,Yz,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(savgol22,Yz,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)


IPLS=ipls(Xv,Yv,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)



#VETTEN



Xv=ObjectsperSample[,3:114]
Yv=ObjectsperSample[,119]

snv=snv(Xv)
savgol2=savgol(snv,2,17,2)
savgol1=savgol(snv,1,17,2)
emsc=EMSC(Xz, degree = 6)
savgol22=savgol(Xv,2,17,2)
savgol11=savgol(Xv,1,17,2)

IPLS=ipls(snv,Yv,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(savgol2,Yv,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(savgol1,Yv,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(emsc$corrected,Yv,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(savgol22,Yv,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)

IPLS=ipls(savgol11,Yv,glob.ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 10),exclcols = NULL,exclrows = NULL, int.ncomp = 11, int.num = 5,int.width = NULL,int.limits = NULL,int.niter = NULL, ncomp.selcrit = "min",  method = "forward",x.test = NULL,y.test = NULL,silent = FALSE)
IPLS$int.limits
summary(IPLS)




#Fat content

#Data Split



set.seed(1)

dt = sort(sample(nrow(ObjectsperSample), nrow(ObjectsperSample)*.7))
train<-ObjectsperSample[dt,]
test<-ObjectsperSample[-dt,]

dim(train)
dim(test)

# create training set and test set

xtrain=train[,3:114]
xtest=test[,3:114]
ytrain=as.data.frame(train[,119])
ytest=as.data.frame(test[,119])

dim(xtrain)
dim(xtest)
dim(ytest)
dim(ytrain)

xtrain=snv(xtrain)
xtest=snv(xtest)
xtrain=savgol(xtrain,2,17,2)
xtest=savgol(xtest,2,17,2)

#xtrain=EMSC(xtrain, degree = 6)
#xtest=EMSC(xtest, degree = 6)
#vs=(covsel(xtrain, ytrain, 15, scaly = TRUE, weights = NULL))
#variables=vs$sel[1]
#variables[1]


## columns subset
#v <- unlist(variables)
#v = sort(v)
#xtrain =xtrain[,c(v)]
#xtest=xtest[,c(v)]


#Pretreatment

#xtrain=snv(xtrain)
#xtest=snv(xtest)

#xtrain=savgol(xtrain, 1, 17, 2)
#xtest=savgol(xtest,1,17,2)

#xtest=snv(xtest)
#xtrain=snv(xtrain)


#xtest=savgol((snv(xtest)),2,17,2)
#xtrain=savgol((snv(xtrain)),2,17,2)



## Divide Training set into cal and val 
n=93

OneTable= cbind(xtrain,ytrain)


indices2= sample(n, 40)
Tuning =  OneTable[indices2,]
Cal  =  OneTable[-indices2,]                  

Ycal=Cal[,113]                    
Xcal=Cal[,1:112]                        



Ycaldf= as.data.frame(Ycal) 

Ytuning=Tuning[,113]                    
Xtuning=Tuning[,1:112]                        

dim(Xtuning)                      

Ytuningdf= as.data.frame(Ytuning) 


#Cross validation


nlvdis <- c(20) ; diss <- "mahal"

nlv <- 0:20
res <- gridscorelv(Xcal, Ycal, Xtuning, Ytuning, fun =plskern, score = rmsep, nlv = nlv, verb = TRUE)
u <- res[res$y1 == min(res$y1), ][1, , drop = FALSE]
plotscore(res$nlv, res$y1,main = "ERR", xlab = "Nb. LVs", ylab = "Value")
u$nlv

#Fit and compare regression models

pls=plsnipals(xtrain$corrected, ytrain, weights = NULL, nlv = u$nlv)

predPLS <- predict(pls, xtest$corrected)$pred
rmsepPLS=rmsep(predPLS, ytest)
msepPLS=msep(predPLS, ytest)
sepPLS=sep(predPLS, ytest)
biasPLS=bias(predPLS, ytest)
cor2PLS=cor2(predPLS, ytest)
r2PLS=r2(predPLS, ytest)
r2PLS

msepPLS
sepPLS
biasPLS
rmsepPLS

cor2PLS
results=as.data.frame(cbind(r2PLS,msepPLS, sepPLS, biasPLS, rmsepPLS,u$nlv ))
row.names(results)="f"
colnames(results)=c("R2Pred","MSEP", "SEP", "BiAS", "RMSEP", "Lvs")

results




#Verzadigde

#We use a Table without NA values

Objectsper <- ObjectsperSample [complete.cases(ObjectsperSample$GLUCIDE), ]

set.seed(1)

dt = sort(sample(nrow(Objectsper), nrow(Objectsper)*.7))
train<-Objectsper[dt,]
test<-Objectsper[-dt,]

dim(train)
dim(test)

# create training set and test set

xtrain=train[,3:114]
xtest=test[,3:114]
ytrain=as.data.frame(train[,119])
ytest=as.data.frame(test[,119])


#xtrain=snv(xtrain)
#xtest=snv(xtest)
#xtrain=savgol(xtrain, 2,17,2)
#xtest=savgol(xtest, 2,17,2)

dim(xtrain)
dim(xtest)
dim(ytest)
dim(ytrain)

#xtrain=detrend(xtrain, degree = 2)
#xtest=detrend(xtest,degree=2)
#xtest=snv(xtest)
#xtrain=savgol(xtrain,1,17,2)
#xtest=savgol(xtest,1,17,2)

#xtrain=EMSC(xtrain, degree = 6)
#xtest=EMSC(xtest, degree = 6)
#vs=(covsel(xtrain, ytrain, 15, scaly = TRUE, weights = NULL))
#variables=vs$sel[1]
#variables[1]


## columns subset
#v <- unlist(variables)
#v = sort(v)
#xtrain =xtrain[,c(v)]
#xtest=xtest[,c(v)]


#Pretreatment

#xtrain=snv(xtrain)
#xtest=snv(xtest)

#xtrain=savgol(xtrain, 1, 17, 2)
#xtest=savgol(xtest,1,17,2)

#xtest=snv(xtest)
#xtrain=snv(xtrain)


#xtest=savgol((snv(xtest)),2,17,2)
#xtrain=savgol((snv(xtrain)),2,17,2)



## Divide Training set into cal and val 
n=93

OneTable= cbind(xtrain,ytrain)


indices2= sample(n, 30)
Tuning =  OneTable[indices2,]
Cal  =  OneTable[-indices2,]                  

Ycal=Cal[,113]                    
Xcal=Cal[,1:112]                        



Ycaldf= as.data.frame(Ycal) 

Ytuning=Tuning[,113]                    
Xtuning=Tuning[,1:112]                        

dim(Xtuning)                      

Ytuningdf= as.data.frame(Ytuning) 


#Cross validation


nlvdis <- c(20) ; diss <- "mahal"

nlv <- 0:20
res <- gridscorelv(Xcal, Ycal, Xtuning, Ytuning, fun =plskern, score = rmsep, nlv = nlv, verb = TRUE)
u <- res[res$y1 == min(res$y1), ][1, , drop = FALSE]
plotscore(res$nlv, res$y1,main = "ERR", xlab = "Nb. LVs", ylab = "Value")
u$nlv

#Fit and compare regression models

pls=plsnipals(xtrain, ytrain, weights = NULL, nlv = u$nlv)

predPLS <- predict(pls, xtest)$pred
rmsepPLS=rmsep(predPLS, ytest)
msepPLS=msep(predPLS, ytest)
sepPLS=sep(predPLS, ytest)
biasPLS=bias(predPLS, ytest)
cor2PLS=cor2(predPLS, ytest)
r2PLS=r2(predPLS, ytest)
r2PLS

msepPLS
sepPLS
biasPLS
rmsepPLS

cor2PLS
results=as.data.frame(cbind(r2PLS,msepPLS, sepPLS, biasPLS, rmsepPLS,u$nlv ))
row.names(results)="f"
colnames(results)=c("R2Pred","MSEP", "SEP", "BiAS", "RMSEP", "Lvs")

results



# Feature Selection by CovSel




i =5
to = 40
results = data.frame()
results = data.frame(matrix(nrow = 5, ncol = 0)) 

# Set the total number of observations
n <- 93

while (i < to) {
  nvar <- i
  
  # Covsel algorithm 
  vs <- covsel(xtrain, ytrain, nvar, scaly = TRUE, weights = NULL)
  variables <- vs$sel[1]
  
  # Columns subset
  v <- unlist(variables)
  v <- sort(v)
  
  # Select the important variables from the original matrices
  CovSelTrain <- xtrain[, v]
  CovSelTest <- xtest[, v]
  
  OneTable <- cbind(CovSelTrain, ytrain)
  noSplit <- TRUE 
  
  while (noSplit) { 
    indices <- sample(n, 30)
    validation <- OneTable[indices, ]
    tuning <- OneTable[-indices, ]   
    
    Yvalidation <- validation[, (nvar + 1)]                    
    Xvalidation <- validation[, -(nvar + 1)] 
    Ytuning <- tuning[, (nvar + 1)]                    
    Xtuning <- tuning[, -(nvar + 1)] 
    
    # Train a PLSDA model 
    nlv <- nvar - 1
    res <- gridscorelv(Xtuning, Ytuning, Xvalidation, Yvalidation, fun = plskern, score = rmsep, nlv = nlv, verb = TRUE)
    u <- res[res$y1 == min(res$y1), ][1, , drop = FALSE]
    
    pls <- plsnipals(xtrain, ytrain, weights = NULL, nlv = u$nlv)
    predPLS <- predict(pls, xtest)$pred
    
    noSplit <- (var(predPLS) == 0)
  }
  
  rmsepPLS <- rmsep(predPLS, ytest)
  msepPLS <- msep(predPLS, ytest)
  sepPLS <- sep(predPLS, ytest)
  biasPLS <- bias(predPLS, ytest)
  cor2PLS <- cor2(predPLS, ytest)
  r2PLS <- r2(predPLS, ytest)
  
  # Convert the numeric vectors into character vectors
  r2PLS <- as.character(r2PLS)
  msepPLS <- as.character(msepPLS)
  sepPLS <- as.character(sepPLS)
  biasPLS <- as.character(biasPLS)
  rmsepPLS <- as.character(rmsepPLS)
  
  # Create a matrix from the vectors
  #comparison_matrix <- matrix(c(r2PLS, msepPLS, sepPLS, biasPLS, rmsepPLS), ncol = 5)
  
  # Transpose the matrix
  #comparison_matrix <- t(comparison_matrix)
  
  # Convert the transposed matrix into a data frame
  #comparativeTable <- as.data.frame(comparison_matrix)
  
  # Set row names for the data frame
  #rownames(comparativeTable) <- c("r2PLS", "msepPLS", "sepPLS", "biasPLS", "rmsepPLS")
  
  # Set column names for the data frame
  #colnames(comparativeTable) <- paste("Variable", 1:ncol(comparativeTable))
  
  
  # Append the results to the main data frame
  results <- rbind(results, comparativeTable)
  
  # Increment the counter
  #i <- i + 1
 
  rownames(comparativeTable) <- c("r2PLS", "msepPLS", "sepPLS", "biasPLS", "rmsepPLS")
  colnames(comparativeTable)<- (nvar)
  results = cbind(results,comparativeTable)
  results
  i = i+1
  
  
}


# Print the results
print(results)


