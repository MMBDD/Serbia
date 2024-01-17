library(mdatools)
library(pcv)
library(rchemo)

#Global Model with all varieties together


file="BriosoaSIMCA.csv"
file2="ProvineSIMCA.csv"
file3="CappriciaSIMCA.csv"

# Import

XBI=read.csv(file, row.names=1, check.names = FALSE)
XP=read.csv(file2, row.names=1, check.names = FALSE)

XC=read.csv(file3, row.names=1, check.names = FALSE)

#Prepare calibration and validation sets
#Brioso

xtest=XBI[1:50,5:116]
ytest=XBI[1:50,4]


xtrain=snv(xtrain)

xtune=XBI[50:76,5:116]
ytune=XBI[50:76,4]


xtune=snv(xtune)


xtrain=XBI[77:159,5:116]
ytrain=XBI[77:159,4]
xtrain=snv(xtest)


#Cappricia


xtrain=XC[1:49,3:114]
ytrain=XC[1:49,2]


xtrain=snv(savgol(xtrain,2,17,2))

xtune=XC[50:77,3:114]
ytune=XC[50:77,2]


xtune=snv(xtune)


xtest=XC[78:166,3:114]
ytest=XC[78:166,2]
xtest=snv(xtest)




#Provine

xtrain=XP[,5:116]
ytrain=XP[,4]


xtrain=snv(xtrain)

xtune=XP[50:76,5:116]
ytune=XP[50:76,4]


xtune=snv(xtune)


xtest=XP[77:159,5:116]
ytest=XP[77:159,4]
xtest=snv(xtest)
###############################################################################################################

ytr = factor(ytrain, labels = c("Healthy"))
yte = factor(ytest, labels = c("Diseased"))
ytu = factor(ytune, labels = c("Healthy"))
# make a plot with spectra
par(mfrow = c(1, 2))
mdaplot(xtrain, type = "l", cgroup = ytr, main = "Train")
mdaplot(xtest, type = "l", cgroup = yte, main = "Test")



############################################

#SIMCA


S= simca(xtrain, "Healthy", ncomp = 3)
summary(S)
plot(S)


layout(matrix(c(1, 3, 2, 3), ncol = 2))
plotSensitivity(S, show.labels = TRUE)
plotMisclassified(S, show.labels = TRUE)
plotPredictions(S, show.labels = TRUE)


#////////////////////////////////////////////////////////////////////////////////////////
  
library(pcv)
Xpv = pcvpca(as.matrix(xtrain), 4, center = TRUE, scale = TRUE, cv = list("ven", 4))


m = simca(xtrain, "Healthy", ncomp = 3, x.test = xtrain)  


par(mfrow = c(1, 2))
plotSensitivity(m, show.line = c(NA, 0.95))
plotVariance(m, type = "h", show.labels = TRUE)


m = selectCompNum(m, 2)
summary(m)

res = predict(m, xtest, yte)
summary(res)



par(mfrow = c(2, 2))
plotSpecificity(res, show.labels = TRUE)
plotSensitivity(res, show.labels = TRUE)
plotMisclassified(res, show.labels = TRUE)
plotPredictions(res)

show(getConfusionMatrix(res))
