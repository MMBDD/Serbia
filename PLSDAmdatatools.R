library(mdatools)
library(pcv)


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
YB=as.data.frame(YB[,12])
joinB=cbind(XB,YB)
joinB=(joinB[!(row.names(joinB) %in% row_names_df_to_remove1),])
XB=joinB[,1:112]
YB=as.data.frame(joinB[,113])

#Import Cappricia

XC=read.csv(file3, row.names=1, check.names = FALSE)
YC=read.csv(file4, row.names=1)
YC=as.data.frame(YC[,12])
joinC=cbind(XC,YC)


row_names_df_to_remove<-c("i2T11S5", "i2T5S3", "i1T1S3", "i1T1S2")
row_names_df_to_remover<-c("i1T1S3", "i1T1S2")
joinC=as.data.frame(joinC[!(row.names(joinC) %in% row_names_df_to_remove),])
joinC=as.data.frame(joinC[!(row.names(joinC) %in% row_names_df_to_remover),])

XC=joinC[,1:112]
YC=as.data.frame(joinC[,113])

#Import Provine

XP=read.csv(file5, row.names=1, check.names = FALSE)
XP=XP[,4:115]
YP=read.csv(file6, row.names=1)
YP=as.data.frame(YP[,12])

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

#Join  X and Y to build only one table

join=cbind(Y,X)

nv=112


#Split dataset 70/30 in a representative way for each class

JoinSepals <- split(join, join$Class)

xtrain = data.frame()

ytrain = as.integer()
xtest = data.frame()
ytest = as.integer()


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

###############################################################################################################

ytr = factor(ytrain, labels = c("Deaseased", "Healthy"))
yte = factor(ytest, labels = c("Deaseased", "Healthy"))

# make a plot with spectra
par(mfrow = c(2, 1))
mdaplot(xtrain, type = "l", cgroup = ytr, main = "Train")
mdaplot(xtest, type = "l", cgroup = yte, main = "Test")

# PLSDA
m = plsda(xtrain, ytr, ncomp = 40, cv = 1)
summary(m)
plot(m)

# I can not find calres
getConfusionMatrix(m$calres)  #I can not find similar information to calres in my model

par(mfrow = c(1, 2))
plotPredictions(m)


par(mfrow = c(3, 2))
plotMisclassified(m, nc = 2)

plotSensitivity(m, nc = 2)

plotSpecificity(m, nc = 2)


par(mfrow = c(1, 2))
plotRegcoeffs(m, ncomp = 3)


#Prediction on new data
res = predict(m, xtest, yte)
summary(res)

par(mfrow = c(1, 1))
plotPredictions(res)

par(mfrow = c(1, 2))
plotXResiduals(res)
plotYVariance(res)



############################################

#SIMCA


S= simca(xtrain, "Healthy", ncomp = 3)
summary(S)
plot(S)


layout(matrix(c(1, 3, 2, 3), ncol = 2))
plotSensitivity(S, show.labels = TRUE)
plotMisclassified(S, show.labels = TRUE)
plotPredictions(S, show.labels = TRUE)

Xpv = pcvpca(as.matrix(xtrain), 4, center = TRUE, scale = TRUE, cv = list("ven", 4))


S = simca(xtrain, "Diseased", ncomp = 3, x.test = Xpv)

par(mfrow = c(1, 2))
plotSensitivity(S, show.line = c(NA, 0.95))
plotVariance(S, type = "h", show.labels = TRUE)

S = selectCompNum(S, 1)
summary(S)




res = predict(S, xtest, yte)
summary(res)


par(mfrow = c(2, 2))
plotSpecificity(res, show.labels = TRUE)
plotSensitivity(res, show.labels = TRUE)
plotMisclassified(res, show.labels = TRUE)
plotPredictions(res)

show(res$c.pred[45:55, 1:3, 1])
show(getConfusionMatrix(res))
show(round(res$p.pred[45:55, 1:3, 1], 4))

par(mfrow = c(2, 1))
plotProbabilities(res, cgroup = yte)
plotProbabilities(res, ncomp = 2, cgroup = yte)

######################################
#SVM




