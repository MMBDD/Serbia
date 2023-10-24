#install.packages("remotes")
#install.packages("rchemo")
#remotes::install_github("mlesnoff/rchemo", dependencies = TRUE, build_vignettes = TRUE)

library(psych)
library(mt)
library(rchemo)
library(kernlab)
library(rrcovHD)
library(isotree)

## This code is made to create different tables with different outlier detection techniques
##Set working directory

#setwd("C:/Users/berto011/OneDrive - Wageningen University & Research/Proyectos/ANTARES-TOMATES/Cappricia/ANTARES-TOMATES")

# The aim of this R code is to remove outliers at pixel level and to average the remaining inliers
# This code was created to build one single file where each row is one sepal



##Import data. Please choose from varieties "Cappricia, Brioso, and Provine"
file = "Cappricia.csv"
#file= "Brioso.csv"
#file="Provine.csv"


Tomatos=read.csv(file, row.names=1, check.names = FALSE)

# As you can see, there are several number of rows for each sepal
# Each row is one pixel
# We will average all rows that belong to the same sepal, in order to have only one spectrum per sepal.

# You will have the possibility to choose different outlier detection methods

# Let's start by splitting the data set into different images

Images<- split(Tomatos, Tomatos$Plate_num)


#Now, it is time to select the outlier detection method:
# "PCA", "PCA+IF", "IF", "PCDIST", "SIMCA", "univariate", "SVM".

# Moreover, you can select the way to average the inliers or remaining pixels:
#1 calculating the average ("average")
#2 adding  the Standard Deviation as another feature, and then calculating the average ("averageSD")
#3 computing the Standard Normal Variate of spectra, then calculating the average ("SNV")


threshold=0.60
ODMethod="PCDIST"
AvgMethod ="SNV" 

#results data frame
dfResults = data.frame()

#for each image
for (image in Images)
{
  
  LocalImages =  split(image, image$Fruit_num)
  #for each fruit
  for ( Fruit in LocalImages)
  {
    Sepals = split(Fruit,Fruit$Sepal_num)
    #for each sepal
    for ( Sepal in Sepals)
    {
      #definitions
      SepalDef = Sepal[1,1:3]
      #data matrix
      #if Data set is provine:   SepalData =Sepal[,6:117]
      SepalData =Sepal[,4:115]
      # we make sure that every collumn is numeric.
      #SepalData <- SepalData %>% mutate_if(is.character, as.numeric)
      if (ODMethod=="PCA")
      {  
        pca=pca.outlier(SepalData, center = TRUE, scale=TRUE,conf.level = 0.90) 
        outliers=pca$outlier
      }else if (ODMethod=="IF")
      {
        modelIF <- isolation.forest(SepalData, ndim=5, ntrees=100, nthreads=1)
        scores <- predict(modelIF, SepalData)
        threshold=0.50
        outliersTmp <- (ifelse(scores > threshold, 0, 1))
        outliers = which(outliersTmp %in% 0)
      }else if (ODMethod =="PCDIST")
      {
        suppressWarnings({ 
          obj <- OutlierPCDist(SepalData)
        })
        outliers = which(getFlag(obj) %in% 0)
        getCutoff(obj)
      }
      
      else if (ODMethod=="SIMCA")
      {
        suppressWarnings({ 
          obj <- OutlierPCDist(SepalData)
        })
        outliers = which(getFlag(obj) %in% 0)
        if (length(outliers)>0 )
        {
          SepalDataNoOut = SepalData[-outliers,]
        }else
        {
          SepalDataNoOut = SepalData
        }
        SIMCA=RSimca(SepalDataNoOut)

        pr <- predict(SIMCA, newdata=SepalData)
        SIMCA@pcaobj
        pr@classification
        
        outliers = which( pr@classification %in% 1)
        
        
      }else if (ODMethod=="univariate")
      {
        rowmean = rowMeans(SepalData)
        quartiles <- quantile(rowmean, probs=c(.25, .75), na.rm = FALSE)
        IQR <- IQR(rowmean)
        dfrowmean = data.frame()
        Lower <- quartiles[1] - 1.5*IQR
        Upper <- quartiles[2] + 1.5*IQR 
        
        dfrowmean = rbind(dfrowmean, rowmean)
        dfrowmean = t(dfrowmean)
        dfrowmean$row_num <- seq.int(nrow(dfrowmean)) 
        
        dfrowmean$outlier <- rowmean > Lower & rowmean < Upper
        out = data.frame()
        out = rbind(out,dfrowmean$row_num, dfrowmean$outlier)
        out = t(out)
        colnames(out) <- c("v1","v2")
        out = as.data.frame(out)
        outliers =  out[out$v2 == 0,1 ]
        
      }else if (ODMethod=="SVM")
      {
        pca=pca.outlier(SepalData, center = TRUE, scale=TRUE,conf.level = 0.95) 
        outliersPreProcess=pca$outlier
        if ( length( outliersPreProcess)>0 ) 
        {
          SepalWOOutlier=SepalData[-outliers,]
        }else
        {
          SepalWOOutlier=SepalData 
        }
        model_svm <- ksvm(as.matrix(SepalWOOutlier), type="one-svc",nu = 0.5)
        pred_svm <- predict(model_svm, as.matrix(SepalData), type="decision")
        #results_svm <- data.frame(Model = "One-Class SVM",AUROC = auc(-pred_svm, y))
        threshold <- quantile(pred_svm, 0.90)
        outliersTmp <- (ifelse(pred_svm > threshold, 0, 1))
        outliers = which(outliersTmp %in% 0)
        # outliers=
        
      }else if(ODMethod=="PCA+IF")
      {
        pca=prcomp(SepalData, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL, rank.=10)
        pcascores <- pca$x
        modelIF <- isolation.forest( pcascores, ndim=5, ntrees=100, nthreads=1)
        scores <- predict(modelIF,  pcascores)
        threshold=0.50
        outliersTmp <- (ifelse(scores > threshold, 0, 1))
        outliers = which(outliersTmp %in% 0)
      }
      else
      {
        #error
      }
      
      
      
      if (length(outliers)>0 )
      {
        #substact outliers
        SepalWOOutlier=SepalData[-outliers,]
      }else
      {
        SepalWOOutlier=SepalData
      }
      if (AvgMethod=="average")
      {
        
        colmean = colMeans(SepalWOOutlier)
        #apply colMean
        colmean = as.data.frame(colmean)
        colmean =  t(colmean)
        #we build a line result
        result = cbind(SepalDef,colmean)
        
        
        
      }else if (AvgMethod=="SNV")
      {
        SNV = as.data.frame(snv(SepalWOOutlier))
        
        colmean = colMeans(SNV)
        #apply colMean
        colmean = as.data.frame(colmean)
        colmean =  t(colmean)
        #we build a line result
        result = cbind(SepalDef,colmean)
        
      }else if (AvgMethod=="averageSD")
      {
        StD=SD(SepalWOOutlier)
        SepalWOOutlier2=rbind(SepalWOOutlier,StD)
        colmean2 = colMeans(SepalWOOutlier2)
        #apply colMean
        colmean2 = as.data.frame(colmean2)
        colmean2 =  t(colmean2)
        #we build a line result
        result = cbind(SepalDef,colmean2)
      }
      #we add the result line to the result data frame
      dfResults = rbind(dfResults,result)
    }
  }
}

# Export the Result table
# Please, write the file name in accordance to the methods you have chosen
# Save the results in a project folder

write.csv(dfResults, "c:/temp/pruebaNN.csv")



