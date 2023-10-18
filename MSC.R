
library(mdatools)
library(rchemo)
library (EMSC)

file="PROVINEa.csv"

#Import data
data=read.csv(file, row.names=1, check.names=FALSE)

data=data[,4:115]

#Change data to class matrix
data=as.matrix(data)

#Apply Multiplicative Scatter Correction
corrdata=prep.msc(data)

#Apply Standard Normal Variate
snv=snv(data)

#Apply Extended Multiplicative Scatter correction
emsc=EMSC(data, degree = 6)

#Plot the results

par(mar = c(6, 6, 6, 6))
par(mfrow = c(1, 3)) 
col = "red"
plotsp(data, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "Before MSC", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=3, cex.lab=2, cex.axis=1.5)
plotsp(corrdata, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After MSC", lwd = 2, type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=3, cex.lab=2, cex.axis=1.5)
#plotsp(snv, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After SNV", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=3, cex.lab=2, cex.axis=1.5)
plotsp(emsc$corrected, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After EMSC", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=3, cex.lab=2, cex.axis=1.5)



par(mar = c(6, 6, 6, 6))
par(mfrow = c(1, 2)) 
col = "red"
#plotsp(data, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "Before MSC", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=3, cex.lab=2, cex.axis=1.5)
plotsp(corrdata, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After MSC", lwd = 2, type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=2, cex.lab=1.5, cex.axis=1)
#plotsp(snv, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After SNV", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=3, cex.lab=2, cex.axis=1.5)
plotsp(emsc$corrected, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After EMSC, degree=6", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=2, cex.lab=1.5, cex.axis=1)

emsc6=EMSC(data, degree = 6)
#emsc4=EMSC(data, degree = 4)
emsc2=EMSC(data, degree = 2)

par(mar = c(6, 6, 6, 6))
par(mfrow = c(1, 2)) 
col = "red"
plotsp(emsc6$corrected, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After EMSC, degree=6", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=2, cex.lab=1.5, cex.axis=1)
#plotsp(emsc4$corrected, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After EMSC, degree=4", lwd = 2, type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=2, cex.lab=1.5, cex.axis=1)
plotsp(emsc2$corrected, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After EMSC, degree=2", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=2, cex.lab=1.5, cex.axis=1.5)
#plotsp(emsc$corrected, ylab = "Reflectance",xlab = "Wavelength (nm)",main = "After EMSC, degree=6", type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE, cex.main=2, cex.lab=1.5, cex.axis=1)




