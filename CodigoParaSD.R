

library(prospectr)
library(caret)
library(rchemo)
library(ggplot2)



file="DatosparaSD.csv"
  
X=read.csv(file, row.names=1, check.names = FALSE, sep=";")

fifteen <- subset(X,X$Temp==15)
twentyfive <- subset(X,X$Temp==25)
thirty <- subset(X,X$Temp==30)
#time <- ifelse(fifteen$Date %in% "11/17/2023"), 
                 # 1,
                  #ifelse(fifteen$Date %in% "11/20/2023"),
                   #      2,
                    #     ifelse(fifteen$Date %in% "11/22/2023"), 
                     #           3,ifelse(fifteen$Date %in% "11/27/2023"),
                      #          4, ifelse(fifteen$Date %in% "11/29/2023"),
                       #                   5, ifelse(fifteen$Date %in% "12-1-2023"),
                        #                            6, ifelse(fifteen$Date %in% "12-6-2023"),
                         #                                     7,
                          #                                    ifelse(fifteen$Date %in% "12-8-2023"),
                           #                                          8,
                            #                                         ifelse(fifteen$Date %in% "12-13-2023"),
                             #                                               9,
                              #                                              ifelse(fifteen$Date %in% "12-15-2023"),
                               #                                                    10,
                                #                                                   ifelse(fifteen$Date %in% "12-19-2023"),
                                 #                                                  11)



# Plot

## Plots

par(mar = c(6, 6, 6, 6))

col  = "red"

plotsp(fifteen, ylab = "15",xlab = "Masses",type = "l", col = col, zeroes = FALSE, labels = FALSE, add = FALSE)

title(main = "CG-MS spectra",line = NA, outer = FALSE)

ggplot(fifteen) +
  geom_line(aes(x = fifteen$Day, y= fifteen$`45`),
             size=rel(3.0))

             