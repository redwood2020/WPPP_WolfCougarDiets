# Clear workspace
rm(list=ls())

# Set Working Directory
setwd("C:/Users/laure/OneDrive/Desktop/Ch1_DietAnalysis/ClusterDietData")
getwd()

#### Functions
AICc.PERMANOVA <- function(adonis.model) {
  
  # check to see if object is an adonis model...
  
  if (!(adonis.model$aov.tab[1,1] >= 1))
    stop("object not output of adonis {vegan} ")
  
  # Ok, now extract appropriate terms from the adonis model
  # Calculating AICc using residual sum of squares (RSS) since I don't think that adonis returns something I can use as a liklihood function...
  
  RSS <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "SumsOfSqs"]
  MSE <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "MeanSqs"]
  
  k <- ncol(adonis.model$model.matrix)# + 1 # add one for error variance
  
  nn <- nrow(adonis.model$model.matrix)
  
  # AIC : 2*k + n*ln(RSS)
  # AICc: AIC + [2k(k+1)]/(n-k-1)
  
  # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
  # https://www.researchgate.net/post/What_is_the_AIC_formula;
  # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
  
  # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
  # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
  
  AIC <- 2*k + nn*log(RSS)
  AIC.g <- 2*k + nn * (1 + log( 2 * pi * RSS / nn))
  AIC.MSE <- 2*k + nn * log(MSE)
  AIC.pi <- k + nn*(1 + log( 2*pi*RSS/(nn-k) )   )
  AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
  AICc.MSE <- AIC.MSE + (2*k*(k + 1))/(nn - k - 1)
  AICc.pi <- AIC.pi + (2*k*(k + 1))/(nn - k - 1)
  
  output <- list("AIC" = AIC, "AIC.g" = AIC.g, "AICc" = AICc,
                 "AIC.MSE" = AIC.MSE, "AICc.MSE" = AICc.MSE,
                 "AIC.pi" = AIC.pi, "AICc.pi" = AICc.pi, "k" = k)
  
  return(output)   
  
}
AICc.PERMANOVA2 <- function(adonis2.model) {
  
  # check to see if object is an adonis2 model...
  
  if (is.na(adonis2.model$SumOfSqs[1]))
    stop("object not output of adonis2 {vegan} ")
  
  # Ok, now extract appropriate terms from the adonis model Calculating AICc
  # using residual sum of squares (RSS or SSE) since I don't think that adonis
  # returns something I can use as a likelihood function... maximum likelihood
  # and MSE estimates are the same when distribution is gaussian See e.g.
  # https://www.jessicayung.com/mse-as-maximum-likelihood/;
  # https://towardsdatascience.com/probability-concepts-explained-maximum-likelihood-estimation-c7b4342fdbb1
  # So using RSS or MSE estimates is fine as long as the residuals are
  # Gaussian https://robjhyndman.com/hyndsight/aic/ If models have different
  # conditional likelihoods then AIC is not valid. However, comparing models
  # with different error distributions is ok (above link).
  
  
  RSS <- null$SumOfSqs[ length(null$SumOfSqs) - 1 ]
  MSE <- RSS / adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
  
  nn <- adonis2.model$Df[ length(adonis2.model$Df) ] + 1
  
  k <- nn - adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
  
  
  # AIC : 2*k + n*ln(RSS/n)
  # AICc: AIC + [2k(k+1)]/(n-k-1)
  
  # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
  # https://www.statisticshowto.datasciencecentral.com/akaikes-information-criterion/ ;
  # https://www.researchgate.net/post/What_is_the_AIC_formula;
  # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf;
  # https://medium.com/better-programming/data-science-modeling-how-to-use-linear-regression-with-python-fdf6ca5481be 
  
  # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
  # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
  
  AIC <- 2*k + nn*log(RSS/nn)
  AIC.g <- 2*k + nn * (1 + log( 2 * pi * RSS / nn))
  AIC.MSE <- 2*k + nn * log(MSE)
  AIC.pi <- k + nn*(1 + log( 2*pi*RSS/(nn-k) )   )
  AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
  AICc.MSE <- AIC.MSE + (2*k*(k + 1))/(nn - k - 1)
  AICc.pi <- AIC.pi + (2*k*(k + 1))/(nn - k - 1)
  
  output <- list("AIC" = AIC, "AICc" = AICc, "AIC.g" = AIC.g, 
                 "AIC.MSE" = AIC.MSE, "AICc.MSE" = AICc.MSE,
                 "AIC.pi" = AIC.pi, "AICc.pi" = AICc.pi, "k" = k, "N" = nn)
  
  return(output)   
  
}

#####

library(tidyverse)

data<-read.csv("PERMANOVA_Query_Categories.csv")
head(data)
str(data)

forage<-data[,6:18]
head(forage)
# Forage data is currently the total duration of video in each forage group
# To convert to a proportion of each deer's total foraging time in each forage group, do
# First convert NA values to zero (0)
forage[is.na(forage)] <- 0
forage<-forage/rowSums(forage)

ind.log<-log(forage+1)
ind.nmds<-vegan::metaMDS(ind.log,"bray",k=2,autotransform = F,trymax = 1000, na.rm=T)

species<-data$Species
sa<-data$StudyArea
sex<-data$Sex
species.sa<-paste(data$Species, data$StudyArea)
species.sex<-paste(data$Species, data$Sex)

# Cougar - Okanogan vs Northeast
plot(ind.nmds, main =  "Cougar")
vegan::ordihull(ind.nmds,groups=(species.sa== "Cougar Okanogan"),draw="polygon",
                col=c("skyblue2"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Cougar Northeast"),draw="polygon",
                col=c("tomato2"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)
with(data, legend("bottomleft", legend = c("Okanogan","Northeast"), bty = "n",
                    pch = c(3,2)))

# Wolf - Okanogan vs Northeast
plot(ind.nmds, main =  "Wolf")
vegan::ordihull(ind.nmds,groups=(species.sa== "Wolf Okanogan"),draw="polygon",
                col=c("skyblue2"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Wolf Northeast"),draw="polygon",
                col=c("tomato2"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)

plot(vegan::anosim(vegan::vegdist(forage),species.sa, distance = "bray"), main = "ANOSIM - Cluster Data")
vegan::anosim(vegan::vegdist(forage),species.sa, distance = "bray")

null<-vegan::adonis2(forage~1,data=forage)
speciesmod<-vegan::adonis2(forage~species,data=forage)
samod<-vegan::adonis2(forage~sa,data=forage)
linmod<-vegan::adonis2(forage~species+sa,data=forage)
full<-vegan::adonis2(forage~species+sa+species:sa,data=forage)

tabs<-rbind(
unlist(AICc.PERMANOVA2(null)) ,
unlist(AICc.PERMANOVA2(speciesmod)),
unlist(AICc.PERMANOVA2(samod)),
unlist(AICc.PERMANOVA2(linmod)),
unlist(AICc.PERMANOVA2(full)))

tabs

# Northeast - Cougar vs Wolf
plot(ind.nmds, main =  "Northeast (Cougar vs. Wolf)")
vegan::ordihull(ind.nmds,groups=(species.sa== "Cougar Northeast"),draw="polygon",
                col=c("green3"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Wolf Northeast"),draw="polygon",
                col=c("cornflowerblue"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)

# Okanogan - Cougar vs Wolf
plot(ind.nmds, main =  "Okanogan (Cougar vs. Wolf)")
vegan::ordihull(ind.nmds,groups=(species.sa== "Cougar Okanogan"),draw="polygon",
                col=c("green3"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Wolf Okanogan"),draw="polygon",
                col=c("cornflowerblue"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)

# Cougar - Male vs Female
plot(ind.nmds, main =  "Cougar (Male vs Female)")
vegan::ordihull(ind.nmds,groups=(species.sex== "Cougar Male"),draw="polygon",
                col=c("blue"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sex== "Cougar Female"),draw="polygon",
                col=c("darkgoldenrod1"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)

# Wolf - Male vs Female
plot(ind.nmds, main =  "Wolf (Male vs Female)")
vegan::ordihull(ind.nmds,groups=(species.sex== "Wolf Male"),draw="polygon",
                col=c("blue"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sex== "Wolf Female"),draw="polygon",
                col=c("darkgoldenrod1"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)


# General Plot Call
plot.new()
plot.window(xlim = c(-2.0,2.0), ylim = c(-2.0,2.0))
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")


#Cougar - Okanogan vs Northeast
with(forage, points(ind.nmds$points[species.sa== "Cougar Okanogan",],
                       pch = 3))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "Cougar Okanogan",1]),
                          mean(ind.nmds$points[species.sa== "Cougar Okanogan",2])),
                    pch = 13,lwd=1,cex=2))
with(forage, points(ind.nmds$points[species.sa== "Cougar Northeast",],
                    pch = 2))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "Cougar Northeast",1]),
                          mean(ind.nmds$points[species.sa== "Cougar Northeast",2])),
                    pch = 14,cex=2))
op <- par(cex = 1)
with(forage, legend("topright", col=c("skyblue2", "tomato2"), legend = c("Okanogan","Northeast"), bty = "n",
                      pch = c(3,2)))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main="Cougar (Okanogan vs. Northeast)")
box()

#Wolf - Okanogan vs Northeast
with(forage, points(ind.nmds$points[species.sa== "Wolf Okanogan",],
                    pch = 3))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "Wolf Okanogan",1]),
                          mean(ind.nmds$points[species.sa== "Wolf Okanogan",2])),
                    pch = 13,cex=2))
with(forage, points(ind.nmds$points[species.sa== "Wolf Northeast",],
                    pch = 2))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "Wolf Northeast",1]),
                          mean(ind.nmds$points[species.sa== "Wolf Northeast",2])),
                    pch = 14,cex=2))
op <- par(cex = 1)
with(forage, legend("topright", col=c("skyblue2", "tomato2"), legend = c("Okanognan","Northeast"), bty = "n",
                    pch = c(3,2)))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main="Wolf (Okanogan vs. Northeast)")
box()

#Northeast - Cougar vs Wolf
with(forage, points(ind.nmds$points[species.sa== "Cougar Northeast",],
                    pch = 3))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "Cougar Northeast",1]),
                          mean(ind.nmds$points[species.sa== "Cougar Northeast",2])),
                    pch = 13,lwd=1,cex=2))
with(forage, points(ind.nmds$points[species.sa== "Wolf Northeast",],
                    pch = 2))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "Wolf Northeast",1]),
                          mean(ind.nmds$points[species.sa== "Wolf Northeast",2])),
                    pch = 14,cex=2))
op <- par(cex = 1)
with(forage, legend("topright", col=c("green3", "cornflowerblue"), legend = c("Cougar","Wolf"), bty = "n",
                    pch = c(3,2)))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main="Northeast (Cougar vs. Wolf)")
box()

#Okanogan - Cougar vs Wolf
with(forage, points(ind.nmds$points[species.sa== "Cougar Okanogan",],
                    pch = 3))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "Cougar Okanogan",1]),
                          mean(ind.nmds$points[species.sa== "Cougar Okanogan",2])),
                    pch = 13,lwd=1,cex=2))
with(forage, points(ind.nmds$points[species.sa== "Wolf Okanogan",],
                    pch = 2))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "Wolf Okanogan",1]),
                          mean(ind.nmds$points[species.sa== "Wolf Okanogan",2])),
                    pch = 14,cex=2))
op <- par(cex = 1)
with(forage, legend("topright", col=c("green3", "cornflowerblue"), legend = c("Cougar","Wolf"), bty = "n",
                    pch = c(3,2)))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main="Okanogan (Cougar vs. Wolf)")
box()

#Cougar - Male vs Female
with(forage, points(ind.nmds$points[species.sex== "Cougar Male",],
                    pch = 3))
with(forage, points(cbind(mean(ind.nmds$points[species.sex== "Cougar Male",1]),
                          mean(ind.nmds$points[species.sex== "Cougar Male",2])),
                    pch = 13,lwd=1,cex=2))
with(forage, points(ind.nmds$points[species.sex== "Cougar Female",],
                    pch = 2))
with(forage, points(cbind(mean(ind.nmds$points[species.sex== "Cougar Female",1]),
                          mean(ind.nmds$points[species.sex== "Cougar Female",2])),
                    pch = 14,cex=2))
op <- par(cex = 1)
with(forage, legend("topright", col=c("blue", "darkgoldenrod1"),  legend = c("Male","Female"), bty = "n",
                    pch = c(3,2)))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main="Cougar (Male vs. Female)")
box()

#Wolf - Male vs Female
with(forage, points(ind.nmds$points[species.sex== "Wolf Male",],
                    pch = 3))
with(forage, points(cbind(mean(ind.nmds$points[species.sex== "Wolf Male",1]),
                          mean(ind.nmds$points[species.sex== "Wolf Male",2])),
                    pch = 13,lwd=1,cex=2))
with(forage, points(ind.nmds$points[species.sex== "Wolf Female",],
                    pch = 2))
with(forage, points(cbind(mean(ind.nmds$points[species.sex== "Wolf Female",1]),
                          mean(ind.nmds$points[species.sex== "Wolf Female",2])),
                    pch = 14,cex=2))
op <- par(cex = 1)
with(forage, legend("topright", col=c("blue", "darkgoldenrod1"), legend = c("Male","Female"), bty = "n",
                    pch = c(3,2)))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main="Wolf (Male vs. Female)")
box()

