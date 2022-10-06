# REDUNDANCY PLOTS
# https://rstudio-pubs-static.s3.amazonaws.com/259028_dbe846c67e144065b4c2bcd012fb130d.html

#Clear workspace
rm(list = ls())

#Set working directory
setwd("C:/Users/Lauren/Documents/PhD/Publications/WesternCarnivoreDietReview")
getwd()

#Load required packages
library(ggplot2)
library(vegan)

#Read in the full %FO diet dataset
FO <- read.csv("WeightedAverage.csv")

#Cut just needed columns
FO.cut <-FO[,c(2,13:17)]
str(FO.cut)
head(FO.cut)

#Visualize the variables
pairs(FO.cut, lower.panel = NULL, col = as.numeric(iris$Species))

#Run principal components analysis (PCA) and use -1 to drop the species column
FO.pca <- princomp(FO.cut[,-1])
biplot(FO.pca)

#Now do a vegan plot
FO.vegan <- rda(FO.cut[,-1])
biplot(FO.vegan)     

#Add some options to the biplot
biplot(FO.vegan,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))

#Create a vector of species names
spp.names <- levels(FO.cut$Species)

#Add hulls around each species group
ordiellipse(FO.vegan, group = FO.cut$Species, col = c(1,2,3,4,5), kind = "sd", conf = 0.95)

legend("topright",
       col = c(1,2,3,4,5), 
       lty = 1,
       lwd = 3,
       legend = spp.names)

