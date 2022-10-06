##PIANKA'S INDEX

#Clear workspace
rm(list = ls())

#Install and load packages
library(EcoSimR)

#Set working directory
setwd("C:/Users/Lauren/Documents/PhD/Publications/WesternCarnivoreDietReview")
getwd()

#Example with simulated data (https://www.rdocumentation.org/packages/EcoSimR/versions/0.1.0/topics/pianka)
m <- matrix(rpois(80, 1), nrow = 10) #rows are species, columns are prey groups
m
pianka(m)

#Read in carnivore diet data
diets <- read.csv("Piankas_Index.csv")

#View the data
diets

#Create a new matrix called "pianka" from the diet data
pianka <- diets

#Strip off the column headers
colnames(pianka) <- NULL

#View the resulting data
pianka

#Create new matrix without the row names
pianka <- pianka[,2:6]
pianka <- as.matrix(pianka, nrow=5)

#View the new data frame. There should not be any row or column headers.
pianka

#Run the pianka() function on the all five species in the data
pianka(pianka)

#Subset out certain species pairs or groups to test in the order they appear in the "diets" matrix
#Do this by adjusting the nubmers in the list c()
#Rows: 1=bobcat, 2=cougar, 3=coyote, 4=grey wolf, 5=mexican wolf
#Pianka's Index numbers are from 0 to 1. The closer to 1, the more the overlap

#Example comparing cougars to coyotes
coug_v_coy <- pianka[c(2,3),]
coug_v_coy
pianka(coug_v_coy) #0.543

#Example comparing cougars to grey wolves
coug_v_wolf <- pianka[c(2,4),]
coug_v_wolf
pianka(coug_v_wolf) #0.977

#Example comaring cougar, bobcats, and coyotes
coug_bob_coy <- pianka[c(1,2,3),]
coug_bob_coy
pianka(coug_bob_coy) #0.652

#Example comaring cougar, bobcats, coyotes, and wolves (excluding mex. wolves)
all <- pianka[c(1,2,3,4),]
all
pianka(all) #0.630

#########################################

#Coyote vs bobcat
coy_v_bob<- pianka[c(1,3),]
coy_v_bob
pianka(coy_v_bob) #0.543

#Coyote vs wolf
coy_v_wolf<- pianka[c(3,4),]
coy_v_wolf
pianka(coy_v_wolf) #0.485

#Coyote vs cougar
coug_v_coy <- pianka[c(2,3),]
coug_v_coy
pianka(coug_v_coy) #0.543

#Cougar vs bobcat
coug_v_bob<- pianka[c(2,1),]
coug_v_bob
pianka(coug_v_bob) #0.438

#Cougar vs wolf
coug_v_wolf <- pianka[c(2,4),]
coug_v_wolf
pianka(coug_v_wolf) #0.977

#Bobcat vs wolf
bob_v_wolf<- pianka[c(1,4),]
bob_v_wolf
pianka(bob_v_wolf) #0.365

#Wolf vs mexican wolf
wolf_v_m.wolf<- pianka[c(4,5),]
wolf_v_m.wolf
pianka(wolf_v_m.wolf) #0.995
