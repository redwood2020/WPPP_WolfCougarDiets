# Clear workspace
rm(list=ls())

getwd()

`#### Functions
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

data<-read.csv("Data/diets_long_PERMANOVA.csv")
head(data)
str(data)

# note that I manually added values to the "prey_simple_deerspp" and "prey_simple_unkdeer" columns in Excel to create the diets_long_PERMANOVA.csv from the diets_long.csv - code this in for future (9/7/22)

# add in season, cluster/non-cluster, and carcass found/not found
scatlog<-read.csv("Data/ScatLog_MetabarcodingSamples.csv")
head(scatlog)
str(scatlog)

# create a new column based on the "Found_Along" column named "Cluster" - Yes = found "AtCluster", No = otherwise
table(scatlog$Found_Along)
scatlog <- scatlog %>%
  mutate(Cluster = recode(Found_Along, 'AtCluster' = 'Yes', 'AtCluster, Off-trail' = 'Yes', 'AtCluster, Road'='Yes', 'AtCluster, Trail' = 'Yes', 'Follow' = 'No', 'Off-trail' = 'No', 'Other' = 'No', 'Road' = 'No', 'Trail' = 'No'))
with(scatlog, table(Metabarcoding_Species, Cluster))
table(scatlog$Cluster)
head(scatlog)
str(scatlog)

# reclassify "Season" into just 'Summer' (Summer + Fall) and 'Winter' (Winter + Spring)
table(scatlog$Season)
scatlog <- scatlog %>%
  mutate(Season = recode(Season, 'Winter' = 'Winter', 'Spring' = 'Winter', 'Summer'='Summer', 'Fall' = 'Summer'))
with(scatlog, table(Season, Metabarcoding_Species))
table(scatlog$Season)
table(scatlog$Date_Sent_OSU)
head(scatlog)

# create column for carcass found/not found at cluster in the scatlog database so we know how many scats came from clusters with prey (true feeding sites) vs from clusters without prey (presumed resting sites)
# subset the scatlog to only scats that were sent for metabarcoding analysis
scatlog_meta <- scatlog[!(scatlog$Date_Sent_OSU == ""), ]
# create a new column for cluster ID in the scat log
scatlog_meta <- scatlog_meta %>% mutate(Cluster_ID = gsub('(.*)-\\w+', '\\1', Scat_ID), .after = Scat_ID)
# read in the cluster database files
clusters_reg <- read.csv("Data/WolfCougarClusters_Database.csv")
clusters_winter <- read.csv("Data/Winter2020_Database.csv")
head(clusters_reg)
head(clusters_winter)
dim(clusters_reg)
dim(clusters_winter)
# merge the two databases into one
clusters <- dplyr::bind_rows(clusters_reg, clusters_winter)
#assign carcass found/not found to each scat ID by cluster ID
scatlog_meta <- merge(scatlog_meta, clusters[, c("Cluster_ID", "CarcassFound","PreySpecies_Final")], by="Cluster_ID", all.x = TRUE)
with(scatlog_meta, table(Metabarcoding_Species, Metabarcoding_Status))

# creat a second database of only wolf and cougar scats with confirmed depositor that also contain vertebrate prey
slm_cc <- subset(scatlog_meta, Metabarcoding_Species == 'Wolf' | Metabarcoding_Species == 'Cougar')
slm_cc <- subset(slm_cc, Metabarcoding_Status == "Pred_w_prey")
# samples sizes for all scats with depostior confirmed regardless of prey contents
with(scatlog_meta, table(Season, StudyArea, Metabarcoding_Species))
# sample sizes for amplified scats with prey (subset of previous)
with(slm_cc, table(Season, StudyArea, Metabarcoding_Species))

# now calculate percent frequency of occurrence by scat sample
FO_coug <- data %>% filter(Depositor_DNA == "Puma concolor", StudyArea == "Okanogan") %>% count(SampleID, prey_simple_deerspp) %>% group_by(SampleID) %>% mutate(n = prop.table(n)) %>% ungroup() %>%
  pivot_wider(names_from = prey_simple_deerspp, values_from = n, names_prefix = '') %>% replace(is.na(.), 0)
# take the resulting table and turn it into 1's if the item was present and 0's otherwise
FO_coug_bi <- FO_coug %>% mutate_if(is.numeric, ~1 * (. > 0))
FO_coug_bi <- subset (FO_coug_bi, select = -SampleID)
FO_coug_spp <- sapply(FO_coug_bi,sum)
# now divide the column sums by the total number of scats (rows) to get the percent frequency of occurrence
FO_coug_spp <- FO_coug_spp/as.integer(nrow(FO_coug_bi))
FO_coug_spp <- stack(FO_coug_spp)
FO_coug_spp
hist(FO_coug_spp)
  
# ok so next step is to make this into a larger table with columns for "Depositor_DNA", "StudyArea", "PreySpp", "%FO"
# maybe make the above into a function - DepostiorSpecies, StudyArea, PreyCategoryType (ie prey_simple_deerspp or prey_simple_unkdeer)

#######################################

# rename the "Metabarcoding_ID" column to "SampleID" to match other database
scatlog <- rename(scatlog, SampleID = Metabarcoding_ID)
# add in the "Season" column from 'scatlog' to 'data'
data <- merge(data, scatlog[, c("SampleID", "Season")], by="SampleID")
# add in the "Cluster" column from 'scatlog' to 'data'
data <- merge(data, scatlog[, c("SampleID", "Cluster")], by="SampleID")
# add in the "Carcass" column from 'scatlog' to 'data'
# data <- merge(data, scatlog[, c("SampleID", "Carcass")], by="SampleID")
head(data)

# rearrange columns so "Season" and "Cluster" come after "StudyArea"
col_order <- c("SampleID", "StudyArea", "Season", "Cluster", "Depositor_Field", "Depositor_DNA", "diet_item", "prey_item", "prey_simple_deerspp", "prey_simple_unkdeer")
data.scat <- data[, col_order]
head(data.scat)
table(data.scat$diet_item)

# look at some summary tables
with(data.scat, table(Depositor_DNA, Season, StudyArea))

# convert to wide format
library(data.table)
# cover to wide format
data.wide.unkdeer <- data.scat %>% 
  group_by(SampleID) %>%
  pivot_wider(names_from = prey_simple_unkdeer, values_from = prey_simple_unkdeer)
data.wide.unkdeer <- as.data.table(data.wide.unkdeer) # need to convert to class data.table for this to work
# samples with multiple prey items are still showing up as multiple rows - collapse to one row per SampleID
View(data.wide.unkdeer)

data.wide.deerspp <- data.scat %>%
  group_by(SampleID) %>%
  pivot_wider(names_from = prey_simple_deerspp, values_from = prey_simple_deerspp)
data.wide.deerspp <- as.data.table(data.wide.deerspp)

# View(data.scat)
# View(data.wide.unkdeer) # with all deer as "deerunkspp" so we can compare wolf/cougar to bear (which does not have deer species differentiated)
# View(data.wide.deerspp) # with deer species divided into "muledeer", "whitetaileddeer", and "deerunkspp" for the wolf vs cougar comparison

str(data.wide.unkdeer)

#collapse
#unkdeer (where we only have deerunkspp/Odocoileus spp.)
data.wide.unkdeer.factors <- data.wide.unkdeer[,1:9]
data.wide.unkdeer.binary <- +!is.na(data.wide.unkdeer[,10:20])
data.wide.unkdeer <- cbind(data.wide.unkdeer.factors, data.wide.unkdeer.binary)
setDT(data.wide.unkdeer)
data.wide.unkdeer <- data.wide.unkdeer[ , .(deerunkspp = sum(deerunkspp), other = sum(other), fish = sum(fish), small_mammal = sum(small_mammal), bird = sum(bird), livestock = sum(livestock), elk = sum(elk), lagomorph = sum(lagomorph), carnivore = sum(carnivore), moose = sum(moose), med_mammal = sum(med_mammal)),
    by = .(SampleID, StudyArea, Season, Cluster, Depositor_Field, Depositor_DNA)]# sum rows by SampleID - sum numeric value only
# View(data.wide.unkdeer)

# for the "wide" prey item columns, convert NA values to 0 and non-NA values to 1, so from column 10 and up
#deerspp (where we have muledeer, whitetaileddeer, and deerunkspp / Odocoileus spp.)
data.wide.deerspp.factors <- data.wide.deerspp[,1:9]
data.wide.deerspp.binary <- +!is.na(data.wide.deerspp[,10:22])
data.wide.deerspp <- cbind(data.wide.deerspp.factors, data.wide.deerspp.binary)
str(data.wide.deerspp)
setDT(data.wide.deerspp)
data.wide.deerspp <- data.wide.deerspp[ , .(deerunkspp = sum(deerunkspp), other = sum(other), fish = sum(fish), small_mammal = sum(small_mammal), whitetaileddeer = sum(whitetaileddeer), bird = sum(bird), livestock = sum(livestock), elk = sum(elk), lagomorph = sum(lagomorph), muledeer = sum(muledeer), carnivore = sum(carnivore), moose = sum(moose), med_mammal = sum(med_mammal)),  by = .(SampleID, StudyArea, Season, Cluster, Depositor_Field, Depositor_DNA)] 

# now convert prey 0/1 values to fractions of the total (so that a scat with both deer and coyote will be 0.5 deer, 0.5 coyote)
#unkdeer
unkdeer.factors <- data.wide.unkdeer[,1:6]
unkdeer.rows <- data.wide.unkdeer[,7:17] # we dropped the columsn for diet_item, prey_item, and prey_simple_unkdeer so we reduce the column call by a value of 3
unkdeer.rows <- unkdeer.rows/rowSums(unkdeer.rows)
data.wide.unkdeer <- cbind(unkdeer.factors, unkdeer.rows)
# View(data.wide.unkdeer)

#deerspp
deerspp.factors <- data.wide.deerspp[,1:6]
deerspp.rows <- data.wide.deerspp[,7:19] # we dropped the columsn for diet_item, prey_item, and prey_simple_unkdeer so we reduce the column call by a value of 3
deerspp.rows <- deerspp.rows/rowSums(deerspp.rows)
data.wide.deerspp <- cbind(deerspp.factors, deerspp.rows)
# View(data.wide.deerspp)

############################################################
# SCAT DATA NOW READY FOR NMDS PLOTS
############################################################

# first, pick which dataset to work with (unkdeer or deerspp)
data <- data.wide.unkdeer # can change to data.wide.deerspp 

# set up for nmds
ind.log <- data[,7:17] # 7:17 for unkdeer, 7:19 for deerspp
##### ind.log<-log(forage+1) <- from Shannon's script but I think log(data+1) doesn't make sense here
par(ask = TRUE) # set plot
ind.nmds<-vegan::metaMDS(ind.log,"bray",k=2,autotransform = F,trymax = 1500, plot = TRUE, previous.best=TRUE, na.rm=T) # set trymax to 1000 for real run - shortened to increase processing time


# assign variables of interest
species<-data$Depositor_DNA
sa<-data$StudyArea
cluster<-data$Cluster
season<-data$Season
species.sa<-paste(data$Depositor_DNA, data$StudyArea)

# look at data structure and assign class
str(data)
# data$StudyArea <- as.factor(data$StudyArea)
# data$Season <- as.factor(data$Season)
# data$Cluster <- as.factor(data$Cluster)
# data$Depositor_Field <- as.factor(data$Depositor_Field)
# data$Depositor_DNA <- as.factor(data$Depositor_DNA)
# str(data)

# subset as desired  (e.g. one study area, one season, etc.)
# data <- subset(data, Season == "Summer") # if only looking at summer data
# View(data.new)
# length(data.new$SampleID)
str(data)

# Cougar - Okanogan vs Northeast
plot(ind.nmds, main =  "Cougar")
vegan::ordihull(ind.nmds,groups=(species.sa== "Puma concolor Northeast"),draw="polygon",
                col=c("tomato2"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Puma concolor Okanogan"),draw="polygon",
                col=c("skyblue2"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.001)
with(data, legend("bottomleft", legend = c("Okanogan","Northeast"), bty = "n",
                    pch = c(3,2)))

# Wolf - Okanogan vs Northeast
plot(ind.nmds, main =  "Wolf")
vegan::ordihull(ind.nmds,groups=(species.sa== "Canis lupus Northeast"),draw="polygon",
                col=c("tomato2"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Canis lupus Okanogan"),draw="polygon",
                col=c("skyblue2"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)

# Bear - Okanogan vs Northeast
plot(ind.nmds, main =  "Black Bear")
vegan::ordihull(ind.nmds,groups=(species.sa== "Ursus americanus Northeast"),draw="polygon",
                col=c("tomato2"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Ursus americanus Okanogan"),draw="polygon",
                col=c("skyblue2"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)

plot(vegan::anosim(vegan::vegdist(forage),species.sa, distance = "bray"), main = "ANOSIM - Scat Data")
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

# Northeast - Cougar vs Wolf vs Bear
plot(ind.nmds, main =  "Northeast (Cougar vs. Wolf vs. Black Bear)")
vegan::ordihull(ind.nmds,groups=(species.sa== "Wolf Northeast"),draw="polygon",
                col=c("cornflowerblue"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "BlackBear Northeast"),draw="polygon",
                col=c("gold"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Cougar Northeast"),draw="polygon",
                col=c("green3"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)

# Okanogan - Cougar vs Wolf vs Bear
plot(ind.nmds, main =  "Okanogan (Cougar vs. Wolf vs. Bear)")
vegan::ordihull(ind.nmds,groups=(species.sa== "Cougar Okanogan"),draw="polygon",
                col=c("green3"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Wolf Okanogan"),draw="polygon",
                col=c("cornflowerblue"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "BlackBear Okanogan"),draw="polygon",
                col=c("gold"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)

# Northeast - Cougar vs Wolf
plot(ind.nmds, main =  "Northeast (Cougar vs. Wolf)")
vegan::ordihull(ind.nmds,groups=(species.sa== "Wolf Northeast"),draw="polygon",
                col=c("cornflowerblue"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Cougar Northeast"),draw="polygon",
                col=c("green3"), show.groups=T,label=F)
vegan::orditorp(ind.nmds,display="species",col="black",air=0.01)

# Okanogan - Cougar vs Wolf
plot(ind.nmds, main =  "Okanogan (Cougar vs. Wolf)")
vegan::ordihull(ind.nmds,groups=(species.sa== "Cougar Okanogan"),draw="polygon",
                col=c("green3"), show.groups=T,label=F)
vegan::ordihull(ind.nmds,groups=(species.sa== "Wolf Okanogan"),draw="polygon",
                col=c("cornflowerblue"), show.groups=T,label=F)
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

#Bear - Okanogan vs Northeast
with(forage, points(ind.nmds$points[species.sa== "BlackBear Okanogan",],
                    pch = 3))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "BlackBear Okanogan",1]),
                          mean(ind.nmds$points[species.sa== "BlackBear Okanogan",2])),
                    pch = 13,cex=2))
with(forage, points(ind.nmds$points[species.sa== "BlackBear Northeast",],
                    pch = 2))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "BlackBear Northeast",1]),
                          mean(ind.nmds$points[species.sa== "BlackBear Northeast",2])),
                    pch = 14,cex=2))
op <- par(cex = 1)
with(forage, legend("topright", col=c("skyblue2", "tomato2"), legend = c("Okanognan","Northeast"), bty = "n",
                    pch = c(3,2)))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main="Black Bear (Okanogan vs. Northeast)")
box()

#Northeast - Cougar vs Wolf vs Bear
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
with(forage, points(ind.nmds$points[species.sa== "BlackBear Northeast",],
                    pch = 0))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "BlackBear Northeast",1]),
                          mean(ind.nmds$points[species.sa== "BlackBear Northeast",2])),
                    pch = 12,cex=2))
op <- par(cex = 1)
with(forage, legend("topright", col=c("green3", "cornflowerblue", "gold"), legend = c("Cougar","Wolf","Black Bear"), bty = "n",
                    pch = c(3,2,0)))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main="Northeast (Cougar vs. Wolf vs. Black Bear)")
box()

#Okanogan - Cougar vs Wolf vs Black Bear
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
with(forage, points(ind.nmds$points[species.sa== "BlackBear Okanogan",],
                    pch = 0))
with(forage, points(cbind(mean(ind.nmds$points[species.sa== "BlackBear Okanogan",1]),
                          mean(ind.nmds$points[species.sa== "BlackBear Okanogan",2])),
                    pch = 12,cex=2))
op <- par(cex = 1)
with(forage, legend("topright", col=c("green3", "cornflowerblue", "gold"), legend = c("Cougar","Wolf", "Black Bear"), bty = "n",
                    pch = c(3,2,0)))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main="Okanogan (Cougar vs. Wolf vs. Black Bear)")
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
