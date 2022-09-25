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
with(scatlog_meta, table(CarcassFound, Cluster, Metabarcoding_Species))

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
  
# ok so next step is to make this into a larger table with columns for "Depositor_DNA", "StudyArea", "PreySpp", "%FO"
# maybe make the above into a function - DepostiorSpecies, StudyArea, PreyCategoryType (ie prey_simple_deerspp or prey_simple_unkdeer)
# first try to turn the %FO code into a function
# note that "preygrouptype" is either "prey_simple_deerspp" or "prey_simple_unkdeer" - the former differentiates between Odocoileus spp and the later groups all deer into "Odocoileus" / deerunkspp

# test criteria
depositorspp = "Puma concolor"
studyarea = "Okanogan"
season = "Summer"
  

# rename "Metabarcoding_ID" column in the scatlog_meta file to "Sample_ID to match the corresponding column in the "data" file
slm_cc <- rename(slm_cc, SampleID = Metabarcoding_ID)
# add a column for season to the "data" file by pulling from the "slm_cc" file
data <- merge(data, slm_cc[, c("SampleID", "Season", "Cluster")], by.x = "SampleID")
# rearrage column order so "Season" comes after "StudyArea"
col_order <- c("SampleID", "StudyArea", "Season", "Cluster", "Depositor_Field", "Depositor_DNA", "diet_item", "prey_item", "prey_simple_deerspp", "prey_simple_unkdeer")
data <- data[, col_order]

# attempt to reassign "deerunkspp" to either "muledeer" or "whitetaileddeer" in the "prey_simple_deerspp" column of the "data" dataframe - assign proportionally within species, study area, and season
## first calculate proportions of MD and WTD species within each group
MD <- data %>% group_by(Depositor_DNA, StudyArea, Season) %>% summarise( priv_perc = sum(prey_simple_deerspp == "muledeer", na.rm=T) / sum(prey_simple_deerspp == "whitetaileddeer" | prey_simple_deerspp == "muledeer", na.rm=T) )

WTD <- data %>% group_by(Depositor_DNA, StudyArea, Season) %>% summarise( priv_perc = sum(prey_simple_deerspp == "whitetaileddeer", na.rm=T) / sum(prey_simple_deerspp == "whitetaileddeer" | prey_simple_deerspp == "muledeer", na.rm=T) )

# create a new column only to indicate rows where "prey_simple_deerspp" = "deerunkspp" and automatically set to "deerunkspp" and NA otherwise
data <- data %>% mutate(Status = case_when(
  prey_simple_deerspp == "deerunkspp" ~ "deerunkspp",    prey_simple_deerspp != "deerunkspp" ~ "unknown"
))
head(data)

# now select your carnivore/studyarea/season combo
carnivore = "Canis lupus"
studyarea = "Northeast"
season = "Winter"

pct_d <- WTD[WTD$Depositor_DNA == carnivore & WTD$StudyArea == studyarea & WTD$Season == season, 4]
as.numeric(pct_d)
1-pct_d

# get number of values of "deerunkspp" in Status, assign n
nd <- table(data$Status)[1] 
# create a list of random binary variables where pct_d represents the proportion of white-taileddeer, and 1-pct_d represent the number of mule deer
WTD
data$random_prop_deer <- data[ , 'random_prop_deer'] <- NA
nd_names <- c()
head(data)
nd_names <- sample(c("whitetaileddeer", "muledeer"),
       size = nd, 
       prob = c(pct_d, (1-pct_d)), replace = TRUE)

# now create a new column called "rand_prop_deer" that fills in all the "deerunkspp" in the "preys_simple_deerspp" column with either "whitetaileddeer" or "muledeer" based on proportion of other known deer for that carnivore-studyarea-season group 
deerloop <- function(carnivore, studyarea, season) { 
  WTD <- data %>% group_by(Depositor_DNA, StudyArea, Season) %>%   
    summarise( priv_perc = sum(prey_simple_deerspp == "whitetaileddeer", na.rm=T) / sum(prey_simple_deerspp == "whitetaileddeer" | prey_simple_deerspp == "muledeer", na.rm=T) )
  pct_d <- WTD[WTD$Depositor_DNA == carnivore & WTD$StudyArea == studyarea & WTD$Season == season, 4]
  data_sub <- filter(data, Depositor_DNA == carnivore, StudyArea == studyarea, Season == season)

  for (i in 1:length(data_sub$Status)) {
  nd_names[i] <- sample(c("whitetaileddeer", "muledeer"),
                     size = length(data_sub$Status[i]), 
                     prob = c(pct_d, (1-pct_d)), replace = TRUE)
  data_sub$rand_prop_deer[i]<-ifelse(data_sub$Status[i] == "deerunkspp",nd_names[i], data_sub$prey_simple_deerspp[i])
  }
  return(data_sub)
}

out_CN_W <- deerloop(carnivore = "Canis lupus", studyarea = "Northeast", season = "Winter")
# check that it filtered correctly
with(out, table(StudyArea, Season, Depositor_DNA))

# ok now use the function to generate outputs for each  of the eight groups (carnivore, studyarea, season) and then rowbind to merge
# wolves
out_CN_S <- deerloop(carnivore = "Canis lupus", studyarea = "Northeast", season = "Summer")
out_CN_W <- deerloop(carnivore = "Canis lupus", studyarea = "Northeast", season = "Winter")
out_CO_S <- deerloop(carnivore = "Canis lupus", studyarea = "Okanogan", season = "Summer")
out_CO_W <- deerloop(carnivore = "Canis lupus", studyarea = "Okanogan", season = "Winter")
# cougars
out_PN_S <- deerloop(carnivore = "Puma concolor", studyarea = "Northeast", season = "Summer")
out_PN_W <- deerloop(carnivore = "Puma concolor", studyarea = "Northeast", season = "Winter")
out_PO_S <- deerloop(carnivore = "Puma concolor", studyarea = "Okanogan", season = "Summer")
out_PO_W <- deerloop(carnivore = "Puma concolor", studyarea = "Okanogan", season = "Winter")

# now combine into a new dataframe
data_new <- rbind(out_CN_S, out_CN_W, out_CO_S, out_CO_W, out_PN_S, out_PN_W, out_PO_S, out_PO_W)

data <- data_new #reassign to "data" because I'm too lazy to change the varible in the following code
head(data)


# create percent frequency of occurrence (Pct_FO) function
Pct_FO <- function(data, depositorspp, studyarea, season) {
  library(dplyr)
  library(tidyr)
  FO <- data %>% dplyr::filter(Depositor_DNA == depositorspp, StudyArea == studyarea, Season == season) %>% count(SampleID, rand_prop_deer) %>% group_by(SampleID) %>% mutate(n = prop.table(n)) %>% ungroup() %>%
    pivot_wider(names_from = rand_prop_deer, values_from = n, names_prefix = '') %>% replace(is.na(.), 0)
  # take the resulting table and turn it into 1's if the item was present and 0's otherwise
  FO_bi <- FO %>% mutate_if(is.numeric, ~1 * (. > 0))
  FO_bi <- subset (FO_bi, select = -SampleID)
  FO_spp <- sapply(FO_bi,sum)
  # now divide the column sums by the total number of scats (rows) to get the percent frequency of occurrence
  FO_spp <- FO_spp/as.integer(nrow(FO_bi))
  FO_spp <- stack(FO_spp)
  FO_spp <- as.data.frame(FO_spp)
  colnames(FO_spp) <- c("Pct_FO","PreyItem")
  FO_spp <- FO_spp[c("PreyItem", "Pct_FO")]
  FO_spp <- FO_spp %>% mutate(StudyArea = studyarea, .before = PreyItem) %>% mutate(DepositorSpp = depositorspp, .before = StudyArea) %>% mutate(Season = season, .after = StudyArea)
  assign( (paste0(sub(" .*", "", depositorspp), "_", studyarea, "_", season)) , FO_spp)
  
}

# create dataframes for each species/study area combo

# cougars / Puma concolor
PN_S <- Pct_FO(data, "Puma concolor", "Northeast", "Summer")
PN_W <- Pct_FO(data, "Puma concolor", "Northeast", "Winter")
PO_S <- Pct_FO(data, "Puma concolor", "Okanogan", "Summer")
PO_W <- Pct_FO(data, "Puma concolor", "Okanogan", "Winter")

# wolves / Canis lupus
CN_S <- Pct_FO(data, "Canis lupus", "Northeast", "Summer")
CN_W <- Pct_FO(data, "Canis lupus", "Northeast", "Winter")
CO_S <- Pct_FO(data, "Canis lupus", "Okanogan", "Summer")
CO_W <- Pct_FO(data, "Canis lupus", "Okanogan", "Winter")

PN <- Pct_FO(data, "Puma concolor", "Northeast")
PO <- Pct_FO(data, "Puma concolor", "Okanogan")
CN <- Pct_FO(data, "Canis lupus", "Northeast")
CO <- Pct_FO(data, "Canis lupus", "Okanogan")
# combine dataframes
FO_spp_all <- rbind(PN_S, PN_W)
# plot
# https://community.rstudio.com/t/help-with-making-plot-with-multiple-columns/50763/2
# by depositor species
ggplot(FO_spp_all, aes(DepositorSpp, Pct_FO, fill = PreyItem)) + geom_col(position = "dodge")

str(FO_spp_all$PreyItem)
# order legend prey items (no more 'deerunkspp')
FO_spp_all$PreyItem <- factor(FO_spp_all$PreyItem, levels=c('bird', 'elk', 'carnivore', 'lagomorph', 'med_mammal', 'moose', 'muledeer', 'small_mammal', 'whitetaileddeer'))

# combine dataframes
FO_spp_all <- rbind(CO_S, CO_W)
# by season
p <- ggplot(FO_spp_all, aes(Season, Pct_FO, fill = PreyItem, label = PreyItem)) + geom_col(position = position_dodge2(width = 0.9, preserve = "single")) + scale_fill_manual(values = c("bird" = "aquamarine3", "elk" = "brown", "carnivore" = "purple","lagomorph" = "darkseagreen","med_mammal" = "chocolate1", "moose" = "burlywood4", "muledeer" = "deepskyblue", "small_mammal" = "chartreuse2", "whitetaileddeer" = "darkgreen", "livestock" = "black")) + geom_text(position = position_dodge2(width = 0.9, preserve = "single"), angle = 90, vjust=0.35, hjust= -0.05)
#change plot title ("Cougar - Northeast") based on selected FO_spp_all combined dataframes
p4 <- p + scale_y_continuous(limits = c(0, 1.5)) + ggtitle("Wolf - Okanogan") + theme_bw() + theme(plot.title=element_text(hjust=0.5)) + ylab("% Frequency of Occurrence") +  guides(fill=guide_legend(title="Prey Item"))

# create a multiplot with gridarrange()
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

# how to create multiple plots in r
# https://www.datamentor.io/r-programming/subplot/
# https://stackoverflow.com/questions/1249548/side-by-side-plots-with-ggplot2
# how to add unused levels to a legend
# https://stackoverflow.com/questions/40313680/r-is-there-a-way-to-add-unused-data-levels-in-a-ggplot-legend
# fill color options
# https://r-graph-gallery.com/img/graph/42-colors-names.png

?ggplot

# create Pianka's index of niche overlap example
# package 'pgirmess'
# https://rdrr.io/cran/pgirmess/man/piankabioboot.html
library(pgirmess)
data(preybiom)
attach(preybiom)
jackal<-preybiom[site=="Y" & sp=="C",5:6]
genet<-preybiom[site=="Y" & sp=="G",5:6]

piankabio(jackal,genet)
piankabioboot(jackal, genet, B = 1000, probs = c(0.025, 0.975))

# calculate pianka's index for our data (quick and dirty version)
g1 <- PO_W[,4:5]
g2 <- CO_W[,4:5]
# add rows to dataframes to match categories in the PreyItem
g1<-rbind(g1, data.frame(PreyItem = c("bird", "elk", "lagomorph", "med_mammal", "moose", "small_mammal", "other"), Pct_FO = c(0,0,0,0,0,0,0)))
g2<-rbind(g2, data.frame(PreyItem = c("bird", "elk", "carnivore", "lagomorph", "med_mammal", "moose", "small_mammal","whitetaileddeer"), Pct_FO = c(0,0,0,0,0,0,0,0)))
# order prey items
# order legend prey items
g2$PreyItem <- factor(g2$PreyItem, levels=c('bird', 'elk', 'carnivore', 'deerunkspp', 'lagomorph', 'med_mammal', 'moose', 'muledeer', 'small_mammal', 'whitetaileddeer','livestock','other'))
# reorder dataframe
g2 <- with(g2, g2[order(PreyItem),])
# set factor levels 
g1$PreyItem <- factor(g1$PreyItem, levels=c('bird', 'elk', 'carnivore', 'deerunkspp', 'lagomorph', 'med_mammal', 'moose', 'muledeer', 'small_mammal', 'whitetaileddeer','livestock','other'))
g1 <- with(g1, g1[order(PreyItem),])

piankabio(g1, g2)
piankabioboot(coug_ne_summ, coug_ne_wint, B = 1000, probs = c(0.025, 0.975))

# now calculate Shannon's index
data(preybiom)
shannonbio(preybiom[,5:6])
shannonbio(g1)

#######################################
# Sex and Age Data
#######################################
# now let's try to work with the sex and age data for important species groups, which are:
# mule deer for wolves and cougars in the Okanogan
# white-tailed deer for wolves and cougars in the Northeast
# moose for wolves (and possibly cougars) in the Northeast




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
