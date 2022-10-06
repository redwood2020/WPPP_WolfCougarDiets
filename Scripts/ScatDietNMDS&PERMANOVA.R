# Clear workspace
rm(list = ls())

getwd()

`#### Functions
AICc.PERMANOVA <- function(adonis.model) {
  # check to see if object is an adonis model...
  
  if (!(adonis.model$aov.tab[1, 1] >= 1))
    stop("object not output of adonis {vegan} ")
  
  # Ok, now extract appropriate terms from the adonis model
  # Calculating AICc using residual sum of squares (RSS) since I don't think that adonis returns something I can use as a liklihood function...
  
  RSS <-
    adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "SumsOfSqs"]
  MSE <-
    adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "MeanSqs"]
  
  k <-
    ncol(adonis.model$model.matrix)# + 1 # add one for error variance
  
  nn <- nrow(adonis.model$model.matrix)
  
  # AIC : 2*k + n*ln(RSS)
  # AICc: AIC + [2k(k+1)]/(n-k-1)
  
  # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
  # https://www.researchgate.net/post/What_is_the_AIC_formula;
  # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
  
  # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
  # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
  
  AIC <- 2 * k + nn * log(RSS)
  AIC.g <- 2 * k + nn * (1 + log(2 * pi * RSS / nn))
  AIC.MSE <- 2 * k + nn * log(MSE)
  AIC.pi <- k + nn * (1 + log(2 * pi * RSS / (nn - k)))
  AICc <- AIC + (2 * k * (k + 1)) / (nn - k - 1)
  AICc.MSE <- AIC.MSE + (2 * k * (k + 1)) / (nn - k - 1)
  AICc.pi <- AIC.pi + (2 * k * (k + 1)) / (nn - k - 1)
  
  output <- list(
    "AIC" = AIC,
    "AIC.g" = AIC.g,
    "AICc" = AICc,
    "AIC.MSE" = AIC.MSE,
    "AICc.MSE" = AICc.MSE,
    "AIC.pi" = AIC.pi,
    "AICc.pi" = AICc.pi,
    "k" = k
  )
  
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
  
  
  RSS <- null$SumOfSqs[length(null$SumOfSqs) - 1]
  MSE <- RSS / adonis2.model$Df[length(adonis2.model$Df) - 1]
  
  nn <- adonis2.model$Df[length(adonis2.model$Df)] + 1
  
  k <- nn - adonis2.model$Df[length(adonis2.model$Df) - 1]
  
  
  # AIC : 2*k + n*ln(RSS/n)
  # AICc: AIC + [2k(k+1)]/(n-k-1)
  
  # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
  # https://www.statisticshowto.datasciencecentral.com/akaikes-information-criterion/ ;
  # https://www.researchgate.net/post/What_is_the_AIC_formula;
  # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf;
  # https://medium.com/better-programming/data-science-modeling-how-to-use-linear-regression-with-python-fdf6ca5481be
  
  # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
  # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
  
  AIC <- 2 * k + nn * log(RSS / nn)
  AIC.g <- 2 * k + nn * (1 + log(2 * pi * RSS / nn))
  AIC.MSE <- 2 * k + nn * log(MSE)
  AIC.pi <- k + nn * (1 + log(2 * pi * RSS / (nn - k)))
  AICc <- AIC + (2 * k * (k + 1)) / (nn - k - 1)
  AICc.MSE <- AIC.MSE + (2 * k * (k + 1)) / (nn - k - 1)
  AICc.pi <- AIC.pi + (2 * k * (k + 1)) / (nn - k - 1)
  
  output <- list(
    "AIC" = AIC,
    "AICc" = AICc,
    "AIC.g" = AIC.g,
    "AIC.MSE" = AIC.MSE,
    "AICc.MSE" = AICc.MSE,
    "AIC.pi" = AIC.pi,
    "AICc.pi" = AICc.pi,
    "k" = k,
    "N" = nn
  )
  
  return(output)
  
}

################################
# Species Breakdowns - Scat Data
################################

library(tidyverse)

data <- read.csv("Data/diets_long_PERMANOVA.csv")
head(data)
str(data)

# note that I manually added values to the "prey_simple_deerspp" and "prey_simple_unkdeer" columns in Excel to create the diets_long_PERMANOVA.csv from the diets_long.csv - code this in for future (9/7/22)

# add in season, cluster/non-cluster, and carcass found/not found
scatlog <- read.csv("Data/ScatLog_MetabarcodingSamples.csv")
head(scatlog)
str(scatlog)

# create a new column based on the "Found_Along" column named "Cluster" - Yes = found "AtCluster", No = otherwise
table(scatlog$Found_Along)
scatlog <- scatlog %>%
  mutate(
    Cluster = recode(
      Found_Along,
      'AtCluster' = 'Yes',
      'AtCluster, Off-trail' = 'Yes',
      'AtCluster, Road' = 'Yes',
      'AtCluster, Trail' = 'Yes',
      'Follow' = 'No',
      'Off-trail' = 'No',
      'Other' = 'No',
      'Road' = 'No',
      'Trail' = 'No'
    )
  )
with(scatlog, table(Metabarcoding_Species, Cluster))
table(scatlog$Cluster)
head(scatlog)
str(scatlog)

# reclassify "Season" into just 'Summer' (Summer + Fall) and 'Winter' (Winter + Spring)
table(scatlog$Season)
scatlog <- scatlog %>%
  mutate(
    Season = recode(
      Season,
      'Winter' = 'Winter',
      'Spring' = 'Winter',
      'Summer' = 'Summer',
      'Fall' = 'Summer'
    )
  )
with(scatlog, table(Season, Metabarcoding_Species))
table(scatlog$Season)
table(scatlog$Date_Sent_OSU)
head(scatlog)

# create column for carcass found/not found at cluster in the scatlog database so we know how many scats came from clusters with prey (true feeding sites) vs from clusters without prey (presumed resting sites)
# subset the scatlog to only scats that were sent for metabarcoding analysis
scatlog_meta <- scatlog[!(scatlog$Date_Sent_OSU == ""),]
# create a new column for cluster ID in the scat log
scatlog_meta <-
  scatlog_meta %>% mutate(Cluster_ID = gsub('(.*)-\\w+', '\\1', Scat_ID),
                          .after = Scat_ID)
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
scatlog_meta <-
  merge(scatlog_meta, clusters[, c(
    "Cluster_ID",
    "CarcassFound",
    "PreySpecies_Final",
    "Sex_Final",
    "AgeClassDeer_Final",
    "AgeEstimate_Final"
  )], by = "Cluster_ID", all.x = TRUE)
###############################
# NEED TO FIX THE ASSIGNMENT FOR CARCASS FOUND/NOT FOUND FOR NA VALUES, IGNORED FOR NOW B/C DOES NOT AFFECT PIANKA'S INDEX
###############################

with(scatlog_meta,
     table(Metabarcoding_Species, Metabarcoding_Status))
with(scatlog_meta,
     table(CarcassFound, Cluster, Metabarcoding_Species))

# create a second database of only wolf and cougar scats with confirmed depositor that also contain vertebrate prey
slm_cc <-
  subset(scatlog_meta,
         Metabarcoding_Species == 'Wolf' |
           Metabarcoding_Species == 'Cougar')
slm_cc <- subset(slm_cc, Metabarcoding_Status == "Pred_w_prey")
# samples sizes for all scats with depostior confirmed regardless of prey contents
with(scatlog_meta, table(Season, StudyArea, Metabarcoding_Species))
# sample sizes for amplified scats with prey (subset of previous)
with(slm_cc, table(Season, StudyArea, Metabarcoding_Species))

# take the resulting table and turn it into 1's if the item was present and 0's otherwise
# FO_coug_bi <- FO_coug %>% mutate_if(is.numeric, ~ 1 * (. > 0))
# FO_coug_bi <- subset (FO_coug_bi, select = -SampleID)
# FO_coug_spp <- sapply(FO_coug_bi, sum)
# # now divide the column sums by the total number of scats (rows) to get the percent frequency of occurrence
# FO_coug_spp <- FO_coug_spp / as.integer(nrow(FO_coug_bi))
# FO_coug_spp <- stack(FO_coug_spp)
# FO_coug_spp

# ok so next step is to make this into a larger table with columns for "Depositor_DNA", "StudyArea", "PreySpp", "%FO"
# maybe make the above into a function - DepostiorSpecies, StudyArea, PreyCategoryType (ie prey_simple_deerspp or prey_simple_unkdeer)
# first try to turn the %FO code into a function
# note that "preygrouptype" is either "prey_simple_deerspp" or "prey_simple_unkdeer" - the former differentiates between Odocoileus spp and the later groups all deer into "Odocoileus" / deerunkspp

# rename "Metabarcoding_ID" column in the scatlog_meta file to "Sample_ID to match the corresponding column in the "data" file
slm_cc <- rename(slm_cc, SampleID = Metabarcoding_ID)
# add a column for season to the "data" file by pulling from the "slm_cc" file
data <-
  merge(data, slm_cc[, c("SampleID", "Season", "Cluster")], by.x = "SampleID")
# rearrage column order so "Season" comes after "StudyArea"
col_order <-
  c(
    "SampleID",
    "StudyArea",
    "Season",
    "Cluster",
    "Depositor_Field",
    "Depositor_DNA",
    "diet_item",
    "prey_item",
    "prey_simple_deerspp",
    "prey_simple_unkdeer"
  )
data <- data[, col_order]

# attempt to reassign "deerunkspp" to either "muledeer" or "whitetaileddeer" in the "prey_simple_deerspp" column of the "data" dataframe - assign proportionally within species, study area, and season
## first calculate proportions of MD and WTD species within each group
MD <-
  data %>% group_by(Depositor_DNA, StudyArea, Season) %>% summarise(
    priv_perc = sum(prey_simple_deerspp == "muledeer", na.rm = T) / sum(
      prey_simple_deerspp == "whitetaileddeer" |
        prey_simple_deerspp == "muledeer",
      na.rm = T
    )
  )

WTD <-
  data %>% group_by(Depositor_DNA, StudyArea, Season) %>% summarise(
    priv_perc = sum(prey_simple_deerspp == "whitetaileddeer", na.rm = T) / sum(
      prey_simple_deerspp == "whitetaileddeer" |
        prey_simple_deerspp == "muledeer",
      na.rm = T
    )
  )

# create a new column only to indicate rows where "prey_simple_deerspp" = "deerunkspp" and automatically set to "deerunkspp" and NA otherwise
data <- data %>% mutate(
  Status = case_when(
    prey_simple_deerspp == "deerunkspp" ~ "deerunkspp",
    prey_simple_deerspp != "deerunkspp" ~ "known"
  )
)
head(data)

# now select your carnivore/studyarea/season combo
# this is a test calculation to be later used inside a function
carnivore = "Canis lupus"
studyarea = "Northeast"
season = "Winter"

pct_d <-
  WTD[WTD$Depositor_DNA == carnivore &
        WTD$StudyArea == studyarea & WTD$Season == season, 4]
as.numeric(pct_d)
1 - pct_d

# get number of values of "deerunkspp" in Status, assign n
nd <- table(data$Status)[1]
# create a list of random binary variables where pct_d represents the proportion of white-taileddeer, and 1-pct_d represent the number of mule deer
WTD
data$random_prop_deer <- data[, 'random_prop_deer'] <- NA
nd_names <- c()
head(data)
nd_names <- sample(
  c("whitetaileddeer", "muledeer"),
  size = nd,
  prob = c(pct_d, (1 - pct_d)),
  replace = TRUE
)

# now create a new column called "rand_prop_deer" that fills in all the "deerunkspp" in the "preys_simple_deerspp" column with either "whitetaileddeer" or "muledeer" based on proportion of other known deer for that carnivore-studyarea-season group
deerloop <- function(carnivore, studyarea, season) {
  WTD <- data %>% group_by(Depositor_DNA, StudyArea, Season) %>%
    summarise(
      priv_perc = sum(prey_simple_deerspp == "whitetaileddeer", na.rm = T) / sum(
        prey_simple_deerspp == "whitetaileddeer" |
          prey_simple_deerspp == "muledeer",
        na.rm = T
      )
    )
  pct_d <-
    WTD[WTD$Depositor_DNA == carnivore &
          WTD$StudyArea == studyarea & WTD$Season == season, 4]
  data_sub <-
    filter(data,
           Depositor_DNA == carnivore,
           StudyArea == studyarea,
           Season == season)
  
  for (i in 1:length(data_sub$Status)) {
    nd_names[i] <- sample(
      c("whitetaileddeer", "muledeer"),
      size = length(data_sub$Status[i]),
      prob = c(pct_d, (1 - pct_d)),
      replace = TRUE
    )
    data_sub$rand_prop_deer[i] <-
      ifelse(data_sub$Status[i] == "deerunkspp",
             nd_names[i],
             data_sub$prey_simple_deerspp[i])
  }
  return(data_sub)
}

out_CN_W <-
  deerloop(carnivore = "Canis lupus",
           studyarea = "Northeast",
           season = "Winter")
# check that it filtered correctly
with(out_CN_W, table(StudyArea, Season, Depositor_DNA))

# ok now use the function to generate outputs for each  of the eight groups (carnivore, studyarea, season) and then rowbind to merge
# wolves
out_CN_S <-
  deerloop(carnivore = "Canis lupus",
           studyarea = "Northeast",
           season = "Summer")
out_CN_W <-
  deerloop(carnivore = "Canis lupus",
           studyarea = "Northeast",
           season = "Winter")
out_CO_S <-
  deerloop(carnivore = "Canis lupus",
           studyarea = "Okanogan",
           season = "Summer")
out_CO_W <-
  deerloop(carnivore = "Canis lupus",
           studyarea = "Okanogan",
           season = "Winter")
# cougars
out_PN_S <-
  deerloop(carnivore = "Puma concolor",
           studyarea = "Northeast",
           season = "Summer")
out_PN_W <-
  deerloop(carnivore = "Puma concolor",
           studyarea = "Northeast",
           season = "Winter")
out_PO_S <-
  deerloop(carnivore = "Puma concolor",
           studyarea = "Okanogan",
           season = "Summer")
out_PO_W <-
  deerloop(carnivore = "Puma concolor",
           studyarea = "Okanogan",
           season = "Winter")

# now combine into a new dataframe
data_new <-
  rbind(out_CN_S,
        out_CN_W,
        out_CO_S,
        out_CO_W,
        out_PN_S,
        out_PN_W,
        out_PO_S,
        out_PO_W)

# test criteria for workshopping function
depositor = "Puma concolor"
studyarea = "Northeast"
season = "Winter"
str(test)
table(test$rand_prop_deer)

# now calculate percent frequency of occurrence by scat sample
library(kwb.utils)
FO_funct <- function(data, preynames, depositor, studyarea, season) {
  data %>%
    # this code ensures that all species categories are included in all dataframes even if values for some species are 0 for all samples
    # commenting out this code results in dataframes that only have columns for species detected for the species/studyarea/season group
    mutate(
      rand_prop_deer = fct_relevel(
        rand_prop_deer,
        "muledeer",
        "whitetaileddeer",
        "moose",
        "elk",
        "bird",
        "carnivore",
        "lagomorph",
        "med_mammal",
        "small_mammal",
        "other",
        "livestock"
      )
    ) %>%
    filter(Depositor_DNA == depositor,
           StudyArea == studyarea,
           Season == season) %>%
    count(SampleID, rand_prop_deer) %>%
    group_by(SampleID) %>%
    mutate(n = prop.table(n)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = rand_prop_deer,
      values_from = n,
      names_prefix = ''
    ) %>%
    replace(is.na(.), 0) %>%
    hsAddMissingCols(., preynames, fill.value = 0)
}

FO_funct(data_new, preynames, depositor="Puma concolor", studyarea="Northeast", season="Summer")
# having a problem where the function drops prey categories with no entries - testing code below to try to keep/add back zero columns for missing prey categories
# https://stackoverflow.com/questions/55024338/add-missing-columns-from-different-data-frame-filled-with-0
library(kwb.utils)
# use the FO_funct from earlier to subset data and create columns by species group
# first define the prey categories for the function call
preynames <- c(
  "muledeer",
  "whitetaileddeer",
  "moose",
  "elk",
  "bird",
  "carnivore",
  "lagomorph",
  "med_mammal",
  "small_mammal",
  "other",
  "livestock"
)

FO_PN_S <-
  FO_funct(
    data_new,
    preynames,
    depositor == 'Puma concolor',
    studyarea == 'Northeast',
    season == 'Summer'
  )
FO_PN_W <- FO_funct(data_new, preynames, "Puma concolor", "Northeast", "Winter")
FO_PO_S <- FO_funct(data_new, preynames, "Puma concolor", "Okanogan", "Summer")
FO_PO_W <- FO_funct(data_new, preynames, "Puma concolor", "Okanogan", "Winter")
FO_CN_S <- FO_funct(data_new, preynames, "Canis lupus", "Northeast", "Summer")
FO_CN_W <- FO_funct(data_new, preynames, "Canis lupus", "Northeast", "Winter")
FO_CO_S <- FO_funct(data_new, preynames, "Canis lupus", "Okanogan", "Summer")
FO_CO_W <- FO_funct(data_new, preynames, "Canis lupus", "Okanogan", "Winter")

# now use this formatted data to calculate the bootstrapped Pianka's index for each studyarea/season pair
# https://search.r-project.org/CRAN/refmans/spaa/html/niche.overlap.boot.html

library(spaa)

# example of niche.overlap.boot from R documentation
data(datasample)

niche.overlap.boot(datasample[, 1:4], method = "pianka")
head(datasample[, 1:4])

# example of niche.overlap.boot.pair from R documentation
niche.overlap.boot.pair(datasample[, 1], datasample[, 2], method = "levins")

# now for our data

niche.overlap.boot(
  FO_PN_S[, 2],
  method = "pianka",
  times = 100,
  quant = c(0.025, 0.975)
)

# small_mammal
nichedata1 <-
  FO_PN_S$whitetaileddeer # Cougar (Puma concolor) Northeast Summer
nichedata2 <-
  FO_CN_S$whitetaileddeer # Wolf (Canis lupus) Northeast Summer
niche.overlap.boot.pair(nichedata1, nichedata2, method = "pianka")

###########################
# STARTING POINT FOR THURSDAY
###########################
# https://stackoverflow.com/questions/32940822/bootstrapping-data-frame-columns-independently-in-r
# replicate(1000, sapply(example.df, function(x)
metric.function(sample(x, replace = TRUE))))

# Ok what you need to do is create your species/studyarea/season dataframes (with rand_prop_deer filled out) and merge into data_new. **** Save this as an output file so you don't have to keep randomly selecting the "deerunkspp" deer **** Or consider re-running the random deer selection as part of the bootstrapping, but I think we don't want to do this because then we're bootstrapping different datasets essentially ####
# Then replicate the full dataset size (all rows) with replacement however many times ) eg 1000 - make sure to turn that column into a factor first so all columns are retained even if there are zero detections of that species in that data subset
# then get proportion categories for prey species ***of all datasets*** to prepare to calculate pianka's index - do this for wolves and cougars within a studyarea/season
# for each pair (wolf vs cougar within studyarea/season), drop any categories that are zero for both (or try with and without to see if it matters)
# calculate Pianka's index for all the bootstrapped datasets - rows are prey species categories and columns are either wolf or cougar
# Find mean and 95% confidence interval of the set of pianka's indices
# link below explains how to calculate a mean and CI from a simple vector
# https://stackoverflow.com/questions/48612153/how-to-calculate-confidence-intervals-for-a-vector

###########################
data <-
  data_new #reassign to "data" because I'm too lazy to change the varible in the following code
head(data_new)

# try to add a column on the end for proportion of scat (so single items have a value of 1, two prey items would each be 0.5, etc.)
data_test <-
  data_new %>% group_by(SampleID) %>% mutate(ScatProp = 1) %>% ungroup() # 1/n() for proportion
head(data_test)
str(data_test)

shan_filter <-
  function(carnivore, studyarea, season) {
    data_test %>% filter(Depositor_DNA == carnivore,
                         StudyArea == studyarea,
                         Season == season)
  }
shan_filter("Puma concolor", "Northeast", "Summer")
data_test[, 13:14]

# create percent frequency of occurrence (Pct_FO) function
Pct_FO <- function(data, depositorspp, studyarea, season) {
  library(dplyr)
  library(tidyr)
  FO <-
    data %>% dplyr::filter(Depositor_DNA == depositorspp,
                           StudyArea == studyarea,
                           Season == season) %>% count(SampleID, rand_prop_deer) %>% group_by(SampleID) %>% mutate(n = prop.table(n)) %>% ungroup() %>%
    pivot_wider(
      names_from = rand_prop_deer,
      values_from = n,
      names_prefix = ''
    ) %>% replace(is.na(.), 0)
  # take the resulting table and turn it into 1's if the item was present and 0's otherwise
  FO_bi <- FO %>% mutate_if(is.numeric, ~ 1 * (. > 0))
  FO_bi <- subset (FO_bi, select = -SampleID)
  FO_spp <- sapply(FO_bi, sum)
  # now divide the column sums by the total number of scats (rows) to get the percent frequency of occurrence
  FO_spp <- FO_spp / as.integer(nrow(FO_bi))
  FO_spp <- stack(FO_spp)
  FO_spp <- as.data.frame(FO_spp)
  colnames(FO_spp) <- c("Pct_FO", "PreyItem")
  FO_spp <- FO_spp[c("PreyItem", "Pct_FO")]
  FO_spp <-
    FO_spp %>% mutate(StudyArea = studyarea, .before = PreyItem) %>% mutate(DepositorSpp = depositorspp, .before = StudyArea) %>% mutate(Season = season, .after = StudyArea)
  assign((paste0(
    sub(" .*", "", depositorspp), "_", studyarea, "_", season
  )) , FO_spp)
  
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

# combine dataframes
FO_spp_all <- rbind(CO_S, CO_W)

# set factor legend order for PreyItem so they come out the same
str(FO_spp_all$PreyItem)
FO_spp_all$PreyItem <-
  factor(
    FO_spp_all$PreyItem,
    levels = c(
      "muledeer",
      "whitetaileddeer",
      "moose",
      "elk",
      "bird",
      "carnivore",
      "lagomorph",
      "med_mammal",
      "small_mammal",
      "other",
      "livestock"
    )
  )
# by season
p <-
  ggplot(FO_spp_all, aes(Season, Pct_FO, fill = PreyItem, label = PreyItem)) + geom_col(position = position_dodge2(width = 0.9, preserve = "single")) + scale_fill_manual(
    values = c(
      "muledeer" = "chartreuse4",
      "whitetaileddeer" = "dodgerblue",
      "moose" = "orangered",
      "elk" = "orange",
      "bird" = "slateblue2",
      "carnivore" = "seagreen3",
      "lagomorph" = "cadetblue2",
      "med_mammal" = "hotpink2",
      "small_mammal" = "gold",
      "other" = "azure3",
      "livestock" = "black"
    )
  )
# + geom_text(position = position_dodge2(width = 0.9, preserve = "single"), angle = 90, vjust=0.35, hjust= -0.05)

#change plot title ("Cougar - Northeast") based on selected FO_spp_all combined dataframes
p4 <-
  p + scale_y_continuous(limits = c(0, 1.0)) + ggtitle("Wolf - Okanogan") + theme_bw() + theme(plot.title =
                                                                                                 element_text(hjust = 0.5)) + ylab("% Frequency of Occurrence") + guides(fill =
                                                                                                                                                                           guide_legend(title = "Prey Item"))
p1
p2
p3
p4

# create a multiplot with gridarrange()
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
# set up to only have one legend for the 4-panel figure
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x)
    x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# p <- ggplot()
# p2 <- ggplot(mtcars, aes(mpg, wt, col=factor(am))) + geom_point()
legend <- get_legend(p2)


grid.arrange(p1, p2, p3, p4, legend,
             layout_matrix = rbind(c(1, 1, 2, 2, 5), c(3, 3, 4, 4, 5)))


# how to create multiple plots in r
# https://www.datamentor.io/r-programming/subplot/
# https://stackoverflow.com/questions/1249548/side-by-side-plots-with-ggplot2
# how to add unused levels to a legend
# https://stackoverflow.com/questions/40313680/r-is-there-a-way-to-add-unused-data-levels-in-a-ggplot-legend
# fill color options
# https://r-graph-gallery.com/img/graph/42-colors-names.png

? ggplot

# create Pianka's index of niche overlap example
# package 'pgirmess'
# https://rdrr.io/cran/pgirmess/man/piankabioboot.html
library(pgirmess)
data(preybiom)
attach(preybiom)
jackal <- preybiom[site == "Y" & sp == "C", 5:6]
genet <- preybiom[site == "Y" & sp == "G", 5:6]

piankabio(jackal, genet)
piankabioboot(jackal, genet, B = 1000, probs = c(0.025, 0.975))

# calculate pianka's index for our data (quick and dirty version)
g1 <- PN_S[, 4:5]
g2 <- CN_S[, 4:5]
# add rows to dataframes to match categories in the PreyItem
g1 <-
  rbind(g1, data.frame(
    PreyItem = c(
      "bird",
      "elk",
      "lagomorph",
      "med_mammal",
      "moose",
      "small_mammal",
      "other"
    ),
    Pct_FO = c(0, 0, 0, 0, 0, 0, 0)
  ))
g2 <-
  rbind(g2, data.frame(
    PreyItem = c(
      "bird",
      "elk",
      "carnivore",
      "lagomorph",
      "med_mammal",
      "moose",
      "small_mammal",
      "whitetaileddeer"
    ),
    Pct_FO = c(0, 0, 0, 0, 0, 0, 0, 0)
  ))
# order prey items
# order legend prey items
g2$PreyItem <-
  factor(
    g2$PreyItem,
    levels = c(
      'bird',
      'elk',
      'carnivore',
      'deerunkspp',
      'lagomorph',
      'med_mammal',
      'moose',
      'muledeer',
      'small_mammal',
      'whitetaileddeer',
      'livestock',
      'other'
    )
  )
# reorder dataframe
g2 <- with(g2, g2[order(PreyItem), ])
# set factor levels
g1$PreyItem <-
  factor(
    g1$PreyItem,
    levels = c(
      'bird',
      'elk',
      'carnivore',
      'deerunkspp',
      'lagomorph',
      'med_mammal',
      'moose',
      'muledeer',
      'small_mammal',
      'whitetaileddeer',
      'livestock',
      'other'
    )
  )
g1 <- with(g1, g1[order(PreyItem), ])

piankabio(g1, g2)
piankabioboot(g1, g2, B = 1000, probs = c(0.025, 0.975))



# now calculate Shannon's index
data(preybiom)
head(preybiom)
str(preybiom)
shannonbio(preybiom[, 5:6])

# calculate shannon's for our data using the shan_filter function created earlier
dat <- shan_filter("Canis lupus", "Okanogan", "Winter")
dat <- as.data.frame(dat[, 13:14])
head(dat)
str(dat)
shannonbio(dat)

# test bootstrapping
data(preybiom)
myboot <- shannonbioboot(preybiom[, 5:6], B = 100)
library(boot)
boot.ci(myboot,
        index = 1,
        type = c("norm", "basic", "perc")) # confidence intervals for H'
boot.ci(myboot,
        index = 2,
        type = c("norm", "basic", "perc")) # confidence intervals for J'

# now calculate the bootstrapped estimate in the difference between two groups
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/pgirmess/pgirmess.pdf
coug_ne_wint <-
  as.data.frame(shan_filter("Puma concolor", "Northeast", "Winter"))
wolf_ne_wint <-
  as.data.frame(shan_filter("Canis lupus", "Northeast", "Winter"))
coug_ne_wint <- coug_ne_wint[, 13:14]
wolf_ne_wint <- wolf_ne_wint[, 13:14]
head(coug_ne_wint)
head(wolf_ne_wint)

shannonbioboot(coug_ne_wint, B = 100)
difshannonbio(coug_ne_wint,
              wolf_ne_summ,
              R = 1000,
              probs = c(0.025, 0.975))

mydata <- shannonbioboot(dat, B = 100)
library(boot)
boot.ci(mydata,
        index = 1,
        type = c("norm", "basic", "perc")) # confidence intervals for H'
boot.ci(mydata,
        index = 2,
        type = c("norm", "basic", "perc")) # confidence intervals for J'

# function "shannon"
x <- c(0.1, 0.5, 0.2, 0.1, 0.1)
sum(x)
shannon(x)

x <- rpois(10, 6)
shannon(x, base = exp(1))

#test shannon's on my data - just need each resource category and the proportion
shan_dat <- CO_W[, 5]
shan_dat <- as.data.frame(shan_dat)
shan_dat <- as.numeric(unlist(shan_dat))
class(shan_dat)
shannon(shan_dat, base = exp(2))

########################################
# Sex and Age Breakdowns - Cluster Data
########################################
# now let's try to work with the sex and age data for important species groups from the **cluster** data, which are:
# mule deer for wolves and cougars in the Okanogan
# white-tailed deer for wolves and cougars in the Northeast
# moose for wolves (and possibly cougars) in the Northeast

# look at the column headings from the cluster data
head(clusters)

# subset to the specific studyarea-prey group we are interested in
clust_prey <- clusters %>%
  filter(StudySite == "Methow") %>%  # Methow or Northest
  # filter(ClusterType == "Cougar") %>%
  filter(PreySpecies_Final == "muledeer")
# filter(Sex_Final == "Male" | Sex_Final == "Female") %>%
# filter(AgeClassDeer_Simple != "unknown")

head(clust_prey)
with(clust_prey,
     table(AgeClassDeer_Simple, Sex_Final, ClusterType))

with(clust_prey, table(Sex_Final, ClusterType))
with(clust_prey,
     table(AgeClassDeer_Simple, Sex_Final, ClusterType))
table(clust_prey, Sex_Final)
#######################################
## Start from here on Tuesday to make NMDS plots
#######################################
scatlog <- slm_cc
# rename the "Metabarcoding_ID" column to "SampleID" to match other database
scatlog <- rename(scatlog, SampleID = Metabarcoding_ID)
# add in the "Season" column from 'scatlog' to 'data'
# data <- merge(data, scatlog[, c("SampleID", "Season")], by="SampleID")
# add in the "Cluster" column from 'scatlog' to 'data'
# data <- merge(data, scatlog[, c("SampleID", "Cluster")], by="SampleID")
# add in the "Carcass" column from 'scatlog' to 'data'
# data <- merge(data, scatlog[, c("SampleID", "Carcass")], by="SampleID")
head(data)
str(data)

# rearrange columns so "Season" and "Cluster" come after "StudyArea"
# col_order <- c("SampleID", "StudyArea", "Season", "Cluster", "Depositor_Field", "Depositor_DNA", "diet_item", "prey_item", "prey_simple_deerspp", "prey_simple_unkdeer")
# data.scat <- data[, col_order]
data.scat <- data
head(data.scat)
table(data.scat$rand_prop_deer)

# look at some summary tables
with(data.scat, table(Depositor_DNA, Season, StudyArea))

# convert to wide format
library(data.table)
# cover to wide format
data.wide.unkdeer <- data.scat %>%
  group_by(SampleID) %>%
  pivot_wider(names_from = prey_simple_unkdeer, values_from = prey_simple_unkdeer)
data.wide.unkdeer <-
  as.data.table(data.wide.unkdeer) # need to convert to class data.table for this to work
# samples with multiple prey items are still showing up as multiple rows - collapse to one row per SampleID
View(data.wide.unkdeer)

data.wide.deerspp <- data.scat %>%
  group_by(SampleID) %>%
  pivot_wider(names_from = rand_prop_deer, values_from = rand_prop_deer)
data.wide.deerspp <- as.data.table(data.wide.deerspp)

# View(data.scat)
# View(data.wide.unkdeer) # with all deer as "deerunkspp" so we can compare wolf/cougar to bear (which does not have deer species differentiated)
# View(data.wide.deerspp) # with deer species divided into "muledeer", "whitetaileddeer", and "deerunkspp" for the wolf vs cougar comparison

str(data.wide.unkdeer)

#collapse
#unkdeer (where we only have deerunkspp/Odocoileus spp.)
data.wide.unkdeer.factors <- data.wide.unkdeer[, 1:9]
data.wide.unkdeer.binary <- +!is.na(data.wide.unkdeer[, 10:20])
data.wide.unkdeer <-
  cbind(data.wide.unkdeer.factors, data.wide.unkdeer.binary)
setDT(data.wide.unkdeer)
data.wide.unkdeer <-
  data.wide.unkdeer[, .(
    deerunkspp = sum(deerunkspp),
    other = sum(other),
    fish = sum(fish),
    small_mammal = sum(small_mammal),
    bird = sum(bird),
    livestock = sum(livestock),
    elk = sum(elk),
    lagomorph = sum(lagomorph),
    carnivore = sum(carnivore),
    moose = sum(moose),
    med_mammal = sum(med_mammal)
  ),
  by = .(SampleID,
         StudyArea,
         Season,
         Cluster,
         Depositor_Field,
         Depositor_DNA)]# sum rows by SampleID - sum numeric value only
# View(data.wide.unkdeer)

# for the "wide" prey item columns, convert NA values to 0 and non-NA values to 1, so from column 10 and up
#deerspp (where we have muledeer, whitetaileddeer, and deerunkspp / Odocoileus spp.)
data.wide.deerspp.factors <- data.wide.deerspp[, 1:12]
data.wide.deerspp.binary <- +!is.na(data.wide.deerspp[, 13:23])
data.wide.deerspp <-
  cbind(data.wide.deerspp.factors, data.wide.deerspp.binary)
str(data.wide.deerspp)
setDT(data.wide.deerspp)
data.wide.deerspp <-
  data.wide.deerspp[, .(
    livestock = sum(livestock),
    whitetaileddeer = sum(whitetaileddeer),
    moose = sum(moose),
    lagomorph = sum(lagomorph),
    carnivore = sum(carnivore),
    bird = sum(bird),
    med_mammal = sum(med_mammal),
    elk = sum(elk),
    muledeer = sum(muledeer),
    small_mammal = sum(small_mammal),
    other = sum(other)
  ),  by = .(SampleID,
             StudyArea,
             Season,
             Cluster,
             Depositor_Field,
             Depositor_DNA)]

# now convert prey 0/1 values to fractions of the total (so that a scat with both deer and coyote will be 0.5 deer, 0.5 coyote)
#unkdeer
unkdeer.factors <- data.wide.unkdeer[, 1:6]
unkdeer.rows <-
  data.wide.unkdeer[, 7:17] # we dropped the columsn for diet_item, prey_item, and prey_simple_unkdeer so we reduce the column call by a value of 3
unkdeer.rows <- unkdeer.rows / rowSums(unkdeer.rows)
data.wide.unkdeer <- cbind(unkdeer.factors, unkdeer.rows)
# View(data.wide.unkdeer)

#deerspp
deerspp.factors <- data.wide.deerspp[, 1:6]
deerspp.rows <-
  data.wide.deerspp[, 7:17] # we dropped the columsn for diet_item, prey_item, and prey_simple_unkdeer so we reduce the column call by a value of 3
deerspp.rows <- deerspp.rows / rowSums(deerspp.rows)
data.wide.deerspp <- cbind(deerspp.factors, deerspp.rows)
# View(data.wide.deerspp)

# Calculate bootstrapped pianka's index
# https://search.r-project.org/CRAN/refmans/spaa/html/niche.overlap.boot.html
library(spaa)
niche.overlap.boot(
  mat,
  method = c(
    "pianka",
    "schoener",
    "petraitis",
    "czech",
    "morisita",
    "levins"
  ),
  times = 999,
  quant = c(0.025, 0.975)
)
# Example 1
data(datasample)
niche.overlap.boot(datasample[, 1:4], method = "pianka")
niche.overlap.boot(datasample[, 1:4], method = "schoener")
niche.overlap.boot(datasample[, 1:4], method = "czech")
niche.overlap.boot(datasample[, 1:4], method = "levins")
# Example 2
### niche.overlap.boot.pair() example
data(datasample)
niche.overlap.boot.pair(datasample[, 1], datasample[, 2], method = "levins")
# EcosimR
install.packages("EcoSimR")


############################################################
# SCAT DATA NOW READY FOR NMDS PLOTS
############################################################

# first, pick which dataset to work with (unkdeer or deerspp)
data <- data.wide.unkdeer # can change to data.wide.deerspp
data <- data.wide.deerspp
head(data)


# NMDS Plot extras in R: Envfit
# https://jkzorz.github.io/2020/04/04/NMDS-extras.html
# set up for nmds
categ <- data[, 2:6]
ind.log <- data[, 7:17] # 7:17 for unkdeer and deerspp
##### ind.log<-log(forage+1) <- from Shannon's script but I think log(data+1) doesn't make sense here

#convert com to a matrix
m_ind.log = as.matrix(ind.log)


par(ask = FALSE) # set plot
# https://rdrr.io/rforge/vegan/man/metaMDS.html
ind.nmds <-
  vegan::metaMDS(
    m_ind.log,
    "bray",
    k = 2,
    autotransform = F,
    trymax = 100,
    plot = TRUE,
    trace = TRUE,
    previous.best = TRUE,
    na.rm = T
  ) # set trymax to 1000 for real run - shortened to increase processing time

# test plot
plot(ind.nmds)
plot(prey_spp)

# creating nice NMDS plots in R
# https://jkzorz.github.io/2019/06/06/NMDS.html
# https://jkzorz.github.io/2019/06/06/NMDS.html
#extract NMDS scores (x and y coordinates)
library(vegan)

# data.scores has MDS1 and MSD2 column headers
data.scores = as.data.frame(scores(ind.nmds$points))
data.scores$Depositor_DNA <- data$Depositor_DNA
data.scores$Season <- data$Season
data.scores$StudyArea <- data$StudyArea
head(data.scores)


# simple plot
hh <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) +
  geom_point(
    data = data.scores,
    aes(colour = Depositor_DNA),
    size = 3,
    alpha = 0.5
  ) +
  scale_colour_manual(values = c("orange", "steelblue")) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 9, colour = "grey30")
  ) +
  labs(colour = "Depositor Spp.")
hh

hh + geom_segment(
  aes(
    x = 0,
    y = 0,
    xend = NMDS1,
    yend = NMDS2
  ),
  data = prey_spp_coord,
  size = 1,
  alpha = 0.5,
  colour = "grey30"
) +
  geom_text(
    data = prey_spp_coord,
    aes(x = NMDS1, y = NMDS2 + 0.04),
    label = row.names(prey_spp_coord),
    colour = "navy",
    fontface = "bold"
  )

prey_spp <-
  envfit(ind.nmds, ind.log, permutations = 999, na.rm = TRUE)
prey_coord_fact = as.data.frame(scores(prey_spp, "factors")) * ordiArrowMul(prey_spp)
prey_spp_coord = as.data.frame(scores(prey_spp, "vectors")) * ordiArrowMul(prey_spp)
str(data.scores)
str(prey_spp_coord)

# plot
gg <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) +
  geom_point(
    data = data.scores,
    aes(colour = Depositor_DNA),
    size = 3,
    alpha = 0.5
  ) +
  scale_colour_manual(values = c("orange", "steelblue"))  +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = NMDS1,
      yend = NMDS2
    ),
    data = prey_spp_coord,
    size = 1,
    alpha = 0.5,
    colour = "grey30"
  ) +
  geom_point(
    data = prey_spp_coord,
    aes(x = NMDS1, y = NMDS2),
    shape = "diamond",
    size = 4,
    alpha = 0.6,
    colour = "navy"
  ) +
  geom_text(
    data = prey_spp_coord,
    aes(x = NMDS1, y = NMDS2 + 0.04),
    label = row.names(prey_spp_coord),
    colour = "navy",
    fontface = "bold"
  ) +
  geom_text(
    data = prey_spp_coord,
    aes(x = NMDS1, y = NMDS2),
    colour = "grey30",
    fontface = "bold",
    label = row.names(prey_spp_coord)
  ) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 9, colour = "grey30")
  ) +
  labs(colour = "Depositor Species")

gg

# assign variables of interest
species <- data$Depositor_DNA
sa <- data$StudyArea
cluster <- data$Cluster
season <- data$Season
species.sa.season <-
  paste(data$Depositor_DNA, data$StudyArea, data$Season)
sa.season <- paste(data$StudyArea, data$Season)

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

hh <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) +
  geom_point(
    data = data.scores,
    aes(colour = Depositor_DNA),
    size = 3,
    alpha = 0.5
  ) +
  scale_colour_manual(values = c("orange", "steelblue")) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 9, colour = "grey30")
  ) +
  labs(colour = "Depositor Spp.")


n1 <-
  hh + geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = NMDS1,
      yend = NMDS2
    ),
    data = prey_spp_coord,
    size = 1,
    alpha = 0.5,
    colour = "grey30"
  ) +
  geom_text(
    data = prey_spp_coord,
    aes(x = NMDS1, y = NMDS2 + 0.04),
    label = row.names(prey_spp_coord),
    colour = "navy",
    fontface = "bold"
  ) +
  ggtitle("Cougar vs Wolf - Northeast - Summer")

# Cougar vs Wolf - Northeast - Summer
plot(ind.nmds, main =  "Cougar vs Wolf - Northeast - Summer")
vegan::ordihull(
  ind.nmds,
  groups = (species.sa.season == "Puma concolor Northeast Summer"),
  draw = "polygon",
  # col=c("tomato2"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa.season == "Canis lupus Northeast Summer"),
  draw = "polygon",
  # col=c("skyblue2"),
  show.groups = T,
  label = F
)
vegan::orditorp(ind.nmds,
                display = "species",
                col = "black",
                air = 0.001)
with(data, points(cbind(
  mean(ind.nmds$points[species.sa.season == "Puma concolor Northeast Summer", 1]),
  mean(ind.nmds$points[species.sa.season == "Puma concolor Northeast Summer", 2])
),
pch = 14, cex = 2))
with(data, points(cbind(
  mean(ind.nmds$points[species.sa.season == "Canis lupus Northeast Summer", 1]),
  mean(ind.nmds$points[species.sa.season == "Canis lupus Northeast Summer", 2])
),
pch = 12, cex = 2))
with(data, legend(
  "topleft",
  # col=c("skyblue2", "tomato2"),
  legend = c("Puma concolor", "Canis lupus"),
  bty = "n",
  pch = c(3, 2)
))

# Cougar vs Wolf - Northeast - Winter
plot(ind.nmds, main =  "Cougar vs Wolf - Northeast - Winter")
vegan::ordihull(
  ind.nmds,
  groups = (species.sa.season == "Puma concolor Northeast Winter"),
  draw = "polygon",
  # col=c("tomato2"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa.season == "Canis lupus Northeast Winter"),
  draw = "polygon",
  # col=c("skyblue2"),
  show.groups = T,
  label = F
)
vegan::orditorp(ind.nmds,
                display = "species",
                col = "black",
                air = 0.001)
with(data, points(cbind(
  mean(ind.nmds$points[species.sa.season == "Puma concolor Northeast Winter", 1]),
  mean(ind.nmds$points[species.sa.season == "Puma concolor Northeast Winter", 2])
),
pch = 14, cex = 2))
with(data, points(cbind(
  mean(ind.nmds$points[species.sa.season == "Canis lupus Northeast Winter", 1]),
  mean(ind.nmds$points[species.sa.season == "Canis lupus Northeast Winter", 2])
),
pch = 12, cex = 2))
with(data, legend(
  "topleft",
  # col=c("skyblue2", "tomato2"),
  legend = c("Puma concolor", "Canis lupus"),
  bty = "n",
  pch = c(3, 2)
))


# Cougar vs Wolf - Okanogan - Summer
plot(ind.nmds, main =  "Cougar vs Wolf - Okanogan - Summer")
vegan::ordihull(
  ind.nmds,
  groups = (species.sa.season == "Puma concolor Okanogan Summer"),
  draw = "polygon",
  # col=c("tomato2"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa.season == "Canis lupus Okanogan Summer"),
  draw = "polygon",
  # col=c("skyblue2"),
  show.groups = T,
  label = F
)
vegan::orditorp(ind.nmds,
                display = "species",
                col = "black",
                air = 0.001)
with(data, points(cbind(
  mean(ind.nmds$points[species.sa.season == "Puma concolor Okanogan Summer", 1]),
  mean(ind.nmds$points[species.sa.season == "Puma concolor Okanogan Summer", 2])
),
pch = 14, cex = 2))
with(data, points(cbind(
  mean(ind.nmds$points[species.sa.season == "Canis lupus Okanogan Summer", 1]),
  mean(ind.nmds$points[species.sa.season == "Canis lupus Okanogan Summer", 2])
),
pch = 12, cex = 2))
with(data, legend(
  "topleft",
  # col=c("skyblue2", "tomato2"),
  legend = c("Puma concolor", "Canis lupus"),
  bty = "n",
  pch = c(3, 2)
))


# Cougar vs Wolf - Okanogan - Winter
plot(ind.nmds, main =  "Cougar vs Wolf - Okanogan - Winter")
vegan::ordihull(
  ind.nmds,
  groups = (species.sa.season == "Puma concolor Okanogan Winter"),
  draw = "polygon",
  # col=c("tomato2"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa.season == "Canis lupus Okanogan Winter"),
  draw = "polygon",
  # col=c("skyblue2"),
  show.groups = T,
  label = F
)
vegan::orditorp(ind.nmds,
                display = "species",
                col = "black",
                air = 0.001)
with(data, points(cbind(
  mean(ind.nmds$points[species.sa.season == "Puma concolor Okanogan Winter", 1]),
  mean(ind.nmds$points[species.sa.season == "Puma concolor Okanogan Winter", 2])
),
pch = 14, cex = 2))
with(data, points(cbind(
  mean(ind.nmds$points[species.sa.season == "Canis lupus Okanogan Winter", 1]),
  mean(ind.nmds$points[species.sa.season == "Canis lupus Okanogan Winter", 2])
),
pch = 12, cex = 2))
with(data, legend(
  "topleft",
  # col=c("skyblue2", "tomato2"),
  legend = c("Puma concolor", "Canis lupus"),
  bty = "n",
  pch = c(3, 2)
))

# create a multiplot with gridarrange()
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
library(gridExtra)
grid.arrange(n1, n2, n3, n4, nrow = 2, ncol = 2)

# Bear - Okanogan vs Northeast
plot(ind.nmds, main =  "Black Bear")
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Ursus americanus Northeast"),
  draw = "polygon",
  col = c("tomato2"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Ursus americanus Okanogan"),
  draw = "polygon",
  col = c("skyblue2"),
  show.groups = T,
  label = F
)
vegan::orditorp(ind.nmds,
                display = "species",
                col = "black",
                air = 0.01)

plot(vegan::anosim(vegan::vegdist(forage), species.sa, distance = "bray"),
     main = "ANOSIM - Scat Data")
vegan::anosim(vegan::vegdist(forage), species.sa, distance = "bray")

null <- vegan::adonis2(forage ~ 1, data = forage)
speciesmod <- vegan::adonis2(forage ~ species, data = forage)
samod <- vegan::adonis2(forage ~ sa, data = forage)
linmod <- vegan::adonis2(forage ~ species + sa, data = forage)
full <- vegan::adonis2(forage ~ species + sa + species:sa, data = forage)

tabs <- rbind(
  unlist(AICc.PERMANOVA2(null)) ,
  unlist(AICc.PERMANOVA2(speciesmod)),
  unlist(AICc.PERMANOVA2(samod)),
  unlist(AICc.PERMANOVA2(linmod)),
  unlist(AICc.PERMANOVA2(full))
)

tabs

# Northeast - Cougar vs Wolf vs Bear
plot(ind.nmds, main =  "Northeast (Cougar vs. Wolf vs. Black Bear)")
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Wolf Northeast"),
  draw = "polygon",
  col = c("cornflowerblue"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "BlackBear Northeast"),
  draw = "polygon",
  col = c("gold"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Cougar Northeast"),
  draw = "polygon",
  col = c("green3"),
  show.groups = T,
  label = F
)
vegan::orditorp(ind.nmds,
                display = "species",
                col = "black",
                air = 0.01)

# Okanogan - Cougar vs Wolf vs Bear
plot(ind.nmds, main =  "Okanogan (Cougar vs. Wolf vs. Bear)")
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Cougar Okanogan"),
  draw = "polygon",
  col = c("green3"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Wolf Okanogan"),
  draw = "polygon",
  col = c("cornflowerblue"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "BlackBear Okanogan"),
  draw = "polygon",
  col = c("gold"),
  show.groups = T,
  label = F
)
vegan::orditorp(ind.nmds,
                display = "species",
                col = "black",
                air = 0.01)

# Northeast - Cougar vs Wolf
plot(ind.nmds, main =  "Northeast (Cougar vs. Wolf)")
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Wolf Northeast"),
  draw = "polygon",
  col = c("cornflowerblue"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Cougar Northeast"),
  draw = "polygon",
  col = c("green3"),
  show.groups = T,
  label = F
)
vegan::orditorp(ind.nmds,
                display = "species",
                col = "black",
                air = 0.01)

# Okanogan - Cougar vs Wolf
plot(ind.nmds, main =  "Okanogan (Cougar vs. Wolf)")
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Cougar Okanogan"),
  draw = "polygon",
  col = c("green3"),
  show.groups = T,
  label = F
)
vegan::ordihull(
  ind.nmds,
  groups = (species.sa == "Wolf Okanogan"),
  draw = "polygon",
  col = c("cornflowerblue"),
  show.groups = T,
  label = F
)
vegan::orditorp(ind.nmds,
                display = "species",
                col = "black",
                air = 0.01)


# General Plot Call
plot.new()
plot.window(xlim = c(-2.0, 2.0), ylim = c(-2.0, 2.0))
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")


#Cougar - Okanogan vs Northeast
with(forage, points(ind.nmds$points[species.sa == "Cougar Okanogan", ],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Okanogan", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Okanogan", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Cougar Northeast", ],
                    pch = 2))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Cougar Northeast", 1]),
  mean(ind.nmds$points[species.sa == "Cougar Northeast", 2])
),
pch = 14, cex = 2))
op <- par(cex = 1)
with(forage, legend(
  "topright",
  col = c("skyblue2", "tomato2"),
  legend = c("Okanogan", "Northeast"),
  bty = "n",
  pch = c(3, 2)
))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main = "Cougar (Okanogan vs. Northeast)")
box()

#Wolf - Okanogan vs Northeast
with(forage, points(ind.nmds$points[species.sa == "Wolf Okanogan", ],
                    pch = 3))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 1]),
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 2])
),
pch = 13, cex = 2))
with(forage, points(ind.nmds$points[species.sa == "Wolf Northeast", ],
                    pch = 2))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Wolf Northeast", 1]),
  mean(ind.nmds$points[species.sa == "Wolf Northeast", 2])
),
pch = 14, cex = 2))
op <- par(cex = 1)
with(forage, legend(
  "topright",
  col = c("skyblue2", "tomato2"),
  legend = c("Okanognan", "Northeast"),
  bty = "n",
  pch = c(3, 2)
))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main = "Wolf (Okanogan vs. Northeast)")
box()

#Bear - Okanogan vs Northeast
with(forage, points(ind.nmds$points[species.sa == "BlackBear Okanogan", ],
                    pch = 3))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "BlackBear Okanogan", 1]),
  mean(ind.nmds$points[species.sa == "BlackBear Okanogan", 2])
),
pch = 13, cex = 2))
with(forage, points(ind.nmds$points[species.sa == "BlackBear Northeast", ],
                    pch = 2))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "BlackBear Northeast", 1]),
  mean(ind.nmds$points[species.sa == "BlackBear Northeast", 2])
),
pch = 14, cex = 2))
op <- par(cex = 1)
with(forage, legend(
  "topright",
  col = c("skyblue2", "tomato2"),
  legend = c("Okanognan", "Northeast"),
  bty = "n",
  pch = c(3, 2)
))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main = "Black Bear (Okanogan vs. Northeast)")
box()

#Northeast - Cougar vs Wolf vs Bear
with(forage, points(ind.nmds$points[species.sa == "Cougar Northeast", ],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Northeast", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Northeast", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Wolf Northeast", ],
                    pch = 2))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Wolf Northeast", 1]),
  mean(ind.nmds$points[species.sa == "Wolf Northeast", 2])
),
pch = 14, cex = 2))
with(forage, points(ind.nmds$points[species.sa == "BlackBear Northeast", ],
                    pch = 0))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "BlackBear Northeast", 1]),
  mean(ind.nmds$points[species.sa == "BlackBear Northeast", 2])
),
pch = 12, cex = 2))
op <- par(cex = 1)
with(forage, legend(
  "topright",
  col = c("green3", "cornflowerblue", "gold"),
  legend = c("Cougar", "Wolf", "Black Bear"),
  bty = "n",
  pch = c(3, 2, 0)
))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main = "Northeast (Cougar vs. Wolf vs. Black Bear)")
box()

#Okanogan - Cougar vs Wolf vs Black Bear
with(forage, points(ind.nmds$points[species.sa == "Cougar Okanogan", ],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Okanogan", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Okanogan", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Wolf Okanogan", ],
                    pch = 2))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 1]),
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 2])
),
pch = 14, cex = 2))
with(forage, points(ind.nmds$points[species.sa == "BlackBear Okanogan", ],
                    pch = 0))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "BlackBear Okanogan", 1]),
  mean(ind.nmds$points[species.sa == "BlackBear Okanogan", 2])
),
pch = 12, cex = 2))
op <- par(cex = 1)
with(forage, legend(
  "topright",
  col = c("green3", "cornflowerblue", "gold"),
  legend = c("Cougar", "Wolf", "Black Bear"),
  bty = "n",
  pch = c(3, 2, 0)
))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main = "Okanogan (Cougar vs. Wolf vs. Black Bear)")
box()

#Northeast - Cougar vs Wolf
with(forage, points(ind.nmds$points[species.sa == "Cougar Northeast", ],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Northeast", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Northeast", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Wolf Northeast", ],
                    pch = 2))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Wolf Northeast", 1]),
  mean(ind.nmds$points[species.sa == "Wolf Northeast", 2])
),
pch = 14, cex = 2))
op <- par(cex = 1)
with(forage, legend(
  "topright",
  col = c("green3", "cornflowerblue"),
  legend = c("Cougar", "Wolf"),
  bty = "n",
  pch = c(3, 2)
))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main = "Northeast (Cougar vs. Wolf)")
box()

#Okanogan - Cougar vs Wolf
with(forage, points(ind.nmds$points[species.sa == "Cougar Okanogan", ],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Okanogan", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Okanogan", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Wolf Okanogan", ],
                    pch = 2))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 1]),
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 2])
),
pch = 14, cex = 2))
op <- par(cex = 1)
with(forage, legend(
  "topright",
  col = c("green3", "cornflowerblue"),
  legend = c("Cougar", "Wolf"),
  bty = "n",
  pch = c(3, 2)
))
axis(side = 1)
axis(side = 2)
title(xlab = "NMDS 1", ylab = "NMDS 2", main = "Okanogan (Cougar vs. Wolf)")
box()
