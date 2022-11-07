# CODE FOR SCAT METABARCODING ANALYSES ---------------------
# data formatting, bootstrapping, %FO plots, NMDS plots, Shannon's Index, Pianka's Index, comparison plots (dietary diversity vs. dietary overlap)

# Clear workspace ------------------------------------------
rm(list = ls())

getwd()

# Species Breakdowns - Scat Data ----------------------------

# read in and format scat data ----------------------------
library(tidyverse)

# load the species found in each scat sample via metabarcoding
data <- read.csv("Data/diets_long_PERMANOVA.csv")
# load the scat log metadata with cluster, carcass, and prey species found at cluster information
scatlog_meta <- read.csv("Data/scatlog_meta_ClusterID_Final.csv")

# create a second database of only wolf and cougar scats with confirmed depositor that also contain vertebrate prey
slm_cc <-
  subset(scatlog_meta,
         Metabarcoding_Species == 'Wolf' |
           Metabarcoding_Species == 'Cougar')
slm_cc <- subset(slm_cc, Metabarcoding_Status == "Pred_w_prey")
# samples sizes for all scats with depositor confirmed regardless of prey contents
with(scatlog_meta, table(Season, StudyArea, Metabarcoding_Species))

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
# check
length(unique(data$SampleID))
length(unique(slm_cc$SampleID))

slm_cc_sampleIDs <- unique(slm_cc$SampleID)
data_sampleIDs <- unique(data$SampleID)
length(slm_cc_sampleIDs) #316
length(data_sampleIDs) #385 
setdiff(slm_cc_sampleIDs, data_sampleIDs) #should be none - character(0)

# sample sizes for Table 1 ----
# sample sizes for amplified scats with prey (subset of previous)
with(slm_cc, table(Season, StudyArea, Metabarcoding_Species))

# sample sizes for Table 2 ----
# check the number of scat found at clusters with prey, without prey, and at non-clusters
with(slm_cc, table(Cluster, CarcassFound, Metabarcoding_Species))

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
length(unique(data$SampleID)) #316

# reassign "deerunkspp" to either "muledeer" or "whitetaileddeer" in the "prey_simple_deerspp" column of the "data" dataframe - assign proportionally within species, study area, and season
## first calculate proportions of MD and WTD species within each group
## then use inside deerloop() function
# MD <-
#   data %>% group_by(Depositor_DNA, StudyArea, Season) %>% summarise(
#     priv_perc = sum(prey_simple_deerspp == "muledeer", na.rm = T) / sum(
#       prey_simple_deerspp == "whitetaileddeer" |
#         prey_simple_deerspp == "muledeer",
#       na.rm = T
#     )
#   )

# WTD <-
#   data %>% group_by(Depositor_DNA, StudyArea, Season) %>% summarise(
#     priv_perc = sum(prey_simple_deerspp == "whitetaileddeer", na.rm = T) / sum(
#       prey_simple_deerspp == "whitetaileddeer" |
#         prey_simple_deerspp == "muledeer",
#       na.rm = T
#     )
#   )

# # create a new column only to indicate rows where "prey_simple_deerspp" = "deerunkspp" and automatically set to "deerunkspp" and NA otherwise
data <- data %>% mutate(
  Status = case_when(
    prey_simple_deerspp == "deerunkspp" ~ "deerunkspp",
    prey_simple_deerspp != "deerunkspp" ~ "known"
  )
)
head(data)

# now select your carnivore/studyarea/season combo
# this is a test calculation to be later used inside a function
# carnivore = "Canis lupus"
# studyarea = "Northeast"
# season = "Winter"
# 
# pct_d <-
#   WTD[WTD$Depositor_DNA == carnivore &
#         WTD$StudyArea == studyarea & WTD$Season == season, 4]
# as.numeric(pct_d)
# 1 - pct_d
# 
# # get number of values of "deerunkspp" in Status, assign n
# nd <- table(data$Status)[1]
# # create a list of random binary variables where pct_d represents the proportion of white-taileddeer, and 1-pct_d represent the number of mule deer
# WTD
# data$random_prop_deer <- data[, 'random_prop_deer'] <- NA
# nd_names <- c()
# head(data)
# nd_names <- sample(
#   c("whitetaileddeer", "muledeer"),
#   size = nd,
#   prob = c(pct_d, (1 - pct_d)),
#   replace = TRUE
# )

# test metrics
carnivore = "Puma concolor"
studyarea = "Northeast"
season = "Summer"

# now create a new column called "rand_prop_deer" that fills in all the "deerunkspp" in the "preys_simple_deerspp" column with either "whitetaileddeer" or "muledeer" based on proportion of other known deer for that carnivore-studyarea-season group
#
#create initial vector
nd_names <- c()
# write "deerloop" function
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
  
  # see if this helps with the factor level error that is halting the for loop that calculates pianka's, %FO, and shannon's
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
    preynames <- factor(preynames)
    data_sub$rand_prop_deer = factor(data_sub$rand_prop_deer, levels = levels(preynames))

  # remove this line witin the loop to just replace deerunkspp, keep to then resample (bootstrap) the entire dataset to a new column
  data_sub$PreyResample <-
    sample(data_sub$rand_prop_deer,
           replace = TRUE)
  return(data_sub)
}

# now create a new column called "rand_prop_deer" that fills in all the "deerunkspp" in the "preys_simple_deerspp" column with either "whitetaileddeer" or "muledeer" based on proportion of other known deer for that carnivore-studyarea-season group
# for this deerloop_cvd option we are looking only at the cervids, so WTD, MD, Elk, and Moose
deerloop_cvd <- function(carnivore, studyarea, season) {
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
  
  # use this line only when calculating for WTD, MD, Elk, Moose, Livestock, and Other - comment out otherwise
  data_sub$rand_prop_deer <- as.factor(data_sub$rand_prop_deer)
  data_sub$rand_prop_deer <- fct_collapse(data_sub$rand_prop_deer, other = c("bird", "lagomorph", "med_mammal", "small_mammal", "carnivore", "elk", "livestock", "other"))

  data_sub <- data_sub[ (data_sub$rand_prop_deer %in% c("muledeer", "whitetaileddeer", "moose") ), ]

  
  # remove this line witin the loop to just replace deerunkspp, keep to then resample (bootstrap) the entire dataset to a new column
  data_sub$PreyResample <-
    sample(data_sub$rand_prop_deer,
           replace = TRUE)
  return(data_sub)
}

#calculate %FO
library(kwb.utils)

# percent FO for all 11 species groups
FO_funct <- function(data, preynames, depositor, studyarea,    season) {
    data %>%
      # this code ensures that all species categories are included in all dataframes even if values for some species are 0 for all samples
      # commenting out this code results in dataframes that only have columns for species detected for the species/studyarea/season group
      mutate(
        rand_prop_deer = fct_relevel(
          PreyResample,
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
      count(SampleID, PreyResample) %>%
      group_by(SampleID) %>%
      mutate(n = prop.table(n)) %>%
      ungroup() %>%
      pivot_wider(
        names_from = PreyResample,
        values_from = n,
        names_prefix = ''
      ) %>%
      replace(is.na(.), 0) %>%
      hsAddMissingCols(., preynames, fill.value = 0) %>%
      mutate_if(is.numeric, ~ 1 * (. > 0)) %>% # make all entries 1 instead of splitting by scat
      select(
        .,
        # reorder so all dataframes have the same column order
        "SampleID",
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
  }

# percent FO for main prey species / cervids only (WTD, MD, Moose)
FO_funct_cvd <- function(data, preynames, depositor, studyarea,  season) {
  data %>%
    # this code ensures that all species categories are included in all dataframes even if values for some species are 0 for all samples
    # commenting out this code results in dataframes that only have columns for species detected for the species/studyarea/season group
      mutate(
        rand_prop_deer = fct_relevel(
          PreyResample,
          "muledeer",
          "whitetaileddeer",
          "moose"
        )
    ) %>%
    filter(Depositor_DNA == depositor,
           StudyArea == studyarea,
           Season == season) %>%
    count(SampleID, PreyResample) %>%
    group_by(SampleID) %>%
    mutate(n = prop.table(n)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = PreyResample,
      values_from = n,
      names_prefix = ''
    ) %>%
    replace(is.na(.), 0) %>%
    hsAddMissingCols(., preynames, fill.value = 0) %>%
    mutate_if(is.numeric, ~ 1 * (. > 0)) %>% # make all entries 1 instead of splitting by scat
    select(
      .,
      # reorder so all dataframes have the same column order
  "SampleID",
  "muledeer",
  "whitetaileddeer",
  "moose")
}

# run 10000 bootstrap iterations with the 11 prey categories-----
# now reassign WTD/MD based on proportion and then bootstrap the dataset once
# do this a set number of times
# create empty vectors for Pianka values
North_Sum_pianka <- c()
North_Wint_pianka <- c()
Okan_Sum_pianka <- c()
Okan_Wint_pianka <- c()
# create empty vectors for Shannon's values
# cougar / Puma concolor / H
PN_S_shan_H <- c()
PN_W_shan_H <- c()
PO_S_shan_H <- c()
PO_W_shan_H <- c()
# cougar / Puma concolor / J
PN_S_shan_J <- c()
PN_W_shan_J <- c()
PO_S_shan_J <- c()
PO_W_shan_J <- c()
# wolf / Canis lupus / H
CN_S_shan_H <- c()
CN_W_shan_H <- c()
CO_S_shan_H <- c()
CO_W_shan_H <- c()
# wolf / Canis lupus / J
CN_S_shan_J <- c()
CN_W_shan_J <- c()
CO_S_shan_J <- c()
CO_W_shan_J <- c()
# create empty vectors for the percent frequency of occurrence 
# first need to assign column names for all 11 prey species
# create an empty dataframe with the correct number of columns
df_spp <- data.frame(matrix(ncol = 11, nrow = 0))
x <- c(
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
colnames(df_spp) <- x

PN_S_FO <- df_spp
PN_W_FO <- df_spp
PO_S_FO <- df_spp
PO_W_FO <- df_spp
CN_S_FO <- df_spp
CN_W_FO <- df_spp
CO_S_FO <- df_spp
CO_W_FO <- df_spp

# COMMENTED OUT THE FOR LOOP - UNCOMMENT TO RE-RUN
for (i in 1:100) {
  # i=1
  # ok now use the function to generate outputs for each  of the eight groups (carnivore, studyarea, season) and then rowbind to merge
  # note that this randomly assigns "whitetailedeeer" or "muledeer" to the "deerunkspp" entries randomly by proportion of that carnivore-studyarea-season subgroup
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

  # now calculate percent frequency of occurrence by scat sample
  
  # having a problem where the function drops prey categories with no entries - testing code below to try to keep/add back zero columns for missing prey categories
  # https://stackoverflow.com/questions/55024338/add-missing-columns-from-different-data-frame-filled-with-0
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
  
  # create dataframes for all combinations that include zero columns for prey categories not recorded
  FO_PN_S <-
    FO_funct(data_new, preynames, 'Puma concolor', 'Northeast', 'Summer')
  FO_PN_W <-
    FO_funct(data_new, preynames, "Puma concolor", "Northeast", "Winter")
  FO_PO_S <-
    FO_funct(data_new, preynames, "Puma concolor", "Okanogan", "Summer")
  FO_PO_W <-
    FO_funct(data_new, preynames, "Puma concolor", "Okanogan", "Winter")
  FO_CN_S <-
    FO_funct(data_new, preynames, "Canis lupus", "Northeast", "Summer")
  FO_CN_W <-
    FO_funct(data_new, preynames, "Canis lupus", "Northeast", "Winter")
  FO_CO_S <-
    FO_funct(data_new, preynames, "Canis lupus", "Okanogan", "Summer")
  FO_CO_W <-
    FO_funct(data_new, preynames, "Canis lupus", "Okanogan", "Winter")
  
  # create a list of the dataframes
  # could save this list below to a .csv for Shannon's index and %FO confidence intervals
  wolfcoug_list <-
    list(FO_PN_S,
         FO_PN_W,
         FO_PO_S,
         FO_PO_W,
         FO_CN_S,
         FO_CN_W,
         FO_CO_S,
         FO_CO_W) %>%
    set_names(
      c(
        'FO_PN_S',
        'FO_PN_W',
        'FO_PO_S',
        'FO_PO_W',
        'FO_CN_S',
        'FO_CN_W',
        'FO_CO_S',
        'FO_CO_W'
      )
    )
  
  library(dplyr)
  library(purrr)
  library(data.table)
  
  # calculating pianka's index using code from above inside a function
  # cougdata <- wolfcoug_list[[1]] #cougar northeast summer
  # wolfdata <- wolfcoug_list[[2]] #wolf northeast summer
  
  piankaind <- function(cougdata, wolfdata) {
    coug_sum <- cougdata %>% select_if(is.numeric) %>% map_dbl(sum)
    coug_sum$PreyItems <- rownames(coug_sum) # make rows columns
    coug_sum <- as.data.frame(coug_sum)
    
    wolf_sum <- wolfdata %>% select_if(is.numeric) %>% map_dbl(sum)
    wolf_sum$PreyItems <- rownames(wolf_sum) # make rows columns
    wolf_sum <- as.data.frame(wolf_sum)
    
    #combine into one dataframe
    P_Index <- rbind(coug_sum, wolf_sum)
    
    # pull vector data for Pianka's index
    coug <- P_Index[1, 1:11]
    wolf <- P_Index[2, 1:11]
    
    # get vectors of PreyItem counts by carnivore
    coug_v <- coug %>% as.numeric(coug)
    wolf_v <- wolf %>% as.numeric(wolf)
    # convert to proportion values
    coug_vp <- coug_v / sum(coug_v)
    wolf_vp <- wolf_v / sum(wolf_v)
    
    # square proportion values
    coug_vp2 <- (coug_vp) ^ 2
    wolf_vp2 <- (wolf_vp) ^ 2
    
    piankaout <-
      sum(coug_vp * wolf_vp) / sqrt(sum(coug_vp2) * sum(wolf_vp2))
    return(piankaout)
  }
  
  
  # output and store final pianka values for each studyarea-season combo
  North_Sum_pianka[i] <-
    piankaind(wolfcoug_list[[1]], wolfcoug_list[[5]]) #0.678 Northeast Summer
  North_Wint_pianka[i] <-
    piankaind(wolfcoug_list[[2]], wolfcoug_list[[6]]) #0.858 Northeast Winter
  Okan_Sum_pianka[i] <-
    piankaind(wolfcoug_list[[3]], wolfcoug_list[[7]]) #0.855 Okanogan Summer
  Okan_Wint_pianka[i] <-
    piankaind(wolfcoug_list[[4]], wolfcoug_list[[8]]) #0.956 Okanogan Winter


  
  # calculate Shannon's diversity index ---------------------
  
  # create dataframes for each species/study area combo
  # all 11 prey categories
  # # cougars / Puma concolor
  PN_S <- wolfcoug_list[[1]]
  PN_W <- wolfcoug_list[[2]]
  PO_S <- wolfcoug_list[[3]]
  PO_W <- wolfcoug_list[[4]]
  
  # wolves / Canis lupus
  CN_S <- wolfcoug_list[[5]]
  CN_W <- wolfcoug_list[[6]]
  CO_S <- wolfcoug_list[[7]]
  CO_W <- wolfcoug_list[[8]]
  
  # # WTD / MD / Moose only (3 prey categories)
  # # cougars / Puma concolor
  # PN_S_O <- wolfcoug_list_other[[1]]
  # PN_W_O <- wolfcoug_list_other[[2]]
  # PO_S_O <- wolfcoug_list_other[[3]]
  # PO_W_O <- wolfcoug_list_other[[4]]
  # 
  # # wolves / Canis lupus
  # CN_S_O <- wolfcoug_list_other[[5]]
  # CN_W_O <- wolfcoug_list_other[[6]]
  # CO_S_O <- wolfcoug_list_other[[7]]
  # CO_W_O <- wolfcoug_list_other[[8]]
  
  # now divide the column sums by the total number of scats (number of rows) to get the percent frequency of occurrence
  # write a function to calculate the %FO to get bootstrapped data prepared for calculation of shannon's index
  FO_spp_shan <-
    function(data_shan,
             depositorspp,
             studyarea,
             season) {
      # # test criteria
      # data_shan <- PN_S
      # depositorspp <- 'Puma concolor'
      # studyarea <- 'Northeast'
      # season <- 'Summer'
      
      # first calcuate
      FO_shan <- subset(data_shan, select = -SampleID)
      FO_spp_shan <- sapply(FO_shan, sum)
      
      FO_spp_shan <- FO_spp_shan / as.integer(nrow(FO_shan))
      FO_spp_shan <- stack(FO_spp_shan)
      FO_spp_shan <- as.data.frame(FO_spp_shan)
      colnames(FO_spp_shan) <- c("Pct_FO", "PreyItem")
      FO_spp_shan <- FO_spp_shan[c("PreyItem", "Pct_FO")]
      FO_spp_shan <-
        FO_spp_shan %>% mutate(StudyArea = studyarea, .before = PreyItem) %>% mutate(DepositorSpp = depositorspp, .before = StudyArea) %>% mutate(Season = season, .after = StudyArea)
      
    }
  
  # now run the shannon filter on all our prey datasets
  # all 11 prey categories
  # # cougars / Puma concolor
  PN_S_shan <-
    FO_spp_shan(PN_S, 'Puma concolor', 'Northeast', 'Summer')
  PN_W_shan <-
    FO_spp_shan(PN_W, 'Puma concolor', 'Northeast', 'Winter')
  PO_S_shan <-
    FO_spp_shan(PO_S, 'Puma concolor', 'Okanogan', 'Summer')
  PO_W_shan <-
    FO_spp_shan(PO_W, 'Puma concolor', 'Okanogan', 'Winter')
  
  # wolves / Canis lupus
  CN_S_shan <- FO_spp_shan(CN_S, 'Canis lupus', 'Northeast', 'Summer')
  CN_W_shan <- FO_spp_shan(CN_W, 'Canis lupus', 'Northeast', 'Winter')
  CO_S_shan <- FO_spp_shan(CO_S, 'Canis lupus', 'Okanogan', 'Summer')
  CO_W_shan <- FO_spp_shan(CO_W, 'Canis lupus', 'Okanogan', 'Winter')
  
  # now save the %FO values for later plotting
  # take dataframe, transpose
  PN_S_FO[i,] <- PN_S_shan$Pct_FO
  PN_W_FO[i,] <- PN_W_shan$Pct_FO
  PO_S_FO[i,] <- PO_S_shan$Pct_FO
  PO_W_FO[i,] <- PO_W_shan$Pct_FO
  CN_S_FO[i,] <- CN_S_shan$Pct_FO
  CN_W_FO[i,] <- CN_W_shan$Pct_FO
  CO_S_FO[i,] <- CO_S_shan$Pct_FO
  CO_W_FO[i,] <- CO_W_shan$Pct_FO
  
  # WTD / MD / Moose only (3 prey categories)
  # cougars / Puma concolor
  # PN_S_O_shan <-
  #   FO_spp_shan(PN_S_O, 'Puma concolor', 'Northeast', 'Summer')
  # PN_W_O_shan <-
  #   FO_spp_shan(PN_W_O, 'Puma concolor', 'Northeast', 'Winter')
  # PO_S_O_shan <-
  #   FO_spp_shan(PO_S_O, 'Puma concolor', 'Okanogan', 'Summer')
  # PO_W_O_shan <-
  #   FO_spp_shan(PO_W_O, 'Puma concolor', 'Okanogan', 'Winter')
  # 
  # # wolves / Canis lupus
  # CN_S_O_shan <-
  #   FO_spp_shan(CN_S_O, 'Canis lupus', 'Northeast', 'Summer')
  # CN_W_O_shan <-
  #   FO_spp_shan(CN_W_O, 'Canis lupus', 'Northeast', 'Winter')
  # CO_S_O_shan <-
  #   FO_spp_shan(CO_S_O, 'Canis lupus', 'Okanogan', 'Summer')
  # CO_W_O_shan <-
  #   FO_spp_shan(CO_W_O, 'Canis lupus', 'Okanogan', 'Winter')
  
  # function "shannon"
  # ? shannon()
  #test shannon's on my data - just need each resource category and the proportion
  library(pgirmess)
  
  # create a function
  shannon_index <- function(data) {
    shan_dat <- data[, 5]
    shan_dat <- as.data.frame(shan_dat)
    shan_dat <- as.numeric(unlist(shan_dat))
    class(shan_dat)
    shannon(shan_dat, base = exp(2))
    
  }
  
  # now use the formatted data to calculate shannon's index
  # all 11 prey species groups
  # cougars / Puma concolor
  PN_S_shan_val <- shannon_index(PN_S_shan)
  PN_W_shan_val <- shannon_index(PN_W_shan)
  PO_S_shan_val <- shannon_index(PO_S_shan)
  PO_W_shan_val <- shannon_index(PO_W_shan)
  
  # wolves / Canis lupus
  CN_S_shan_val <- shannon_index(CN_S_shan)
  CN_W_shan_val <- shannon_index(CN_W_shan)
  CO_S_shan_val <- shannon_index(CO_S_shan)
  CO_W_shan_val <- shannon_index(CO_W_shan)
  
  # WTD/MD/Moose only (3 prey groups)
  # cougars / Puma concolor
  # PN_S_O_shan_val <- shannon_index(PN_S_O_shan)
  # PN_W_O_shan_val <- shannon_index(PN_W_O_shan)
  # PO_S_O_shan_val <- shannon_index(PO_S_O_shan)
  # PO_W_O_shan_val <- shannon_index(PO_W_O_shan)
  # 
  # # wolves / Canis lupus
  # CN_S_O_shan_val <- shannon_index(PN_S_O_shan)
  # CN_W_O_shan_val <- shannon_index(PN_W_O_shan)
  # CO_S_O_shan_val <- shannon_index(PO_S_O_shan)
  # CO_W_O_shan_val <- shannon_index(PO_W_O_shan)
  
  # now get Shannons diversity - H [1] and evenness - J [2] values
  
  # # first create dataframe for the values
  # AllPrey_Shannon <- data.frame()
  # WildPrey_Shannon <- data.frame()
  # 
  # # fill out column values
  # AllPrey_Shannon %>% mutate(
  #   Depositor_DNA = c(
  #     'Puma concolor',
  #     'Puma concolor',
  #     'Puma concolor',
  #     'Puma concolor',
  #     'Canis lupus',
  #     'Canis lupus',
  #     'Canis lupus',
  #     'Canis lupus'
  #   )
  # )
  
  
  # AllPrey_Shannon$StudyArea <-
  #   c(
  #     'Northeast',
  #     'Northeast',
  #     'Okanogan',
  #     'Okanogan',
  #     'Northeast',
  #     'Northeast',
  #     'Okanogan',
  #     'Okanogan'
  #   )
  # AllPrey_Shannon$Season <-
  #   c('Summer',
  #     'Winter',
  #     'Summer',
  #     'Winter',
  #     'Summer',
  #     'Winter',
  #     'Summer',
  #     'Winter')
  
  # all 11 prey species groups
  # cougars / Puma concolor
  # # H
  PN_S_shan_H[i] <- PN_S_shan_val[1]
  PN_W_shan_H[i] <- PN_W_shan_val[1]
  PO_S_shan_H[i] <- PO_S_shan_val[1]
  PO_W_shan_H[i] <- PO_W_shan_val[1]
  # # J
  PN_S_shan_J[i] <- PN_S_shan_val[2]
  PN_W_shan_J[i] <- PN_W_shan_val[2]
  PO_S_shan_J[i] <- PO_S_shan_val[2]
  PO_W_shan_J[i] <- PO_W_shan_val[2]
  
  # wolves / Canis lupus
  # # H
  CN_S_shan_H[i] <- CN_S_shan_val[1]
  CN_W_shan_H[i] <- CN_W_shan_val[1]
  CO_S_shan_H[i] <- CO_S_shan_val[1]
  CO_W_shan_H[i] <- CO_W_shan_val[1]
  # # J
  CN_S_shan_J[i] <- CN_S_shan_val[2]
  CN_W_shan_J[i] <- CN_W_shan_val[2]
  CO_S_shan_J[i] <- CO_S_shan_val[2]
  CO_W_shan_J[i] <- CO_W_shan_val[2]
  
  # WTD/MD/Moose only (3 prey groups)
  # cougars / Puma concolor
  # PN_S_O_shan_val
  # PN_W_O_shan_val
  # PO_S_O_shan_val
  # PO_W_O_shan_val
  # 
  # # wolves / Canis lupus
  # CN_S_O_shan_val
  # CN_W_O_shan_val
  # CO_S_O_shan_val
  # CO_W_O_shan_val
  # 
  print(i)
}

# check on the bootstrapping outputs ----
# bootstrapped pianka value vectors
North_Sum_pianka
North_Wint_pianka
Okan_Sum_pianka
Okan_Wint_pianka

# bootstrapped %FO vectors by species
PN_S_FO
PN_W_FO
PO_S_FO
PO_W_FO
CN_S_FO
CN_W_FO
CO_S_FO
CO_W_FO

# boostrapped shannon's index values (diversity - H and evenness - J)
# cougars
# # H
PN_S_shan_H
PN_W_shan_H
PO_S_shan_H
PO_W_shan_H
# # J
PN_S_shan_J
PN_W_shan_J
PO_S_shan_J
PO_W_shan_J
# wolves
# # H
CN_S_shan_H
CN_W_shan_H
CO_S_shan_H
CO_W_shan_H
# # J
CN_S_shan_J
CN_W_shan_J
CO_S_shan_J
CO_W_shan_J

# save session workspace and outputs
# save.image(file="Oct21.22_bs10000.RData")
# reload session with bootstrapped output -----
# load(file="Oct21.22_bs10000.RData")

# plot Pianka's index, Shannon's index, %FO from bootstrapping ----
# # make a table of the bootstrap results
# drop the _O at the end for pianka's calculated across the 11 prey species groups
# add the _O at the end (eg North_Sum_O) for reduced prey categories (eg WTD, MD, Moose, Livestock)
# plot pianka's index ----
pianka_plot <-
  data.frame(
    group = c(
      'Northeast Summer',
      'Northeast Winter',
      'Okanogan Summer',
      'Okanogan Winter'
    ),
    pianka_mean = c(
      mean(North_Sum_pianka),
      mean(North_Wint_pianka),
      mean(Okan_Sum_pianka),
      mean(Okan_Wint_pianka)
    ),
    # lower = c(
    #   mean(North_Sum_O) - 1.96*sd(North_Sum_O),
    #   mean(North_Wint_O)  - 1.96*sd(North_Wint_O),
    #   mean(Okan_Sum_O) - 1.96*sd(Okan_Sum_O),
    #   mean(Okan_Wint_O) - 1.96*sd(Okan_Wint_O)
    # ),
    # upper = c(
    #   mean(North_Sum_O) + 1.96*sd(North_Sum_O),
    #   mean(North_Wint_O)  + 1.96*sd(North_Wint_O),
    #   mean(Okan_Sum_O) + 1.96*sd(Okan_Sum_O),
    #   mean(Okan_Wint_O) + 1.96*sd(Okan_Wint_O)
    # )
    lower = c(
      quantile(North_Sum_pianka, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(North_Wint_pianka, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(Okan_Sum_pianka, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(Okan_Wint_pianka, probs = c(0, 0.025, 0.975, 1))[2]
      
    ),
    upper = c(
      quantile(North_Sum_pianka, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(North_Wint_pianka, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(Okan_Sum_pianka, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(Okan_Wint_pianka, probs = c(0, 0.025, 0.975, 1))[3]
      
    )
    
  )

    # values for Pianka'a table (Tbl. B1) in supplement ----
        pianka_plot$season <- c("Summer", "Winter", "Summer", "Winter")
        pianka_plot

    # plot for Pianka's results plot (Fig. 4) in paper ----
    # plot the results
    ggplot(pianka_plot, aes(x=group, y=pianka_mean, color=season)) +
      # geom_bar(
      #   aes(x = group, y = pianka_mean),
      #   stat = "identity",
      #   fill = c("firebrick", "skyblue", "firebrick", "skyblue"),
      #   alpha = 0.7
      # ) +
      geom_point(size = 5, color = c("firebrick", "skyblue", "firebrick", "skyblue")) +
      geom_errorbar(
        aes(x = group, ymin = lower, ymax = upper),
        width = 0.4,
        colour = c("firebrick", "skyblue", "firebrick", "skyblue"),
        alpha = 0.9,
        size = 1.3
      ) + 
      theme_bw()+
      theme(axis.text = element_text(size = 15)) +  # x axis subgroups
      theme(axis.title = element_text(size = 15)) + # y axis text
      theme(text = element_text(size = 15)) + # all plot text including numbers
      ggtitle("") + ylab("Pianka's Index (Mean)") + xlab("")


# plot the %FO results from the bootstrap estimates (w/ error bars)
plot(PN_S_FO)
plot(CN_S_FO)
        
############################
# START HERE!
############################        
        
        
        

# plot Shannon's H - diversity index values ----
shannon_H_plot <-
  data.frame(
    group = c(
      'NE Sum (Cougar)',
      'NE Sum (Wolf)',
      'NE Win (Cougar)',
      'NE Win (Wolf)',
      'OK Sum (Cougar)',
      'Ok Sum (Wolf)',
      'OK Win (Cougar)',
      'OK Win (Wolf)'
    ),
    shannon_mean = c(
      mean(PN_S_shan_H),
      mean(CN_S_shan_H),
      mean(PN_W_shan_H),
      mean(CN_W_shan_H),
      mean(PO_S_shan_H),
      mean(CO_S_shan_H),
      mean(PO_W_shan_H),
      mean(CO_W_shan_H)
    ),
    lower = c(
      quantile(PN_S_shan_H, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(CN_S_shan_H, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(PN_W_shan_H, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(CN_W_shan_H, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(PO_S_shan_H, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(CO_S_shan_H, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(PO_W_shan_H, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(CO_W_shan_H, probs = c(0, 0.025, 0.975, 1))[2]
      
    ),
    upper = c(
      quantile(PN_S_shan_H, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(CN_S_shan_H, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(PN_W_shan_H, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(CN_W_shan_H, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(PO_S_shan_H, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(CO_S_shan_H, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(PO_W_shan_H, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(CO_W_shan_H, probs = c(0, 0.025, 0.975, 1))[3]
      
    )
    
  )

# values for Shannon's H table (Tbl. B2) in supplement ----
shannon_H_plot

# plot for Shannon's H results plot (Fig. 5) in paper ----
# plot the results
ggplot(shannon_H_plot) +
  geom_bar(
    aes(x = group, y = shannon_mean),
    stat = "identity",
    fill = c("firebrick", "firebrick", "skyblue", "skyblue", "firebrick", "firebrick", "skyblue", "skyblue"),
    alpha = 0.7
  ) +
  geom_errorbar(
    aes(x = group, ymin = lower, ymax = upper),
    width = 0.4,
    colour = "orange",
    alpha = 0.9,
    size = 1.3
  ) + ggtitle("") + ylab("Shannon's H Index (Mean) \n") + xlab("")

# plot Shannon's J - evenness index values ----
shannon_J_plot <-
  data.frame(
    group = c(
      'NE Sum (Cougar)',
      'NE Sum (Wolf)',
      'NE Win (Cougar)',
      'NE Win (Wolf)',
      'OK Sum (Cougar)',
      'Ok Sum (Wolf)',
      'OK Win (Cougar)',
      'OK Win (Wolf)'
    ),
    shannon_mean = c(
      mean(PN_S_shan_J),
      mean(CN_S_shan_J),
      mean(PN_W_shan_J),
      mean(CN_W_shan_J),
      mean(PO_S_shan_J),
      mean(CO_S_shan_J),
      mean(PO_W_shan_J),
      mean(CO_W_shan_J)
    ),
    lower = c(
      quantile(PN_S_shan_J, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(CN_S_shan_J, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(PN_W_shan_J, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(CN_W_shan_J, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(PO_S_shan_J, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(CO_S_shan_J, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(PO_W_shan_J, probs = c(0, 0.025, 0.975, 1))[2],
      quantile(CO_W_shan_J, probs = c(0, 0.025, 0.975, 1))[2]
      
    ),
    upper = c(
      quantile(PN_S_shan_J, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(CN_S_shan_J, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(PN_W_shan_J, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(CN_W_shan_J, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(PO_S_shan_J, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(CO_S_shan_J, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(PO_W_shan_J, probs = c(0, 0.025, 0.975, 1))[3],
      quantile(CO_W_shan_J, probs = c(0, 0.025, 0.975, 1))[3]
      
    )
    
  )

# values for Shannon's table (XXXX) in ???? ----
shannon_J_plot

# plot for Pianka's results plot (Fig. 4) in paper ----
# plot the results
ggplot(shannon_J_plot) +
  geom_bar(
    aes(x = group, y = shannon_mean),
    stat = "identity",
    fill = c("firebrick", "firebrick", "skyblue", "skyblue", "firebrick", "firebrick", "skyblue", "skyblue"),
    alpha = 0.7
  ) +
  geom_errorbar(
    aes(x = group, ymin = lower, ymax = upper),
    width = 0.4,
    colour = "orange",
    alpha = 0.9,
    size = 1.3
  ) + ggtitle("") + ylab("Shannon's J Index (Mean) \n") + xlab("")


# now open and use the "overlap_diversity_indices.R" code ----
# this is used to make the Shannons vs Pianka's plots for Shannon's H and J
# then continue below

# COMMMENTED OUT BELOW SINCE THIS IS PIANKA'S AND SHANNONS'S FOR WTD, MD, ELK ONLY ----
# # run bootstrap iterations on WTD, MD, Elk, Moose, Livestock and everything else combined into "Other" 
# # now reassign WTD/MD based on proportion and then bootstrap the dataset once
# # do this a set number of times
# # output a vector of pianka values
# North_Sum_O <- c()
# North_Wint_O <- c()
# Okan_Sum_O <- c()
# Okan_Wint_O <- c()
# 
# for (i in 1:25) {
#   # ok now use the function to generate outputs for each  of the eight groups (carnivore, studyarea, season) and then rowbind to merge
#   # wolves
#   out_CN_S <-
#     deerloop(carnivore = "Canis lupus",
#              studyarea = "Northeast",
#              season = "Summer")
#   out_CN_W <-
#     deerloop(carnivore = "Canis lupus",
#              studyarea = "Northeast",
#              season = "Winter")
#   out_CO_S <-
#     deerloop(carnivore = "Canis lupus",
#              studyarea = "Okanogan",
#              season = "Summer")
#   out_CO_W <-
#     deerloop(carnivore = "Canis lupus",
#              studyarea = "Okanogan",
#              season = "Winter")
#   # cougars
#   out_PN_S <-
#     deerloop(carnivore = "Puma concolor",
#              studyarea = "Northeast",
#              season = "Summer")
#   out_PN_W <-
#     deerloop(carnivore = "Puma concolor",
#              studyarea = "Northeast",
#              season = "Winter")
#   out_PO_S <-
#     deerloop(carnivore = "Puma concolor",
#              studyarea = "Okanogan",
#              season = "Summer")
#   out_PO_W <-
#     deerloop(carnivore = "Puma concolor",
#              studyarea = "Okanogan",
#              season = "Winter")
#   
#   # now combine into a new dataframe
#   data_new <-
#     rbind(out_CN_S,
#           out_CN_W,
#           out_CO_S,
#           out_CO_W,
#           out_PN_S,
#           out_PN_W,
#           out_PO_S,
#           out_PO_W)
#   
#   # now calculate percent frequency of occurrence by scat sample
#   
#   # having a problem where the function drops prey categories with no entries - testing code below to try to keep/add back zero columns for missing prey categories
#   # https://stackoverflow.com/questions/55024338/add-missing-columns-from-different-data-frame-filled-with-0
#   # use the FO_funct from earlier to subset data and create columns by species group
#   # first define the prey categories for the function call
#   preynames <- c("muledeer",
#                  "whitetaileddeer",
#                  "moose")
#   
#   # create dataframes for all combinations that include zero columns for prey categories not recorded
#   FO_PN_S <-
#     FO_funct(data_new, preynames, 'Puma concolor', 'Northeast', 'Summer')
#   FO_PN_W <-
#     FO_funct(data_new, preynames, "Puma concolor", "Northeast", "Winter")
#   FO_PO_S <-
#     FO_funct(data_new, preynames, "Puma concolor", "Okanogan", "Summer")
#   FO_PO_W <-
#     FO_funct(data_new, preynames, "Puma concolor", "Okanogan", "Winter")
#   FO_CN_S <-
#     FO_funct(data_new, preynames, "Canis lupus", "Northeast", "Summer")
#   FO_CN_W <-
#     FO_funct(data_new, preynames, "Canis lupus", "Northeast", "Winter")
#   FO_CO_S <-
#     FO_funct(data_new, preynames, "Canis lupus", "Okanogan", "Summer")
#   FO_CO_W <-
#     FO_funct(data_new, preynames, "Canis lupus", "Okanogan", "Winter")
#   
#   # create a list of the dataframes
#   # could save this list below to a .csv for Shannon's index and %FO confidence intervals
#   wolfcoug_list_other <-
#     list(FO_PN_S,
#          FO_PN_W,
#          FO_PO_S,
#          FO_PO_W,
#          FO_CN_S,
#          FO_CN_W,
#          FO_CO_S,
#          FO_CO_W) %>%
#     set_names(
#       c(
#         'FO_PN_S',
#         'FO_PN_W',
#         'FO_PO_S',
#         'FO_PO_W',
#         'FO_CN_S',
#         'FO_CN_W',
#         'FO_CO_S',
#         'FO_CO_W'
#       )
#     )
#   
#   library(dplyr)
#   library(purrr)
#   library(data.table)
#   
#   # calculating pianka's index using code from above inside a function
#   # cougdata <- wolfcoug_list[[1]] #cougar northeast summer
#   # wolfdata <- wolfcoug_list[[2]] #wolf northeast summer
#   #
#   
#   piankaind <- function(cougdata, wolfdata) {
#     coug_sum <- cougdata %>% select_if(is.numeric) %>% map_dbl(sum)
#     coug_sum$PreyItems <- rownames(coug_sum) # make rows columns
#     coug_sum <- as.data.frame(coug_sum)
#     
#     wolf_sum <- wolfdata %>% select_if(is.numeric) %>% map_dbl(sum)
#     wolf_sum$PreyItems <- rownames(wolf_sum) # make rows columns
#     wolf_sum <- as.data.frame(wolf_sum)
#     
#     #combine into one dataframe
#     P_Index <- rbind(coug_sum, wolf_sum)
#     
#     # pull vector data for Pianka's index
#     coug <- P_Index[1, 1:ncol(P_Index)]
#     wolf <- P_Index[2, 1:ncol(P_Index)]
#     
#     
#     # get vectors of PreyItem counts by carnivore
#     coug_v <- coug %>% as.numeric(coug)
#     wolf_v <- wolf %>% as.numeric(wolf)
#     # convert to proportion values
#     coug_vp <- coug_v / sum(coug_v)
#     wolf_vp <- wolf_v / sum(wolf_v)
#     
#     # square proportion values
#     coug_vp2 <- (coug_vp) ^ 2
#     wolf_vp2 <- (wolf_vp) ^ 2
#     
#     piankaout <-
#       sum(coug_vp * wolf_vp) / sqrt(sum(coug_vp2) * sum(wolf_vp2))
#     return(piankaout)
#   }
#   
#   North_Sum_O[i] <-
#     piankaind(wolfcoug_list_other[[1]], wolfcoug_list_other[[5]]) #0.678 Northeast Summer
#   North_Wint_O[i] <-
#     piankaind(wolfcoug_list_other[[2]], wolfcoug_list_other[[6]]) #0.858 Northeast Winter
#   Okan_Sum_O[i] <-
#     piankaind(wolfcoug_list_other[[3]], wolfcoug_list_other[[7]]) #0.855 Okanogan Summer
#   Okan_Wint_O[i] <-
#     piankaind(wolfcoug_list_other[[4]], wolfcoug_list_other[[8]]) #0.956 Okanogan Winter
#   
# }
# 
# # save as a microsoft word .xlsx file with a tab for each scenario 
# # library(xlsx)
# # file <- paste("Output/wolfcoug_list_binaryprey_wtd_md_moose_only.xlsx", sep = "")
# # # write to xlsx
# # write.xlsx(wolfcoug_list_other[[1]], file, sheetName = "FO_PN_S")
# # write.xlsx(wolfcoug_list_other[[2]], file, sheetName = "FO_PN_W", append = TRUE)
# # write.xlsx(wolfcoug_list_other[[3]], file, sheetName = "FO_PO_S", append = TRUE)
# # write.xlsx(wolfcoug_list_other[[4]], file, sheetName = "FO_PO_W", append = TRUE)
# # write.xlsx(wolfcoug_list_other[[5]], file, sheetName = "FO_CN_S", append = TRUE)
# # write.xlsx(wolfcoug_list_other[[6]], file, sheetName = "FO_CN_W", append = TRUE)
# # write.xlsx(wolfcoug_list_other[[7]], file, sheetName = "FO_CO_S", append = TRUE)
# # write.xlsx(wolfcoug_list_other[[8]], file, sheetName = "FO_CO_W", append = TRUE)
# 
# # calculate Pianka's index 
# 
# # make a table of the bootstrap results
# # drop the _O at the end for pianka's calculated across the 11 prey species groups
# # add the _O at the end (eg North_Sum_O) for reduced prey categories (eg WTD, MD, Moose, Livestock)
# pianka_plot <-
#   data.frame(
#     group = c(
#       'Northeast Summer',
#       'Northeast Winter',
#       'Okanogan Summer',
#       'Okanogan Winter'
#     ),
#     pianka_mean = c(
#       mean(North_Sum_pianka),
#       mean(North_Wint_pianka),
#       mean(Okan_Sum_pianka),
#       mean(Okan_Wint_pianka)
#     ),
#     # lower = c(
#     #   conf_int(North_Sum_pianka, 0.95)[1],
#     #   conf_int(North_Wint_pianka, 0.95)[1],
#     #   conf_int(Okan_Sum_pianka, 0.95)[1],
#     #   conf_int(Okan_Wint_pianka, 0.95)[1]
#     #
#     # ),
#     # upper = c(
#     #   conf_int(North_Sum_pianka, 0.95)[2],
#     #   conf_int(North_Wint_pianka, 0.95)[2],
#     #   conf_int(Okan_Sum_pianka, 0.95)[2],
#     #   conf_int(Okan_Wint_pianka, 0.95)[2]
#     #
#     # )
#     # lower = c(
#     #   quantile(North_Sum_pianka, probs = c(0,0.025, 0.975, 1))[2],
#     #   quantile(North_Wint_pianka, probs = c(0,0.025, 0.975, 1))[2],
#     #   quantile(Okan_Sum_pianka, probs = c(0,0.025, 0.975, 1))[2],
#     #   quantile(Okan_Wint_pianka, probs = c(0,0.025, 0.975, 1))[2]
#     #
#     # ),
#     # upper = c(
#     #   quantile(North_Sum_pianka, probs = c(0,0.025, 0.975, 1))[3],
#     #   quantile(North_Wint_pianka, probs = c(0,0.025, 0.975, 1))[3],
#     #   quantile(Okan_Sum_pianka, probs = c(0,0.025, 0.975, 1))[3],
#     #   quantile(Okan_Wint_pianka, probs = c(0,0.025, 0.975, 1))[3]
#     #
#     # )
#     # lower = c(
#     #   mean(North_Sum_O) - 1.96*sd(North_Sum_O),
#     #   mean(North_Wint_O)  - 1.96*sd(North_Wint_O),
#     #   mean(Okan_Sum_O) - 1.96*sd(Okan_Sum_O),
#     #   mean(Okan_Wint_O) - 1.96*sd(Okan_Wint_O)
#     # ),
#     # upper = c(
#     #   mean(North_Sum_O) + 1.96*sd(North_Sum_O),
#     #   mean(North_Wint_O)  + 1.96*sd(North_Wint_O),
#     #   mean(Okan_Sum_O) + 1.96*sd(Okan_Sum_O),
#     #   mean(Okan_Wint_O) + 1.96*sd(Okan_Wint_O)
#     # )
#     lower = c(
#       quantile(North_Sum_pianka, probs = c(0, 0.05, 0.95, 1))[2],
#       quantile(North_Wint_pianka, probs = c(0, 0.05, 0.95, 1))[2],
#       quantile(Okan_Sum_pianka, probs = c(0, 0.05, 0.95, 1))[2],
#       quantile(Okan_Wint_pianka, probs = c(0, 0.05, 0.95, 1))[2]
#       
#     ),
#     upper = c(
#       quantile(North_Sum_pianka, probs = c(0, 0.05, 0.95, 1))[3],
#       quantile(North_Wint_pianka, probs = c(0, 0.05, 0.95, 1))[3],
#       quantile(Okan_Sum_pianka, probs = c(0, 0.05, 0.95, 1))[3],
#       quantile(Okan_Wint_pianka, probs = c(0, 0.05, 0.95, 1))[3]
#       
#     )
#     
#   )
# pianka_plot
# 
# # plot the results
# ggplot(pianka_plot) +
#   geom_bar(
#     aes(x = group, y = pianka_mean),
#     stat = "identity",
#     fill = c("firebrick", "skyblue", "firebrick", "skyblue"),
#     alpha = 0.7
#   ) +
#   geom_errorbar(
#     aes(x = group, ymin = lower, ymax = upper),
#     width = 0.4,
#     colour = "orange",
#     alpha = 0.9,
#     size = 1.3
#   ) + ggtitle("") + ylab("Pianka's Index (Mean)") + xlab("")
# 
# # calculate Shannon's diversity index 
# 
# # create dataframes for each species/study area combo
# # all 11 prey categories
# # # cougars / Puma concolor
# PN_S <- wolfcoug_list[[1]]
# PN_W <- wolfcoug_list[[2]]
# PO_S <- wolfcoug_list[[3]]
# PO_W <- wolfcoug_list[[4]]
# 
# # wolves / Canis lupus
# CN_S <- wolfcoug_list[[5]]
# CN_W <- wolfcoug_list[[6]]
# CO_S <- wolfcoug_list[[7]]
# CO_W <- wolfcoug_list[[8]]
# 
# # WTD / MD / Moose only (3 prey categories)
# # cougars / Puma concolor
# PN_S_O <- wolfcoug_list_other[[1]]
# PN_W_O <- wolfcoug_list_other[[2]]
# PO_S_O <- wolfcoug_list_other[[3]]
# PO_W_O <- wolfcoug_list_other[[4]]
# 
# # wolves / Canis lupus
# CN_S_O <- wolfcoug_list_other[[5]]
# CN_W_O <- wolfcoug_list_other[[6]]
# CO_S_O <- wolfcoug_list_other[[7]]
# CO_W_O <- wolfcoug_list_other[[8]]
# 
# # now divide the column sums by the total number of scats (number of rows) to get the percent frequency of occurrence
# # write a function to calculate the %FO to get bootstrapped data prepared for calculation of shannon's index
# FO_spp_shan <-
#   function(data_shan,
#            depositorspp,
#            studyarea,
#            season) {
#     # test criteria
#     data_shan <- PN_S_O
#     depositorspp <- 'Puma concolor'
#     studyarea <- 'Northeast'
#     season <- 'Summer'
#     
#     # first calcuate
#     FO_shan <- subset (data_shan, select = -SampleID)
#     FO_spp_shan <- sapply(FO_shan, sum)
#     
#     FO_spp_shan <- FO_spp_shan / as.integer(nrow(FO_shan))
#     FO_spp_shan <- stack(FO_spp_shan)
#     FO_spp_shan <- as.data.frame(FO_spp_shan)
#     colnames(FO_spp_shan) <- c("Pct_FO", "PreyItem")
#     FO_spp_shan <- FO_spp_shan[c("PreyItem", "Pct_FO")]
#     FO_spp_shan <-
#       FO_spp_shan %>% mutate(StudyArea = studyarea, .before = PreyItem) %>% mutate(DepositorSpp = depositorspp, .before = StudyArea) %>% mutate(Season = season, .after = StudyArea)
#     
#   }
# 
# # now run the shannon filter on all our prey datasets
# # all 11 prey categories
# # # cougars / Puma concolor
# PN_S_shan <-
#   FO_spp_shan(PN_S, 'Puma concolor', 'Northeast', 'Summer')
# PN_W_shan <-
#   FO_spp_shan(PN_W, 'Puma concolor', 'Northeast', 'Winter')
# PO_S_shan <-
#   FO_spp_shan(PO_S, 'Puma concolor', 'Okanogan', 'Summer')
# PO_W_shan <-
#   FO_spp_shan(PO_W, 'Puma concolor', 'Okanogan', 'Winter')
# 
# # wolves / Canis lupus
# CN_S_shan <- FO_spp_shan(CN_S, 'Canis lupus', 'Northeast', 'Summer')
# CN_W_shan <- FO_spp_shan(CN_W, 'Canis lupus', 'Northeast', 'Winter')
# CO_S_shan <- FO_spp_shan(CO_S, 'Canis lupus', 'Okanogan', 'Summer')
# CO_W_shan <- FO_spp_shan(CO_W, 'Canis lupus', 'Okanogan', 'Winter')
# 
# # WTD / MD / Moose only (3 prey categories)
# # cougars / Puma concolor
# PN_S_O_shan <-
#   FO_spp_shan(PN_S_O, 'Puma concolor', 'Northeast', 'Summer')
# PN_W_O_shan <-
#   FO_spp_shan(PN_W_O, 'Puma concolor', 'Northeast', 'Winter')
# PO_S_O_shan <-
#   FO_spp_shan(PO_S_O, 'Puma concolor', 'Okanogan', 'Summer')
# PO_W_O_shan <-
#   FO_spp_shan(PO_W_O, 'Puma concolor', 'Okanogan', 'Winter')
# 
# # wolves / Canis lupus
# CN_S_O_shan <-
#   FO_spp_shan(CN_S_O, 'Canis lupus', 'Northeast', 'Summer')
# CN_W_O_shan <-
#   FO_spp_shan(CN_W_O, 'Canis lupus', 'Northeast', 'Winter')
# CO_S_O_shan <-
#   FO_spp_shan(CO_S_O, 'Canis lupus', 'Okanogan', 'Summer')
# CO_W_O_shan <-
#   FO_spp_shan(CO_W_O, 'Canis lupus', 'Okanogan', 'Winter')
# 
# # function "shannon"
# ? shannon()
# #test shannon's on my data - just need each resource category and the proportion
# shan_dat <- PN_S_shan[, 5]
# shan_dat <- as.data.frame(shan_dat)
# shan_dat <- as.numeric(unlist(shan_dat))
# class(shan_dat)
# shannon(shan_dat, base = exp(1))
# 
# # create a function
# shannon_index <- function(data) {
#   shan_dat <- data[, 5]
#   shan_dat <- as.data.frame(shan_dat)
#   shan_dat <- as.numeric(unlist(shan_dat))
#   class(shan_dat)
#   shannon(shan_dat, base = exp(2))
#   
# }
# 
# # now use the formatted data to calculate shannon's index
# # all 1 prey species groups
# # cougars / Puma concolor
# PN_S_shan_val <- shannon_index(PN_S_shan)
# PN_W_shan_val <- shannon_index(PN_W_shan)
# PO_S_shan_val <- shannon_index(PO_S_shan)
# PO_W_shan_val <- shannon_index(PO_W_shan)
# 
# # wolves / Canis lupus
# CN_S_shan_val <- shannon_index(CN_S_shan)
# CN_W_shan_val <- shannon_index(CN_W_shan)
# CO_S_shan_val <- shannon_index(CO_S_shan)
# CO_W_shan_val <- shannon_index(CO_W_shan)
# 
# # WTD/MD/Moose only (3 prey groups)
# # cougars / Puma concolor
# PN_S_O_shan_val <- shannon_index(PN_S_O_shan)
# PN_W_O_shan_val <- shannon_index(PN_W_O_shan)
# PO_S_O_shan_val <- shannon_index(PO_S_O_shan)
# PO_W_O_shan_val <- shannon_index(PO_W_O_shan)
# 
# # wolves / Canis lupus
# CN_S_O_shan_val <- shannon_index(PN_S_O_shan)
# CN_W_O_shan_val <- shannon_index(PN_W_O_shan)
# CO_S_O_shan_val <- shannon_index(PO_S_O_shan)
# CO_W_O_shan_val <- shannon_index(PO_W_O_shan)
# 
# # now get Shannons diversity - H [1] and evenness - J [2] values
# 
# # first create dataframe for the values
# AllPrey_Shannon <- data.frame()
# WildPrey_Shannon <- data.frame()
# 
# # fill out column values
# AllPrey_Shannon %>% mutate(
#   Depositor_DNA = c(
#     'Puma concolor',
#     'Puma concolor',
#     'Puma concolor',
#     'Puma concolor',
#     'Canis lupus',
#     'Canis lupus',
#     'Canis lupus',
#     'Canis lupus'
#   )
# )
# 
# 
# AllPrey_Shannon$StudyArea <-
#   c(
#     'Northeast',
#     'Northeast',
#     'Okanogan',
#     'Okanogan',
#     'Northeast',
#     'Northeast',
#     'Okanogan',
#     'Okanogan'
#   )
# AllPrey_Shannon$Season <-
#   c('Summer',
#     'Winter',
#     'Summer',
#     'Winter',
#     'Summer',
#     'Winter',
#     'Summer',
#     'Winter')
# 
# # all 1 prey species groups
# # cougars / Puma concolor
# PN_S_shan_val[1]
# PN_W_shan_val
# PO_S_shan_val
# PO_W_shan_val
# 
# # wolves / Canis lupus
# CN_S_shan_val
# CN_W_shan_val
# CO_S_shan_val
# CO_W_shan_val
# 
# # WTD/MD/Moose only (3 prey groups)
# # cougars / Puma concolor
# PN_S_O_shan_val
# PN_W_O_shan_val
# PO_S_O_shan_val
# PO_W_O_shan_val
# 
# # wolves / Canis lupus
# CN_S_O_shan_val
# CN_W_O_shan_val
# CO_S_O_shan_val
# CO_W_O_shan_val
# 

# # END OF FORMATTED CODE ------------------------------------------------

# # Plot %FO for all species ----
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

# test criteria
depositorspp = "Puma concolor"
studyarea = "Northeast"
season = "Summer"

# create percent frequency of occurrence (Pct_FO) function for single run
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

# create %FO plot for bootstrapped data run

# first look at the data (loaded from bs10000 run above)
# bootstrapped %FO vectors by species
PN_S_FO
PN_W_FO
PO_S_FO
PO_W_FO
CN_S_FO
CN_W_FO
CO_S_FO
CO_W_FO

# test criteria
# data_FO_boot <- PN_S_FO
# depositorspp <- "Puma concolor"
# studyarea <- "Northeast"
# season <- "Summer"

# create a quantile function to use inside the loop to get the lower and upper 95% confidence intervals using the quantile function

quant_boot_low <- function(x) {  # lower 95% CIs
  quantile(as.vector(x), 0.025)
}
quant_boot_high <- function(x) { # upper 95% CIs
  quantile(as.vector(x), 0.975)
}

# take the resulting table and turn it into 1's if the item was present and 0's otherwise
Pct_FO_boot <-
  function(data_FO_boot, depositorspp, studyarea, season) {
    library(dplyr)
    library(tidyr)
        # take the resulting table and turn it into 1's if the item was present and 0's otherwise
    FO_bi <- data_FO_boot %>% mutate_if(is.numeric, ~ 1 * (. > 0))
    FO_spp <- colSums(data_FO_boot)
    FO_low <- sapply(data_FO_boot, quant_boot_low)
    names(FO_low) <- preynames
    FO_high <- sapply(data_FO_boot, quant_boot_high)
    names(FO_high) <- preynames
    # rename columns for FO_low and FO_high from muledeer.97.5%, etc. to just muledeer, etc. using preynames()

    # now divide the column sums by the total number of scats (rows) to get the percent frequency of occurrence
    FO_spp <- FO_spp / as.integer(nrow(FO_bi))
    FO_spp <- stack(FO_spp)
    FO_spp <- as.data.frame(FO_spp)
    colnames(FO_spp) <- c("Pct_FO", "PreyItem")
    FO_spp <- FO_spp[c("PreyItem", "Pct_FO")]
    FO_spp$Pct_low <- FO_low # add lower 95% CI
    FO_spp$Pct_high <- FO_high # add higher 95% CI
    FO_spp <-
      FO_spp %>% mutate(StudyArea = studyarea, .before = PreyItem) %>% mutate(DepositorSpp = depositorspp, .before = StudyArea) %>% mutate(Season = season, .after = StudyArea)
    assign((paste0(
      sub(" .*", "", depositorspp), "_", studyarea, "_", season
    )) , FO_spp)
    
  }

# create dataframes for each species/study area bootstrapped combo

# cougars / Puma concolor
PN_S <- Pct_FO_boot(PN_S_FO, "Puma concolor", "Northeast", "Summer")
PN_W <- Pct_FO_boot(PN_W_FO, "Puma concolor", "Northeast", "Winter")
PO_S <- Pct_FO_boot(PO_S_FO, "Puma concolor", "Okanogan", "Summer")
PO_W <- Pct_FO_boot(PO_W_FO, "Puma concolor", "Okanogan", "Winter")

# wolves / Canis lupus
CN_S <- Pct_FO_boot(CN_S_FO, "Canis lupus", "Northeast", "Summer")
CN_W <- Pct_FO_boot(CN_W_FO, "Canis lupus", "Northeast", "Winter")
CO_S <- Pct_FO_boot(CO_S_FO, "Canis lupus", "Okanogan", "Summer")
CO_W <- Pct_FO_boot(CO_W_FO, "Canis lupus", "Okanogan", "Winter")

# PN <- Pct_FO(data, "Puma concolor", "Northeast")
# PO <- Pct_FO(data, "Puma concolor", "Okanogan")
# CN <- Pct_FO(data, "Canis lupus", "Northeast")
# CO <- Pct_FO(data, "Canis lupus", "Okanogan")

# combine dataframes
FO_spp_all <- rbind(PN_S, PN_W, PO_S, PO_W, CN_S, CN_W, CO_S, CO_W)
# plot
# https://community.rstudio.com/t/help-with-making-plot-with-multiple-columns/50763/2
# by depositor species
# 

# combine dataframes, remove rows (prey items) where Pct_FO is 0, and then and run the code below
# Cougar - Northeast
FO_spp_all <- rbind(PN_S, PN_W) %>% .[.$Pct_FO != 0, ]
# Cougar - Okanogan
FO_spp_all <- rbind(PO_S, PO_W) %>% .[.$Pct_FO != 0, ]
# Wolf - Northeast
FO_spp_all <- rbind(CN_S, CN_W) %>% .[.$Pct_FO != 0, ]
# Wolf - Okanogan
FO_spp_all <- rbind(CO_S, CO_W) %>% .[.$Pct_FO != 0, ]

ggplot(FO_spp_all, aes(DepositorSpp, Pct_FO, fill = PreyItem)) + geom_col(position = "dodge")

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
  ggplot(FO_spp_all, aes(Season, Pct_FO, fill = PreyItem, label = PreyItem)) +
    geom_col(position = position_dodge(width = 0.9, preserve = "single")) +
    geom_errorbar(
      aes(ymin = Pct_low, ymax = Pct_high), 
      position = position_dodge(width = 0.9, preserve = "single"), 
      width = 0.2) +
    scale_fill_manual(
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
p1 <-
  p + scale_y_continuous(limits = c(0, 1.0)) +
    ggtitle("Cougar - Northeast") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("% Frequency of Occurrence") +
    # guides(fill = guide_legend(title = "Prey Item"))
    theme(legend.position = "none")
p2 <-
  p + scale_y_continuous(limits = c(0, 1.0)) +
    ggtitle("Cougar - Okanogan") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("% Frequency of Occurrence") +
    # guides(fill = guide_legend(title = "Prey Item")) 
    theme(legend.position = "none")
p3 <-
  p + scale_y_continuous(limits = c(0, 1.0)) +
    ggtitle("Wolf - Northeast") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("% Frequency of Occurrence") +
    # guides(fill = guide_legend(title = "Prey Item")) 
    theme(legend.position = "none")
p4 <-
  p + scale_y_continuous(limits = c(0, 1.0)) +
    ggtitle("Wolf - Okanogan") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("% Frequency of Occurrence") +
    # guides(fill = guide_legend(title = "Prey Item")) 
    theme(legend.position = "none")
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

p_legend <-
  p + scale_y_continuous(limits = c(0, 1.0)) +
  ggtitle("Cougar - Northeast") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("% Frequency of Occurrence") +
  guides(fill = guide_legend(title = "Prey Item"))
  # theme(legend.position = "none")
# p <- ggplot()
# p2 <- ggplot(mtcars, aes(mpg, wt, col=factor(am))) + geom_point()
legend <- get_legend(p_legend)

#### Percent Frequency of Occurrence Plots ####
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
g2 <- with(g2, g2[order(PreyItem),])
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
g1 <- with(g1, g1[order(PreyItem),])

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
full <-
  vegan::adonis2(forage ~ species + sa + species:sa, data = forage)

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
with(forage, points(ind.nmds$points[species.sa == "Cougar Okanogan",],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Okanogan", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Okanogan", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Cougar Northeast",],
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
with(forage, points(ind.nmds$points[species.sa == "Wolf Okanogan",],
                    pch = 3))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 1]),
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 2])
),
pch = 13, cex = 2))
with(forage, points(ind.nmds$points[species.sa == "Wolf Northeast",],
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
with(forage, points(ind.nmds$points[species.sa == "BlackBear Okanogan",],
                    pch = 3))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "BlackBear Okanogan", 1]),
  mean(ind.nmds$points[species.sa == "BlackBear Okanogan", 2])
),
pch = 13, cex = 2))
with(forage, points(ind.nmds$points[species.sa == "BlackBear Northeast",],
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
with(forage, points(ind.nmds$points[species.sa == "Cougar Northeast",],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Northeast", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Northeast", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Wolf Northeast",],
                    pch = 2))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Wolf Northeast", 1]),
  mean(ind.nmds$points[species.sa == "Wolf Northeast", 2])
),
pch = 14, cex = 2))
with(forage, points(ind.nmds$points[species.sa == "BlackBear Northeast",],
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
with(forage, points(ind.nmds$points[species.sa == "Cougar Okanogan",],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Okanogan", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Okanogan", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Wolf Okanogan",],
                    pch = 2))
with(forage, points(cbind(
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 1]),
  mean(ind.nmds$points[species.sa == "Wolf Okanogan", 2])
),
pch = 14, cex = 2))
with(forage, points(ind.nmds$points[species.sa == "BlackBear Okanogan",],
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
with(forage, points(ind.nmds$points[species.sa == "Cougar Northeast",],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Northeast", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Northeast", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Wolf Northeast",],
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
with(forage, points(ind.nmds$points[species.sa == "Cougar Okanogan",],
                    pch = 3))
with(forage, points(
  cbind(mean(ind.nmds$points[species.sa == "Cougar Okanogan", 1]),
        mean(ind.nmds$points[species.sa == "Cougar Okanogan", 2])),
  pch = 13,
  lwd = 1,
  cex = 2
))
with(forage, points(ind.nmds$points[species.sa == "Wolf Okanogan",],
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
