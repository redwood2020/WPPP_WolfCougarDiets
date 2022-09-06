# =============================================================
# Metabarcoding data cleaning and merger with WTD/MD analysis
# Washington Predator-Prey Project
# Lauren Satterfield
# August 2022
# =============================================================
# Script to clean wolf, cougar, and black bear metabarcoding results by condensing all the strings of genetic code down to only those with high confidence matches in Genebank. These metabarcoding genetic results include results from a follow-up PCR to differentiate between wolf and coyote run by Jennifer Allen at the Levi Lab at OSU). Note that vertebrate metabarcoding is highly sensitive but often can only classify down to genus (as is the case with Canis and Odocoileus), and thus the follow-up PCR analyses are needed. Once cleaned, the cleaned dataset is merged with results from Ellie Reese (UW SEFS Genetics Lab) to differentiate between white-tailed deer and mule deer. WTD vs MD analysis results exist for some but not all samples, as some samples did not have enough genetic material for this species-level determination, in which case they were left as "Odocoileus spp.".
# ===============================================================

# clear workspace
rm(list = ls());gc()

# see working directory
getwd()

# load libraries
library(tidyr)
library(dplyr)
library(readr)
library(tidyverse)
library(data.table)
library(reshape2)
library(splitstackshape)

# read in the wolf, cougar, and black bear scat metabarcoding results
meta <- read.csv('../Data/IlluminaRun23_Satterfield_UW__results_Fall2021_v2_wCanisResults.csv', header=T)
meta <- as.data.table(meta)
head(meta)

# # look at contents by scat sample
# diets <- dcast(meta, SampleID ~ Scientific_Name_Simple)
# diets <- as.data.frame(diets)
# head(diets)
# str(diets)
# 
# # convert to binary data
# diet_t <- diets %>% 
#   mutate_at(c(2:48), ~replace(., .!=0, 1)) 

# # filter out controls
"%!in%" = Negate("%in%") # create function to filter out list of controls
meta_nc <- meta %>%
  filter(SampleID %!in% c("NTC22", "NTC23", "OBB2", "OWB1","OWB2"))

 # subset to unique rows by prey item within each scat, create columns for % query cover, % identical matches, relative read count, total read count, and percent of read count
# % query cover - is all of the 12S locus there or not?
# % identical - how much does the DNA match the sequence available
# relative read count - when read count is <0.5-1% of sample, you throw those items out as unreliable 
# we want the row with a) highest % query cover, then b) highest % identity, then c) check the relative read count against the total read count
# read more here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3867762/

# ERASE THIS WHEN DONE

# subset the metabarcoding results to unique prey items in each scat and select the rows with the highest query coverage
diets <- meta_nc %>%
  group_by(SampleID, Scientific_Name) %>%
  slice(which.max(Query_Coverage))

# subset the metabarcoding results to unique prey items in each scat and select the rows with the highest query coverage
diets <- meta_nc %>%
  group_by(SampleID, Scientific_Name, Rep) %>%
  slice(which.max(Query_Coverage)) %>%
  mutate(sum.rep = sum(Counts)) 
View(diets)

# now create a dataframe of sampleID and total sum of read counts per sample/ diet item / rep
diets.tot.rep <- meta_nc %>%
  group_by(SampleID) %>%
  mutate(sum.tot = sum(Counts)) 

diets.tot.rep <- diets.tot.rep %>%
  group_by(SampleID) %>% 
  slice(which.max(sum.tot)) %>%
  select(SampleID, sum.tot)
View(diets.tot.rep)
  
# merge dataframes to create new column for total read count per sample and percent of total read count per rep a, b, and c
diets <- merge(diets, diets.tot.rep, by = "SampleID")

# create new column to calculate the percent of total read counts represented by each rep count
diets <- diets %>%
  mutate(rep.pct.tot = sum.rep/sum.tot)
View(diets)

# now look for instances where the read counts percents (rep.pct.tot) are <0.005 or <0.01 (less than 0.5% or less than 1%) of total read count (sum.tot)
low.reads <- diets[which(diets$rep.pct.tot <= 0.01),]
View(low.reads)

# left off here on 8/31/2022
#########################

# subset to data where Query Coverage and % Identical Matches are less than 100

diets_sub <- diets %>%
  filter(Percent_Identical_Matches<100)

table(meta_nc$Percent_Identical_Matches) 
table(diets$Percent_Identical_Matches) 
head(meta_nc[ which(meta_nc$Percent_Identical_Matches < 91), ])

table(meta_nc$Query_Coverage) 
table(diets$Query_Coverage) 


# figure out which entries are a 100% match in both Percent_Identical_Matches and Query_Coverage
diets_100 <- diets %>% distinct(SampleID, Scientific_Name_Simple, .keep_all = TRUE)
# find rows in full dataset (meta_nc) that are not in the "100%" dataset (diets_100)
# https://stackoverflow.com/questions/3171426/compare-two-data-frames-to-find-the-rows-in-data-frame-1-that-are-not-present-in
questionable <- anti_join(meta_nc, diets_100)

diets <- as.data.frame(diets)
head(diets)
str(diets)

diets <- diets %>%
  group_by(SampleID) %>%
  # Create string listing all items in given Id, separated by comma
  summarise(Items = str_c(Scientific_Name_Simple, collapse = ', '))

head(diets)

# look at results by species
table(meta_nc$Scientific_Name_Simple, meta$Depositor_Field)

# remove spaces before diet items
diets <- separate(diets, 'Items', paste("Diet_Item", 1:6, sep="_"), sep=",", extra="drop")
diets$Diet_Item_2 <- trimws(diets$Diet_Item_2)
diets$Diet_Item_3 <- trimws(diets$Diet_Item_3)
diets$Diet_Item_4 <- trimws(diets$Diet_Item_4)
diets$Diet_Item_5 <- trimws(diets$Diet_Item_5)
diets$Diet_Item_6 <- trimws(diets$Diet_Item_6)
diets <- as.data.frame(diets)
diets

# Create fields for Study Area and Depositor from FieldID, plus a column of NAs for Depositor from DNA
diets <- diets %>%
  mutate(StudyArea=substr(SampleID,1,1),.after = SampleID) %>%
  mutate(Depositor_Field=substr(SampleID,2,2),.after = StudyArea) %>%
  add_column(Depositor_DNA = "NA",.after = "Depositor_Field")
  
diets <- data.table(diets) %>%
  .[StudyArea == "N", StudyArea := "Northeast"] %>%
  .[StudyArea == "O", StudyArea := "Okanogan"] %>%
  .[Depositor_Field == "B", Depositor_Field := "Blackbear"] %>%
  .[Depositor_Field == "W", Depositor_Field := "Wolf"] %>%
  .[Depositor_Field == "C", Depositor_Field := "Cougar"] 

# check the output
head(diets)
tail(diets)

# check to see if last three rows are necessary - they are all empty so remove them
# table(diets$Diet_Item_7)
# table(diets$Diet_Item_8)
# table(diets$Diet_Item_9)
# diets <- diets[,1:10]

# populate Depositor_DNA with carnivore if present in same row

for (i in 1:length(diets$SampleID)) {
  temp <- as.matrix(diets[i])[1,]
  if ("Canis lupus" %in% temp) {
    diets$Depositor_DNA[i] <- "Canis lupus"
  } else if ("Puma concolor" %in% temp) {
    diets$Depositor_DNA[i] <- "Puma concolor"
  } else if ("Ursus americanus" %in% temp) {
    diets$Depositor_DNA[i] <- "Ursus americanus"
  } else if ("Canis latrans" %in% temp) {
    diets$Depositor_DNA[i] <- "Canis latrans"
  } else if ("Lynx rufus" %in% temp) {
    diets$Depositor_DNA[i] <- "Lynx rufus"
  } else if ("Canis" %in% temp) {
    diets$Depositor_DNA[i] <- "Canis spp."
  } else {
    diets$Depositor_DNA[i] <- "Other"
  }
}

# check it
head(diets)
tail(diets)
table(diets$Depositor_DNA)
table(diets$Depositor_Field)

# now remove values from the Diet_Item_X columns if they match the Depostitor_DNA columns
for (i in 1:length(diets$SampleID)) {
  temp <- as.vector(diets[i])[1,]
  print_val <- NULL
  if (temp$Depositor_DNA %in% temp[1,5:10]) {
    col_val <- match(temp[1,5:10], temp$Depositor_DNA)
    col_val <- as.numeric(which(col_val==1) + 4) # add 4 to account for the four columns we ignore
    print_val <- as.character(diets[i, ..col_val])
    diets[i, col_val] <- NA
  } else {
    print(cat(diets$SampleID[i], diets$Depositor_Field[i], diets$Depositor_DNA[i], print_val, "skip", " "), sep=" ")
  }
}

# look at the resulting data
head(diets)
# View(diets)
diets <- as.data.frame(diets)
str(diets)

# look at some summary data tables
# https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
with(diets, table(StudyArea, Depositor_DNA)) 

# shift missing diet values to the left
# https://stackoverflow.com/questions/49079789/using-r-to-shift-values-to-the-left-of-data-frame
diets[] <-  t(apply(diets, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
head(diets)
# View(diets)

# write to .csv
# write.csv(diets, "diets.csv", row.names = FALSE)

# Now merge with WTD vs MD analysis by Ellie Reese to get Odocoileus down to species
# Note that I subset Ellie's results to only those with Med and High confidence 
# To make the code simpler and made the matching column "SampleID"
# First, read in Ellie's data, subset to only Med and High confidence results
# https://stackoverflow.com/questions/1299871/how-to-join-merge-data-frames-inner-outer-left-right
deer_results <- read.csv("Ellie_ScatDeerResults_Preliminary_March7.2022.csv", header=T)
head(deer_results)

#merge with the diets results
diets_deer <- merge(x = diets, y = deer_results, by = "SampleID", all = TRUE)
head(diets_deer)

# subset just columns from diets_long plus Species.ID from Ellie's data
diets_deer <- diets_deer[c("SampleID", "StudyArea", "Depositor_Field", "Depositor_DNA","Diet_Item_1", "Diet_Item_2", "Diet_Item_3", "Diet_Item_4", "Diet_Item_5", "Diet_Item_6", "Species.ID")]
head(diets_deer)
table(diets_deer$Species.ID)

# replace MD and WTD with scientific names in the Species.ID column
diets_deer$Species.ID[diets_deer$Species.ID == "MD"] <- "Odocoileus hemionus"
diets_deer$Species.ID[diets_deer$Species.ID == "WTD"] <- "Odocoileus virginianus"
View(diets_deer)
head(diets_deer)

# now we want to replace all occurrences of Odocoileus spp. with the specific species
# from the Species.ID column where it exists, but do not change it otherwise

for (i in 1:length(diets_deer$SampleID)) {
  if ("Odocoileus spp." %in% diets_deer$Diet_Item_1[i]) {
    diets_deer$Diet_Item_1[i] <- diets_deer$Species.ID[i]
  } else if ("Odocoileus spp." %in% diets_deer$Diet_Item_2[i]) {
    diets_deer$Diet_Item_2[i] <- diets_deer$Species.ID[i]
  } else if ("Odocoileus spp." %in% diets_deer$Diet_Item_3[i]) {
    diets_deer$Diet_Item_3[i] <- diets_deer$Species.ID[i]
  } else if ("Odocoileus spp." %in% diets_deer$Diet_Item_4[i]) {
    diets_deer$Diet_Item_4[i] <- diets_deer$Species.ID[i]
  } else if ("Odocoileus spp." %in% diets_deer$Diet_Item_5[i]) {
    diets_deer$Diet_Item_5[i] <- diets_deer$Species.ID[i]
  } else if ("Odocoileus spp." %in% diets_deer$Diet_Item_6[i]) {
    diets_deer$Diet_Item_6[i] <- diets_deer$Species.ID[i]
  } else {
      print(cat(diets_deer$SampleID[i], diets_deer$Species.ID[i], print_val, "skip", " "), sep=" ")
  }
}
View(diets_deer)

# shift missing diet values to the left again
# https://stackoverflow.com/questions/49079789/using-r-to-shift-values-to-the-left-of-data-frame
diets_deer[] <-  t(apply(diets_deer, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
head(diets_deer)
View(diets_deer)

# check for Diet_Item columns that are all NA and remove if so (don't do this earlier 
# though because we have to move the carnivore depositor ID over, so need at least 5
# diet item columns since the most in one scat is four items plus depositor)
table(diets_deer$Diet_Item_1)
table(diets_deer$Diet_Item_2)
table(diets_deer$Diet_Item_3)
table(diets_deer$Diet_Item_4) # only two samples have four diet items
table(diets_deer$Diet_Item_5) # none have 5 items
table(diets_deer$Diet_Item_6) # none have 6 items

diets_deer <- diets_deer[,1:8]
head(diets_deer)

# write.csv(diets_deer, "diets_deer.csv", row.names = FALSE)

# now create it in long format
colnames(diets_deer)
diets_long <- gather(diets_deer, diet_item, prey_item, Diet_Item_1:Diet_Item_4, factor_key=TRUE) # long format
diets_long <- diets_long %>% drop_na(prey_item)
with(diets_long, diets_long[order(SampleID, prey_item),]) # sort by SampleID, prey_item
tail(diets_long, 30)

# write the long format .csv file
write.csv(diets_long, "diets_long_deer.csv", row.names = FALSE)



