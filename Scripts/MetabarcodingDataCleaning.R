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
meta <- read.csv('./Data/IlluminaRun23_Satterfield_UW__results_Fall2021_v2_wCanisResults.csv', header=T)
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
# View(diets)

# subset the metabarcoding results to unique prey items in each scat and select the rows with the highest query coverage
# diets <- diets %>%
#   mutate(sum.rep = sum(Counts)) 
# View(diets)

# now create a dataframe of sampleID and total sum of read counts per sample/ diet item / rep
# diets.tot.rep <- meta_nc %>%
#   group_by(SampleID) %>%
#   mutate(sum.tot = sum(Counts))
# 
# diets.tot.rep <- diets.tot.rep %>%
#   group_by(SampleID) %>% 
#   slice(which.max(sum.tot)) %>%
#   select(SampleID, sum.tot)
# View(diets.tot.rep)

# merge dataframes to create new column for total read count per sample and percent of total read count per rep a, b, and c
# diets <- merge(diets, diets.tot.rep, by = "SampleID")
# View(diets)

# create new column to calculate the percent of total read counts represented by each rep count
# diets <- diets %>%
#   mutate(rep.pct.tot = sum.rep/sum.tot)
# View(diets)

# now look for instances where the read counts percents (rep.pct.tot) are <0.005 or <0.01 (less than 0.5% or less than 1%) of total read count (sum.tot)
# low.reads <- diets[which(diets$rep.pct.tot <= 0.01),]
# View(low.reads)

# now look for instances where only one species' dna showed up in the scat, meaning that it was "empty" (depositor dna found but no prey dna) or "prey only" (no depositor dna, only prey dna)
unusable <- diets %>%
  group_by(SampleID) %>%
  filter(length(unique(Scientific_Name)) == 1)
# View(unusable)
length(unique(unusable$SampleID)) #160 (so 404 usable)
length(unique(diets$SampleID)) #564 amplified from the total of 606

# remove rows from full dataframe that match with sample IDs in the "unusable" dataframe which don't contain dna from at least one depositor and one prey dna
usable <- diets[!(diets$SampleID %in% unusable$SampleID),]
View(usable)
length(unique(usable$SampleID)) #404

# calculate success rates - 606 samples sent in all
564/606 # 93.1% - percent amplified 
404/606 # 66.7% percent both depostor and prey dna (no empty (depostior only) or prey only)


# left off here on 8/31/2022
#########################

# read in the wolf, cougar, and black bear scat metabarcoding results from Ellie looking at WTD vs MD species
scat.deer <- read.csv('./Data/Scat Metabarcoding - Deer Genetics Results.csv', header=T)
scat.deer <- as.data.table(scat.deer)
head(scat.deer)

# subset scat where the "Confidence" value is "High" or "Medium" -> good enough to include in analysis 
scat.deer <- scat.deer %>%
  filter(Confidence == "High" | Confidence == "Med")
levels(as.factor(scat.deer$Confidence))
# View(scat.deer)
length(unique(scat.deer$Sample))
head(scat.deer)

# rename "Sample" column to "SampleID" to match other database
scat.deer <- rename(scat.deer, SampleID = Sample)

# rename MD (mule deer) to "Odocoileus hemionus" and WTD to "Odocoileus virginianus"
scat.deer$Species.ID <- recode_factor(scat.deer$Species.ID, 'MD' = "Odocoileus hemionus", 'WTD' = "Odocoileus virginianus") # recode all to Odocoileus

# in the usable database, remove a trailing "A" from any samples with the "A" listed at the end - needed so IDs will match with the scat.deer database
usable$SampleID <- sub("A$", "", usable$SampleID)
View(usable)

# now replace "Odocoileus spp." column with appropriate confirmed deer spp from the scat.deer database
# first under replace all "Odocoileus (HapA)" and "Odocoileus HapB" with just "Odocoileus" - the HapA and HapB didn't mean anything (genetic differentiation from Levi Lab where they thought this might be a species distinguisher - it was not)
sort(unique(usable$Scientific_Name)) #first look to make sure there aren't random ones
usable$Scientific_Name <- as.factor(usable$Scientific_Name) #create factor
usable$Scientific_Name <- recode_factor(usable$Scientific_Name, 'Odocoileus (HapA)' = "Odocoileus", 'Odocoileus HapB' = "Odocoileus") # recode all to Odocoileus
View(usable)

# now replace cases where "Scientific Name" = "Odocoileus in the 'usable' database with the specific deer species designation from the 'scat.deer' database
usable.d <- as.data.frame(usable) # convert list to dataframe
class(usable.d)
usable.m <- merge(usable.d, scat.deer[, c("SampleID", "Species.ID")], by="SampleID", all.x = TRUE)
View(usable.m)
usable.m$Scientific_Name <- as.character(usable.m$Scientific_Name)
usable.m$Species.ID <- as.character(usable.m$Species.ID)
str(usable.m)

for (i in 1:length(usable.m$Scientific_Name)) {
  if (usable.m[i,]$Scientific_Name == "Odocoileus" & !is.na(usable.m[i,]$Species.ID)) {
    (usable.m[i,]$Scientific_Name <- usable.m[i,]$Species.ID)
  } else {
    (usable.m[i,]$Scientific_Name <- usable.m[i,]$Scientific_Name)
  }
}

usable.m <- select(usable.m, -Species.ID)

View(usable.m) #now we have Ellie's scat sample results combined with the metabarcoding results!

# Write a database with only the 404 useable samples (no prey only, no depositor only), already filtered at 0.05% relative read counts per rep by Taal/Jenn, with both Canis PCR and Odocoileus PCR results included
# write.csv(usable.m, file = './Data/Metabarcoding_UsableDataset_Canis_Odocoileus.csv', row.names = FALSE)

######################################################################

#to create a different format of data
diets <- usable.m

diets <- diets %>%
  group_by(SampleID) %>%
  # Create string listing all items in given Id, separated by comma
  summarise(Items = str_c(Scientific_Name, collapse = ', '))

head(diets, 12) #NB018 is the first one with Odocoileus down to spp.

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

##########
### Check for the few cases where and Id changed study areas based on error - might be none but check the old metabarcoding log that the undergrad lab interns worked on for any anomalies
##########

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
# View(diets)
table(diets$Depositor_DNA)
table(diets$Depositor_Field)
sum(table(diets$Depositor_DNA)) # gut check - should be 404
sum(table(diets$Depositor_Field)) # gut check - should be 404

# now remove values from the Diet_Item_X columns if they match the Depostitor_DNA columns
for (i in 1:length(diets$SampleID)) {
  temp <- as.vector(diets[i])
  print_val <- NULL
  if (temp$Depositor_DNA %in% temp[5:10]) {
    col_val <- match(temp[5:10], temp$Depositor_DNA)
    col_val <- as.integer(as.numeric(which(col_val==1)) + 4) # add 4 to account for the four columns we ignore (the first four columns so SampleID...DepositorID)
    print_val <- as.character(diets[i, ..col_val])
    diets[i, col_val] <- NA
  } else {
    print(cat(diets$SampleID[i], diets$Depositor_Field[i], diets$Depositor_DNA[i], print_val, "skip", " "), sep=" ")
  }
}

# Only two samples threw an error - NW077 and OW035 which are the two scats that came back as Canis  spp. only (PCR couldn't differentiate  between wolf and coyote)

# look at the resulting data
head(diets)
View(diets)

# test to make sure that no depositors are left in the "diet" columns (note that bear shows up a few times as a diet item in wolf scats)
table(diets$Diet_Item_1)
table(diets$Diet_Item_2)
table(diets$Diet_Item_3)
table(diets$Diet_Item_4)
table(diets$Diet_Item_5)
table(diets$Diet_Item_6) # no prey items in Diet_Item_6, so remove
diets <- diets[,1:9] # removing the Diet_Items_6 column that is empty

diets <- as.data.frame(diets)
str(diets)

# look at some summary data tables
# https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
with(diets, table(StudyArea, Depositor_DNA)) 

# shift missing diet values to the left
# https://stackoverflow.com/questions/49079789/using-r-to-shift-values-to-the-left-of-data-frame
diets[] <-  t(apply(diets, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
head(diets)
length(diets$SampleID) # should be 404
# View(diets)

# write.csv(diets, "./Data/diets_wide.csv", row.names = FALSE)

# now create it in long format
colnames(diets)
diets_long <- gather(diets, diet_item, prey_item, Diet_Item_1:Diet_Item_5, factor_key=TRUE) # long format
diets_long <- diets_long %>% drop_na(prey_item)
with(diets_long, diets_long[order(SampleID, prey_item),]) # sort by SampleID, prey_item
head(diets_long, 30)

# look at some summary tables of field accuracy of Depositor_Field as compared to Depositor_DNA
# note that we filter to "Diet_Item_1" to get one entry per SampleID for the purposes of making the table only
diets_long_unique <- diets_long[which(diets_long$diet_item == "Diet_Item_1"),]
length(unique(diets_long_unique$SampleID)) # should be 404
with(diets_long_unique, table(Depositor_DNA, Depositor_Field))

# we have some non-target scats deposited by coyote (Canis latrans), unknown canid (Canis spp.), and bobcat (Lynx rufus) - remove these entries
# note that we are doing this from the 'diets_long' file since we want to apply this to the full dataset 
diets_long_clean <- diets_long[!(diets_long$Depositor_DNA=="Canis latrans" | diets_long$Depositor_DNA=="Canis spp." | diets_long$Depositor_DNA=="Lynx rufus"),]
View(diets_long_clean) 
# should have removed 19 scats in total for being from non-target spp.
404-19 # new data should have 385 unique scats

# again, look at some summary tables of field accuracy of Depositor_Field as compared to Depositor_DNA
# note that we filter to "Diet_Item_1" to get one entry per SampleID for the purposes of making the table only
diets_long_unique <- diets_long_clean[which(diets_long_clean$diet_item == "Diet_Item_1"),]
length(unique(diets_long_unique$SampleID)) # should be 385
with(diets_long_unique, table(Depositor_DNA, Depositor_Field))
# so we misidentified 43 scats total out of 404 - 19 (non-target depositor, ie, coyote, bobcat, canis spp.) + 23 (one of three target depositors but misidentified)
(404-43)/404 # 89.4% accuracy overall
(385-19)/385 # 95.1% accuracy within genetically confirmed target spp.

# write the long format .csv file
# write.csv(diets_long_clean, "./Data/diets_long.csv", row.names = FALSE)



