# CONTINUED DATA CLEANING FOR SCAT METABARCODING ANALYSIS (PART 2) ---------------------
# data formatting for bootstrapping, %FO plots, NMDS plots, Shannon's Index, Pianka's Index, comparison plots (dietary diversity vs. dietary overlap)

# Clear workspace ------------------------------------------
rm(list = ls())

getwd()

# Species Breakdowns - Scat Data ----------------------------



#    # only need to do the steps below once because I had to do some manual editing to match ClusterID in the scatlog to the cluster database so the data would extract properly

#   # note that I manually added values to the "prey_simple_deerspp" and "prey_simple_unkdeer" columns in Excel to create the diets_long_PERMANOVA.csv from the diets_long.csv - code this in for future (9/7/22)
# 
#   # read in the metabarcoding results from the Scat Log so that we can add in season, cluster/non-cluster, carcass found/not found, etc.
#   
# scatlog <- read.csv("Data/ScatLog_MetabarcodingSamples.csv")
# 
#   # create a new column based on the "Found_Along" column named "Cluster" - Yes = found "AtCluster", No = otherwise
# table(scatlog$Found_Along)
# scatlog <- scatlog %>%
#   mutate(
#     Cluster = recode(
#       Found_Along,
#       'AtCluster' = 'Yes',
#       'AtCluster, Off-trail' = 'Yes',
#       'AtCluster, Road' = 'Yes',
#       'AtCluster, Trail' = 'Yes',
#       'Follow' = 'No',
#       'Off-trail' = 'No',
#       'Other' = 'No',
#       'Road' = 'No',
#       'Trail' = 'No'
#     )
#   )
# with(scatlog, table(Metabarcoding_Species, Cluster))
# table(scatlog$Cluster)
# 
#   # reclassify "Season" into just 'Summer' (Summer + Fall) and 'Winter' (Winter + Spring)
#   # note that seasons are as follows: Summer = 6/1 - 9/14; Fall = 9/15 - 12/14; Winter = 12/15 - 3/14; Spring = 3/15 - 5/31
# table(scatlog$Season)
# scatlog <- scatlog %>%
#   mutate(
#     Season = recode(
#       Season,
#       'Winter' = 'Winter',
#       'Spring' = 'Winter',
#       'Summer' = 'Summer',
#       'Fall' = 'Summer'
#     )
#   )
# with(scatlog, table(Season, Metabarcoding_Species))
# table(scatlog$Season) # check to see if it worked
# table(scatlog$Date_Sent_OSU) # check that all submitted samples were counted by comparing to previous, should sum to the same value

  # read in the cluster database files ------------------------
clusters_reg <- read.csv("Data/WolfCougarClusters_Database.csv")
clusters_winter <- read.csv("Data/Winter2020_Database.csv")
  # check that they read in correctly and have the same number of columns
head(clusters_reg)
head(clusters_winter)
dim(clusters_reg)
dim(clusters_winter)
# merge the two databases into one
clusters <- dplyr::bind_rows(clusters_reg, clusters_winter)

#   # add and manually edit ClusterID column then comment out ---------------------------
#   # create column for carcass found/not found at cluster in the scatlog database so we know how many scats came from clusters with prey (true feeding sites) vs from clusters without prey (presumed resting sites)
#   # only need to do this once - comment out after
#   # subset the scatlog to only scats that were sent for metabarcoding analysis
# scatlog_meta <- scatlog[!(scatlog$Date_Sent_OSU == ""),]
# head(scatlog_meta)
#   # create a new column for cluster ID in the scat log
#   # this code removed the last - and everything after (eg 52F-112-bear1 becomes 52F-112)
# scatlog_meta <-
#   scatlog_meta %>% mutate(Cluster_ID = gsub('(.*)-\\w+', '\\1', Scat_ID),
#                           .after = Scat_ID)
#   # write the .csv to a file
# write.csv(scatlog_meta, file = "scatlog_meta.csv", row.names = FALSE)
#   # manually ensure that the names in the ClusterID column match the cluster names in the cluster database so later database pairing will work
#   # note that some clusters had -A added to them, and some had -CSU added, while others I changed the -CPU to -CSU since the -CSU (camera set up) entry has info on carcass found, etc. even though scats were collected on the CPU (camera pick up) trip
#   # also need to check on two clusters missing from the cluster database: 52F-1362 and MVC207F-729 
# read in the edited scatlog_meta_ClusterID.csv file -------------------
  # and name it scatlog_meta
scatlog_meta <- read.csv("Data/scatlog_meta_ClusterID.csv")

  # assign carcass found/not found to each scat ID by cluster ID
scatlog_meta <-
  merge(scatlog_meta, clusters[, c(
    "Cluster_ID",
    "CarcassFound",
    "PreySpecies_Final",
    "Sex_Final",
    "AgeClassDeer_Final",
    "AgeEstimate_Final"
  )], by = "Cluster_ID", all.x = TRUE)

with(scatlog_meta,
     table(Metabarcoding_Species, Metabarcoding_Status))
with(scatlog_meta,
     table(CarcassFound, Cluster, Metabarcoding_Species))

  # if Cluster = No and CarcassFound = NA, reassign the latter to "No"
scatlog_meta <- within(scatlog_meta, CarcassFound[ is.na(CarcassFound) & Cluster == "No"] <- "No" ) 
  # correct Cluster and CarcassFound columns for some of the obscure clusters
  # assigned Cluster = No and CarcassFound = Yes to opportunistic clusters, ungulate morts, etc. where there was a carcass but it was not a wolf or cougar cluster
  # Note that 48F-197 had a female moose whose tooth aged her at about 13 years old. However, I coded this incident as "No Carcass Found" for the purposes of documenting wolf predation patterns since this strongly appeared to be scavenging, not predation. That said, for the purposes of the diet analysis, where we want to know if scats are collected near a carcass at a cluster or not, this was coded as "Carcass Found" since that is the more applicable interpretation in this context
  # # 48F-197: CarcassFound = Yes, PreySpecies_Final = moose, Sex_Final = Female
scatlog_meta[scatlog_meta$Cluster_ID == "48F-197", "CarcassFound"] <- "Yes"
scatlog_meta[scatlog_meta$Cluster_ID == "48F-197", "PreySpecies_Final"] <- "moose"
scatlog_meta[scatlog_meta$Cluster_ID == "48F-197", "Sex_Final"] <- "Female"
  # # 3848WTD19: True case of an "NA" for cluster/carcass because this was a collared deer mortality
  # ## CHECK ON 52F-1362 - not in cluster database??
scatlog_meta <- within(scatlog_meta, CarcassFound[Cluster_ID == "MVC209M-314"] <- "Yes" ) # was a dead cougar so does not have a prey species listed
scatlog_meta <- within(scatlog_meta, Cluster[Cluster_ID == "MVC209M-243"] <- "Yes" ) 
scatlog_meta <- within(scatlog_meta, Cluster[Cluster_ID == "NEC107F-669"] <- "Yes" ) 
scatlog_meta <- within(scatlog_meta, Cluster[Cluster_ID == "NEC107F-1827"] <- "Yes" )
scatlog_meta <- within(scatlog_meta, Cluster[Cluster_ID == "NEC104F-1515"] <- "Yes" )
scatlog_meta <- within(scatlog_meta, Cluster[Cluster_ID == "NEC103F-427"] <- "Yes" )
scatlog_meta <- within(scatlog_meta, Cluster[Cluster_ID == "NEC103F-326"] <- "Yes" )
scatlog_meta <- within(scatlog_meta, Cluster[Cluster_ID == "NEC101M-1697"] <- "Yes" )
scatlog_meta <- within(scatlog_meta, Cluster[Cluster_ID == "MVC207F-729"] <- "Yes" )
scatlog_meta <- within(scatlog_meta, CarcassFound[ CarcassFound == 'Mortality' & Cluster_ID == "MVC228F-2274-A-CSU"] <- "Yes" )  
scatlog_meta <- within(scatlog_meta, CarcassFound[Cluster_ID == "3848WTD19"] <- "Yes" ) 
scatlog_meta <- within(scatlog_meta, CarcassFound[Cluster_ID == "NEOP-03JULY2019"] <- "Yes" ) 
scatlog_meta <- within(scatlog_meta, CarcassFound[Cluster_ID == "NEOP-16March2019"] <- "Yes" )

# write the edited scatlog_meta_ClusterID_2.csv file to the Data folder ----
# upload to the ScatDietNMDS&PERMANOVA.R file starting with this version as all previous steps are still data prep, and some were done manually
# write.csv(scatlog_meta, "Data/scatlog_meta_ClusterID_2.csv")
# 
# Made some more manual edits so the final file is "scatlog_meta_ClusterID_Final.csv"
