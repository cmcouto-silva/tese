#### SETTINGS ####
####################################################################################################!

# Set working directory
setwd("dataset_S2/pbs")

# Load required libraries
snpsel::pkg.lib (
  data.table, magrittr, snpsel, glue,
  ggplot2, ggrepel, cowplot, gdstat
)

# Load population annotations
popRef <- fread("../data/populations/populations.csv")

# Load pbs pipeline
source("scripts/pbs_pipeline.R")

# Make input
source("scripts/pbs_setup.R")
rm(fam, fam_merge, fam_file, target_ids)

#### MAIN ANALYSIS ####
####################################################################################################!

# Set target groups
eas <- popRef[sel_group == "EastAsia", IID]
mes <- popRef[sel_group == "Mesoamerica", IID]
and <- popRef[sel_group == "Andes", IID]
amz <- popRef[sel_group == "Amazonia", IID]
amz_filt <- popRef[sel_group == "Amazonia" & !grepl("HighSelva", FID), IID]

# AMZ x AND x MES
pbs_pipeline(plink_file = "data/pbs_dataset", focal = amz, close = and, outgroup = mes,
             focal_name = "Amazonia", close_name = "Andes", outgroup_name = "Mesoamerica",
             monomorphic_removal = T, maf = 0.05, nchrobs = 2, filter = "ID",
             method = "rwc")

# AMZ x MES x EAS
pbs_pipeline(plink_file = "data/pbs_dataset", focal = amz, close = mes, outgroup = eas,
             focal_name = "Amazonia", close_name = "Mesoamerica", outgroup_name = "EastAsia",
             monomorphic_removal = T, maf = 0.05, nchrobs = 2, filter = "ID",
             method = "rwc")

# MES x AMZ x EAS
pbs_pipeline(plink_file = "data/pbs_dataset", focal = mes, close = amz, outgroup = eas,
             focal_name = "Mesoamerica", close_name = "Amazonia", outgroup_name = "EastAsia",
             monomorphic_removal = T, maf = 0.05, nchrobs = 2, filter = "ID",
             method = "rwc")

# AMZ (Without Yanesha HighSelva) x AND x EAS
pbs_pipeline(plink_file = "data/pbs_dataset", focal = amz_filt, close = and, outgroup = eas,
             focal_name = "AMZFilt", close_name = "Andes", outgroup_name = "EastAsia",
             monomorphic_removal = T, maf = 0.05, nchrobs = 2, filter = "ID",
             method = "rwc")

# AMZ (Without Yanesha HighSelva) x MES x EAS
pbs_pipeline(plink_file = "data/pbs_dataset", focal = amz_filt, close = mes, outgroup = eas,
             focal_name = "AMZFilt", close_name = "Mesoamerica", outgroup_name = "EastAsia",
             monomorphic_removal = T, maf = 0.05, nchrobs = 2, filter = "ID",
             method
             = "rwc")

# AND x AMZ x MES
pbs_pipeline(plink_file = "data/pbs_dataset", focal = and, close = amz, outgroup = mes,
             focal_name = "Andes", close_name = "Amazonia", outgroup_name = "Mesoamerica",
             monomorphic_removal = T, maf = 0.05, nchrobs = 2, filter = "ID",
             method = "rwc")

# Delete Plink data
unlink("data", recursive = T)
