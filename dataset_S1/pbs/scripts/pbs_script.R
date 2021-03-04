#### SETTINGS ####
####################################################################################################!

# Set working directory
setwd("dataset_S1/pbs")

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

# AMZ x MES x EAS
pbs_pipeline(plink_file = "data/pbs_dataset", focal = "Amazonia", close = "Mesoamerica", outgroup = "EastAsia",
             focal_name = "Amazonia", close_name = "Mesoamerica", outgroup_name = "EastAsia",
             monomorphic_removal = F, maf = 0.05, nchrobs = 2, filter = "FID",
             method = "rwc")

# MES x AMZ x EAS
pbs_pipeline(plink_file = "data/pbs_dataset", focal = "Mesoamerica", close = "Amazonia", outgroup = "EastAsia",
             focal_name = "Mesoamerica", close_name = "Amazonia", outgroup_name = "EastAsia",
             monomorphic_removal = F, maf = 0.05, nchrobs = 2, filter = "FID",
             method = "rwc")

# Delete Plink data
unlink("data", recursive = T)
