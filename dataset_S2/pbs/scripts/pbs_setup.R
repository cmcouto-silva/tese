#### SETTINGS ####
####################################################################################################!

# # Set working directory
# setwd("dataset_S2/pbs")
# 
# # Load required libraries
# snpsel::pkg.lib(data.table, magrittr, snpsel)
# 
# # Load required files
# popRef <- fread("../data/populations/populations.csv")

#### DATASET ####
####################################################################################################!

# Create folder to store dataset
dir.create("data", F)

# Create list of target-individuals
target_ids <- popRef[sel_group != "", .(FID, IID)]
fwrite(target_ids, "data/target_ids.txt", sep = " ", col.names = F)

# Filter target-individuals
plink2 (
  `--bfile` = "../data/phased/dataset_S2", `--keep` = "data/target_ids.txt",
  "--make-bed", `--out` = "data/pbs_dataset"
  )
unlink("data/target_ids.txt")

# Update Family IDs
fam_file <- "data/pbs_dataset.fam"
fam <- fread(fam_file)

fam_merge <- merge(fam[, .(IID = V2)], popRef, by = "IID", sort = F)
fam[, V1 := fam_merge[, sel_group]]

fwrite(fam, fam_file, sep = " ", col.names = F)
