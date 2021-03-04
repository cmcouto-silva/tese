#### SETTINGS ####
####################################################################################################!

# Set working directory
#setwd("dataset_S1/pbs")

# Load required libraries
#snpsel::pkg.lib(data.table, magrittr, snpsel)

# Load required files
#popRef <- fread("../data/populations/populations.csv")

#### DATASET ####
####################################################################################################!

# Create list of target-individuals
target_ids <- popRef[Ecoregion != "" | Pop7Groups == "Est_Asia", .(FID, IID)]
fwrite(target_ids, "target_ids.txt", sep = " ", col.names = F)

# Filter target-individuals
dir.create("data", F)
plink2(`--bfile` = "../data/phased/dataset_S1", `--keep` = "target_ids.txt",
       "--make-bed", `--out` = "data/pbs_dataset")
unlink("target_ids.txt")

# Update Family IDs
fam_file <- "data/pbs_dataset.fam"
fam <- fread(fam_file)

fam_merge <- merge(fam[, .(IID = V2)], popRef, by = "IID", sort = F)
fam[, V1 := fam_merge[, Ecoregion]]
fam[V1 == "", V1 := "EastAsia"]

fwrite(fam, fam_file, sep = " ", col.names = F)
