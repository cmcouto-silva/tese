#### Settings ####
####################################################################################################!

# WD & libraries
setwd("~/cmcouto.silva@usp.br/lab_files/tese/dataset_S2/webgestalt/")
snpsel::pkg.lib(WebGestaltR, feather, data.table, dplyr, magrittr, glue, snpsel)

## -- Arguments -- ##

modes <- c("persnp", "pergene")
methods <- c("pbs_windowed", "xpehh", "ihs")
cutoff_vals <- c(.999, .995, .99, .95)

gene_sets <- c (
  "geneontology_Biological_Process_noRedundant",
  "pathway_KEGG", "pathway_Reactome",
  "phenotype_Human_Phenotype_Ontology"
)  

#### MAIN ####
####################################################################################################!

## -- ORA -- ##
ora <- sapply(modes, function(mode) {
  sapply(methods, function(method) {
    sapply(cutoff_vals, function(cutoff) {
      readRDS(glue("ORA/{mode}/{method}/pvalue_{cutoff}/pvalue_{cutoff}.RDS")) %>%
        sapply(as.data.table, simplify = F)
    }, simplify = F) %>% set_names(cutoff_vals)
  }, simplify = F)
}, simplify = F)

# Function to get statistical significative results
get_stat_results <- function(mode, method) {
  lapply(cutoff_vals, function(cutoff) {
    sapply(gene_sets[2:3], function(gs) {
      ora[[mode]][[method]][[cutoff]][[gs]][FDR <= 0.05]
      
    }, simplify = F)
  }) %>% set_names(cutoff_vals)
}

get_stat_results(mode = "persnp", method = "pbs_windowed")
get_stat_results(mode = "persnp", method = "xpehh")
get_stat_results(mode = "persnp", method = "ihs")

get_stat_results(mode = "pergene", method = "pbs_windowed")
get_stat_results(mode = "pergene", method = "xpehh")
get_stat_results(mode = "pergene", method = "ihs")



# Load GSEA results
gsea <- readRDS("GSEA/gsea.RDS")

# PBS_WINDOWED
sapply(names(gsea$pbs_windowed)[2:3], function(gs) {
  gsea$pbs_windowed[[gs]][FDR < 0.05][, .(geneSet, description, pValue, FDR)]# [pValue < 0.05] [FDR < 0.05]
}, simplify = F)

# XP-EHH
sapply(names(gsea$pbs_windowed)[2:3], function(gs) {
  gsea$xpehh[[gs]][FDR < 0.05][, .(geneSet, description, pValue, FDR)]# [pValue < 0.05] [FDR < 0.05]
}, simplify = F)

# iHS
sapply(names(gsea$pbs_windowed)[2:3], function(gs) {
  gsea$ihs[[gs]][FDR < 0.05][, .(geneSet, description, pValue, FDR)]# [pValue < 0.05] [FDR < 0.05]
}, simplify = F)

