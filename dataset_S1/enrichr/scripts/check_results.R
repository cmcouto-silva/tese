#### SETTINGS ####
####################################################################################################!

# WD & libraries
setwd("~/cmcouto.silva@usp.br/lab_files/tese/dataset_S1/enrichr/results")
snpsel::pkg.lib(WebGestaltR, feather, data.table, dplyr, magrittr, glue, snpsel)

## -- Arguments -- ##

modes <- c("persnp", "pergene")
methods <- c("pbs_windowed", "xpehh", "ihs")
cutoff_vals <- c("0.999", "0.995", "0.99")#, "0.95")

#### MAIN ####
####################################################################################################!

# Load results
results <- sapply(modes, function(mode) {
  sapply(methods, function(method) {
    sapply(cutoff_vals, function(cutoff) {
      fread(glue("{method}_{mode}_pvalue{cutoff}.txt"))
    }, simplify = F) %>% set_names(cutoff_vals)
  }, simplify = F)
}, simplify = F)

# Define common function to get statistical significative results
get_stat_results <- function(dt, database, p_col = "Adjusted P-value", p_cutoff = 0.05) {
  dt[Gene_set==database & get(p_col) <= p_cutoff, .(Term, `Adjusted P-value`)][order(`Adjusted P-value`)]
}

## -- GWAS_Catalog_2019 -- ##

# Check results per SNP
lapply(results$persnp$pbs_windowed, get_stat_results, "GWAS_Catalog_2019")
lapply(results$persnp$xpehh, get_stat_results, "GWAS_Catalog_2019")
lapply(results$persnp$ihs, get_stat_results, "GWAS_Catalog_2019")

# Check results per gene
lapply(results$pergene$pbs_windowed, get_stat_results, "GWAS_Catalog_2019")
lapply(results$pergene$xpehh, get_stat_results, "GWAS_Catalog_2019")
lapply(results$pergene$ihs, get_stat_results, "GWAS_Catalog_2019")

## -- KEGG_2019_Human -- ##

# Check results per SNP
lapply(results$persnp$pbs_windowed, get_stat_results, "KEGG_2019_Human")
lapply(results$persnp$xpehh, get_stat_results, "KEGG_2019_Human")
lapply(results$persnp$ihs, get_stat_results, "KEGG_2019_Human")

# Check results per gene
lapply(results$pergene$pbs_windowed, get_stat_results, "KEGG_2019_Human")
lapply(results$pergene$xpehh, get_stat_results, "KEGG_2019_Human")
lapply(results$pergene$ihs, get_stat_results, "KEGG_2019_Human")
