#### SETUP ####
####################################################################################################!

# Libraries
snpsel::pkg.lib(data.table, magrittr, glue, snpsel)
snpsel::pkg.lib(dplyr, ggplot2, cowplot)
library(kableExtra)
library(openxlsx)

#####  Functions #####

# Source functions
source("scripts/manhattan.R")
source("scripts/efotrait2DT.R")
source("scripts/gene2efotrait.R")

# Verify p-values
get_dt <- function(df) {
  for(method in c("IHS", "XPEHH", "PBS")) {
    if(method %in% names(df)) df %<>% rename("SCORE" = all_of(method))
  }
  dt <- df %>%
    mutate(PVALUE = 1 - frank(SCORE) / (n()+1)) %>%
    mutate(LOG_PVALUE = -log10(PVALUE)) %>%
    as.data.table()
  return(dt)
}

# Get top values
get_top <- function(dt, cutoff, mode = "snp", pbsw_correction = TRUE) {
  
  dt <- copy(dt)
  
  grp <- c("CHR", "POS", "SNP")
  p_cols <- c("PVALUE", "LOG_PVALUE")
  gene_info_cols <- c("GENE", "GENE_ID", "FUNCTION")
  
  if(mode=="snp") {
    
    dt <- dt %>%
      filter(SCORE >= quantile(SCORE, cutoff)) %>%
      {if(pbsw_correction & cutoff > 0 & "W_ID" %in% names(.)) gdstat::exclude_consecutive_windows(., "SCORE") else . } %>%
      {if(pbsw_correction & cutoff > 0 & "W_ID" %in% names(.)) select(., -all_of(c("W_ID", "W_GROUP"))) else .} %>%
      select(all_of(c(grp, gene_info_cols, "SCORE", p_cols))) %>%
      separate_genes(gene_info_cols) %>%
      group_by(GENE) %>%
      slice(which.max(SCORE)) %>%
      arrange(CHR, POS, GENE) %>%
      ungroup() %>%
      group_by_at(setdiff(names(.), gene_info_cols)) %>%
      summarise_all(~ paste(.x, collapse = ",")) %>%
      select(all_of(c(grp, gene_info_cols, "SCORE", p_cols))) %>%
      as.data.table()
  } else {
    grp <- grp[1:length(grp)-1]
    gene_info_cols <- gene_info_cols[1:length(gene_info_cols)-1]
    
    dt <- separate_genes(dt, gene_info_cols)
    dup_genes <- dt[, unique(CHR), by = GENE][duplicated(GENE), unique(GENE)]
    
    dt <- dt %>%
      filter(!GENE %in% dup_genes) %>%
      group_by_at(gene_info_cols) %>%
      summarise (
        CHR = unique(CHR),
        POS = round(mean(POS)),
        SNP = paste0(CHR, ":", POS),
        SCORE = mean(SCORE),
        .groups = 'drop'
      ) %>%
      arrange(CHR, POS, GENE) %>%
      select(all_of(c(grp, gene_info_cols, "SCORE"))) %>%
      mutate(PVALUE = 1 - frank(SCORE) / (n()+1)) %>%
      mutate(LOG_PVALUE = -log10(PVALUE)) %>%
      filter(SCORE >= quantile(SCORE, cutoff)) %>%
      group_by_at(setdiff(names(.), gene_info_cols)) %>%
      summarise_all(~ paste(.x, collapse = ","), .groups = 'drop') %>%
      ungroup() %>%
      select(all_of(c(grp, gene_info_cols, "SCORE", p_cols))) %>%
      as.data.table()
  }
  return(dt)
}

get_top_intersection <- function(dt, cutoff, mode = "persnp", pbsw_correction = TRUE) {
  dt <- get_top(dt, cutoff, mode = mode, pbsw_correction = pbsw_correction)
  strsplit(dt$GENE, ",") %>% unlist %>% unique %>% sort
}

#####  Folder Structure  #####

dir.create("tables/tex", F, T)
dir.create("tables/xlsx", F, T)

dir.create("figures/png", F, T)
dir.create("figures/pdf", F, T)

#### LOAD DATA ####
####################################################################################################!

# Methods
methods <- c("pbs", "pbs_windowed", "ihs", "xpehh")

## -- DATASET S1 -- ##

# Per SNP
ds1_persnp <- sapply(methods, function(method) {
  dt <- feather::read_feather(glue::glue("dataset_S1/results_bin/persnp/{method}.feather")) %>%
    get_dt()
}, simplify = F)

# Per Gene
ds1_pergene <- sapply(methods, function(method) {
  dt <- feather::read_feather(glue::glue("dataset_S1/results_bin/pergene/{method}.feather")) %>%
    get_dt()
}, simplify = F)

# Get outliers
ds1_persnp_top <- lapply(ds1_persnp, get_top, .999, "snp")
ds1_pergene_top <- lapply(ds1_pergene, get_top, .999, "gene")

## -- DATASET S2 -- ##

# Per SNP
ds2_persnp <- sapply(methods, function(method) {
  dt <- feather::read_feather(glue::glue("dataset_S2/results_bin/persnp/{method}.feather")) %>%
    get_dt()
}, simplify = F)

# Per Gene
ds2_pergene <- sapply(methods, function(method) {
  dt <- feather::read_feather(glue::glue("dataset_S2/results_bin/pergene/{method}.feather")) %>%
    get_dt()
}, simplify = F)

# Get outliers
ds2_persnp_top <- lapply(ds2_persnp, get_top, .999, "snp")
ds2_pergene_top <- lapply(ds2_pergene, get_top, .999, "gene")

#### RUN SCRIPTS ####
####################################################################################################!

# Populations
source("scripts/pop_tables.R")

# Set common gene pattern exclusion for scan tables/figures
exclude_gene_pattern = "^$"  #"^LOC|^LINC|^$"

source("scripts/pbs.R")
source("scripts/xpehh.R")
source("scripts/ihs.R")

source("scripts/density_plot.R")
source("scripts/get_intersection.R")

source("scripts/gtex_plot.R")
source("scripts/enrichment.R")
source("scripts/pheno_tables.R")

#### ANNOTATE ####
####################################################################################################!

# # Set common function
# get_genes <- function(genes, rm_pattern = "^$") {
#   genes[!grepl(rm_pattern, genes)] %>%
#     strsplit(",") %>% unlist %>%
#     unique %>% sort
#   }
# 
# # Target methods (remove "pbs" only)
# target_methods <- c("pbs_windowed", "xpehh", "ihs")
# 
# ## -- Dataset S1 -- ##
# 
# # Per SNP
# ds1_persnp_annot <- lapply(ds1_persnp_top[target_methods], function(dt) {
#   efolist <- sapply(dt[, get_genes(GENE)], gene2efotrait, simplify = F)
#   efotable <- efolist2DT(efolist[!is.na(efolist)])
#   list(efolist=efolist, efotable=efotable)
# })
# 
# # Per Gene
# ds1_pergene_annot <- lapply(ds1_pergene_top[target_methods], function(dt) {
#   efolist <- sapply(dt$GENE, gene2efotrait, simplify = F)
#   efotable <- efolist2DT(efolist[!is.na(efolist)])
#   list(efolist=efolist, efotable=efotable)
# })
# 
# ## -- Dataset S2 -- ##
# 
# # Per SNP
# ds2_persnp_annot <- lapply(ds2_persnp_top[target_methods], function(dt) {
#   efolist <- sapply(dt[, get_genes(GENE)], gene2efotrait, simplify = F)
#   efotable <- efolist2DT(efolist[!is.na(efolist)])
#   list(efolist=efolist, efotable=efotable)
# })
# 
# # Per Gene
# ds2_pergene_annot <- lapply(ds2_pergene_top[target_methods], function(dt) {
#   efolist <- sapply(dt$GENE, gene2efotrait, simplify = F)
#   efotable <- efolist2DT(efolist[!is.na(efolist)])
#   list(efolist=efolist, efotable=efotable)
# })
# 
# ## -- Search -- ##
# 
# term <- "tryp"   # "mosqu", "tryp", "immun", "thermo", "lymphocy", "vir", "repro", "thyr"
# 
# tmp <- function(data, len=FALSE) {
#   out <- rbindlist(data) %>% .[, unlist(strsplit(GENES, ","))] %>%
#     unique %>% sort
#   if(len) {
#     return(length(out))
#   } else {
#     return(out)
#   }
# }
# 
# # Dataset S1 - per SNP
# lapply(ds1_persnp_annot, function(annot) {
#   annot$efotable[grepl(term, MAPPED_TRAIT, T)]
# })# %>% tmp(T)
# 
# # Dataset S1 - per Gene
# lapply(ds1_pergene_annot, function(annot) {
#   annot$efotable[grepl(term, MAPPED_TRAIT, T)]
# }) # %>% tmp(F)
# 
# # Dataset S2 - per SNP
# lapply(ds2_persnp_annot, function(annot) {
#   annot$efotable[grepl(term, MAPPED_TRAIT, T)]
# })# %>% tmp(F)
# 
# # Dataset S2 - per Gene
# lapply(ds2_pergene_annot, function(annot) {
#   annot$efotable[grepl(term, MAPPED_TRAIT, T)]
# }) # %>% tmp(F)
