#### SETTINGS ####
####################################################################################################!

setwd("dataset_S2/results_bin")
snpsel::pkg.lib(data.table, dplyr, magrittr, snpsel)
dir.create("pergene", F)
dir.create("persnp", F)

#### FUNCTIONS ####
####################################################################################################!

dt_transform <- function(dt, score_col, add_entrezID = TRUE, ref = NULL) {
  # Rename columns
  if("POSITION" %in% names(dt)) dt <- rename(dt, POS = "POSITION")
  if("LOGPVALUE" %in% names(dt)) dt <- rename(dt, LOG_PVALUE = "LOGPVALUE")
  if("LOGPVALUE_RIGHT" %in% names(dt)) dt <- rename(dt, LOG_PVALUE = "LOGPVALUE_RIGHT")
  # Remove columns
  padj_cols <- grep("PADJ", names(dt), value = T)
  if(length(padj_cols) > 0) dt <- select(dt, -all_of(padj_cols))
  # Reorder columns
  p_cols <- c("PVALUE", "LOG_PVALUE")
  dt <- dt %>% select(all_of(c(setdiff(names(.), p_cols), p_cols)))
  # Remove NAs 
  dt <- dt[!is.na(get(score_col))]
  # Force score as absolute values
  if(score_col == "IHS") dt[, eval(score_col) := abs(get(score_col))]
  # Annotate Entrez ID
  if(add_entrezID) {
    dt[, GENE_ID := annotate_entrezID(GENE, ref)]
    setcolorder(dt, c(names(dt)[1:grep("^GENE$", names(dt))], "GENE_ID"))
  }
  # Set NA as empty gene symbols
  dt[is.na(GENE_ID), GENE_ID := ""]
  # Return data.table
  return(dt)
}

# Function for splitting gene
split_by_gene <- function(dt, score_col) {
  
  gene_info_cols <- c("GENE", "GENE_ID")
  
  dt <- separate_genes(dt, gene_info_cols)
  dup_genes <- dt[, unique(CHR), by = GENE][duplicated(GENE), unique(GENE)]
  
  dt %>%
    filter(!GENE %in% dup_genes) %>%
    group_by_at(gene_info_cols) %>%
    summarise (
      CHR = unique(CHR),
      POS = round(mean(POS)),
      SNP = paste0(CHR, ":", POS),
      SCORE = mean(get(score_col)),
      .groups = 'drop'
    ) %>%
    arrange(CHR, POS, GENE) %>%
    select("CHR", "POS", all_of(gene_info_cols), "SCORE") %>%
    mutate(PVALUE = sapply(SCORE, function(value) {1 - sum(value >= .$SCORE) / (nrow(.)+1)})) %>%
    mutate(LOG_PVALUE = -log10(PVALUE)) %>%
    as.data.table()
  
}

# Annotate ENTREZ ID
annotate_entrezID <- function(genes, ref) {
  dt <- data.table(SYMBOL=genes)
  if(any(grepl(",", genes))) {
    isolated_genes <- !grepl(",", genes)
    dt[isolated_genes, ENTREZID := merge(dt[isolated_genes], ref, by = "SYMBOL", sort = F, all.x = T)[["ENTREZID"]]]
    dt[, ENTREZID := as.character(ENTREZID)]
    dt[!isolated_genes, ENTREZID := sapply(SYMBOL, function(x) {
      dt_tmp <- data.table(SYMBOL=unlist(strsplit(x, split=",")))
      paste(merge(dt_tmp, ref, by = "SYMBOL", sort = F)[['ENTREZID']], collapse = ",")
    })]
    entrezID <- dt[["ENTREZID"]]
  } else {
    entrezID <- merge(dt, ref, by = "SYMBOL", sort = F, all.x = T)[["ENTREZID"]]
  }
  return(entrezID)
}

#### LOAD DATA ####
####################################################################################################!

# Load ref key-genes
ref <- fread("../data/annot/gene_keys.annot")
no_geneid <- ref[is.na(ENTREZID), SYMBOL]

# Load results per snp
pbs <- readRDS("../pbs/results/RDS/Amazonia_Mesoamerica_EastAsia.Rds")$PBS
pbs_windowed <- readRDS("../pbs/results/RDS/Amazonia_Mesoamerica_EastAsia.Rds")$PBS_MEAN
pbs_windowed_contrast <- readRDS("../pbs/results/RDS/Mesoamerica_Amazonia_EastAsia.Rds")$PBS_MEAN

xpehh <- fread("../rehh/results/xpehh_amz_mes.csv")
xpehh_contrast <- fread("../rehh/results/xpehh_mes_amz.csv")

ihs <- fread("../rehh/results/ihs_amz.csv")

#### MAIN ####
####################################################################################################!

# LD-based Analysis
ihs_persnp <- dt_transform(dt = ihs, score_col = "IHS", add_entrezID = T, ref = ref)
xpehh_persnp <- dt_transform(dt = xpehh, score_col = "XPEHH", add_entrezID = T, ref = ref)
xpehh_persnp_contrast <- dt_transform(dt = xpehh_contrast, score_col = "XPEHH", add_entrezID = T, ref = ref)

xpehh_persnp[XPEHH >= quantile(XPEHH, .999), sort(unique(GENE))]
xpehh_persnp_contrast[XPEHH >= quantile(XPEHH, .999), sort(unique(GENE))]

# Fst-based Analysis
pbs_persnp <- dt_transform(dt = pbs, score_col = "PBS", add_entrezID = T, ref = ref)
pbs_windowed_persnp <- dt_transform(dt = pbs_windowed, score_col = "PBS", add_entrezID = T, ref = ref)
pbs_windowed_persnp_contrast <- dt_transform(dt = pbs_windowed_contrast, score_col = "PBS", add_entrezID = T, ref = ref)

## Get Mean Values per Gene
ihs_pergene <- split_by_gene(ihs_persnp, "IHS")
xpehh_pergene <- split_by_gene(xpehh_persnp, "XPEHH")
xpehh_pergene_contrast <- split_by_gene(xpehh_persnp_contrast, "XPEHH")

pbs_pergene <- split_by_gene(pbs_persnp, "PBS")
pbs_windowed_pergene <- split_by_gene(pbs_windowed_persnp, "PBS")
pbs_windowed_pergene_contrast <- split_by_gene(pbs_windowed_persnp_contrast, "PBS")

# Save data per snp
feather::write_feather(pbs_persnp, "persnp/pbs.feather")
feather::write_feather(pbs_windowed_persnp, "persnp/pbs_windowed.feather")
feather::write_feather(pbs_windowed_persnp_contrast, "persnp/pbs_windowed_contrast.feather")

feather::write_feather(ihs_persnp, "persnp/ihs.feather")
feather::write_feather(xpehh_persnp, "persnp/xpehh.feather")
feather::write_feather(xpehh_persnp_contrast, "persnp/xpehh_contrast.feather")

# Save data per gene
feather::write_feather(pbs_pergene, "pergene/pbs.feather")
feather::write_feather(pbs_windowed_pergene, "pergene/pbs_windowed.feather")
feather::write_feather(pbs_windowed_pergene_contrast, "pergene/pbs_windowed_contrast.feather")

feather::write_feather(ihs_pergene, "pergene/ihs.feather")
feather::write_feather(xpehh_pergene, "pergene/xpehh.feather")
feather::write_feather(xpehh_pergene_contrast, "pergene/xpehh_contrast.feather")

