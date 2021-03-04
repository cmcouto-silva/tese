#### Settings ####
####################################################################################################!

# WD & libraries
setwd("dataset_S1/fumagwas")
snpsel::pkg.lib(data.table, dplyr, glue, snpsel)

# Arguments
mode <- "pergene"     # persnp, pergene
method <- "ihs"   # pbs_windowed, xpehh, ihs

#### Pipeline ####
####################################################################################################!

# Create directory for results
outdir <- glue("{mode}/{method}")
dir.create(outdir, F, T)

# Load data
dt <- setDT(feather::read_feather(glue("../results_bin/{mode}/{method}.feather")))
background_genes <- fread("../data/annot/gene_keys.annot")[ENTREZID!=""][, SYMBOL]

# Function to get top genes
get_top_genes <- function(dt, cutoff, mode = "persnp", pbsw_correction = TRUE) {
  
  for(method in c("IHS", "XPEHH", "PBS")) {
    if(method %in% names(dt)) setnames(dt, method, "SCORE")
  }
  
  if(mode=="persnp") {
    genes <- dt %>%
      filter(SCORE >= quantile(SCORE, cutoff)) %>%
      {if(pbsw_correction & cutoff > 0 & "W_ID" %in% names(.)) gdstat::exclude_consecutive_windows(., "SCORE") else . } %>%
      {if(pbsw_correction & cutoff > 0 & "W_ID" %in% names(.)) select(., -all_of(c("W_ID", "W_GROUP"))) else .} %>%
      separate_genes() %>%
      arrange(-SCORE) %>%
      distinct(GENE) %>%
      pull
  } else {
    dt <- separate_genes(dt)
    dup_genes <- dt[, unique(CHR), by = GENE][duplicated(GENE), unique(GENE)]
    genes <- dt[!GENE %in% dup_genes, .(SCORE=mean(SCORE)), GENE][SCORE >= quantile(SCORE, cutoff)][, sort(unique(GENE))]
  }
  return(genes)
}

# top 0.1%
target_genes <- get_top_genes(dt, cutoff = .999, mode = "persnp")

clipr::write_clip(target_genes)
clipr::write_clip(background_genes)

# top 0.5%
target_genes <- get_top_genes(dt, cutoff = .995, mode = "persnp")

clipr::write_clip(target_genes)
# clipr::write_clip(background_genes)

# top 1%
target_genes <- get_top_genes(dt, cutoff = .99, mode = "persnp")

clipr::write_clip(target_genes)
# clipr::write_clip(background_genes)

# # top 5%
# target_genes <- dt[get_top_genes(get(score_col), .95), isolate_genes(GENE)]
# background_genes <- dt[, isolate_genes(GENE)]

# clipr::write_clip(target_genes)
# clipr::write_clip(background_genes)

# Make a table
cutoff_vals <- c("top01", "top05", "top1")
fuma_lst <- list()

for (cutoff in cutoff_vals) {
  fuma_file <- glue("{mode}/{method}/{cutoff}.zip")
  fuma_lst[[rm.extension(fuma_file)]] <- fread(unzip(fuma_file, "GS.txt"))[
    Category %in% c("GWAScatalog","Chemical_and_Genetic_pertubation")
  ][adjP <= 0.05]
  unlink("GS.txt")
}

# Name list
names(fuma_lst) <- basename(names(fuma_lst))

# Check Results
lapply(fuma_lst, function(dt) dt[, .(Category, GeneSet, p, adjP)]) 

# Write Results
WriteXLS::WriteXLS(fuma_lst, glue("{mode}/{method}/fumagwas_tbl.xls"),
                   row.names = F, AutoFilter = T, BoldHeaderRow = T, AdjWidth = T)
clean()
