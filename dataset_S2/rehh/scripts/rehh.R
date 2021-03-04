#### SETTINGS ####
####################################################################################################!

# Setting up working directory
setwd("dataset_S2/rehh")

# Load required packages
snpsel::pkg.lib(data.table, magrittr, glue, snpsel, rehh)

## --  Setup Arguments -- ##

# VCF File
vcf_file <- "../data/phased_adj_alleles/dataset_S2.vcf.gz"
vcf_subfile <- "data/rehh_dataset"

# Target Chromosomes
chrs <- 1:22

# Gene Dist Reference (for annotation)
gene_dist_ref <- "../data/annot/gene_dist.annot"

# Target Populations
populations <- fread("../data/populations/populations.csv")
pop_subset <- populations[sel_group!="", paste0(FID, ":", IID)]
# Get list of target-populations
target_pops_names <- populations[sel_group!="", unique(sel_group)]
target_pops <- sapply(target_pops_names, function(pop) {
  populations[sel_group==pop, paste0(FID, ":", IID)]
}) %>% set_names(target_pops_names)

#### MAIN ####
####################################################################################################!

## -- Main -- ##

# Subset data to faster analysis
dir.create("data", F)
subset_vcf(vcf_file, populations[sel_group!="", paste0(FID, ":", IID)], out = vcf_subfile)

# ## -- With Parallelization -- ##
# 
# # Set n cores
# cores <- parallel::detectCores() - 1 # for making haplohh data
# threads <- parallel::detectCores() - 1 # for making scanhh data
# 
# # Register cluster
# clst <- makeCluster(threads)
# registerDoSNOW(clst)
# 
# # Make progress bar
# pb <- txtProgressBar(max = length(chrs), style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# 
# # Load data into haplohh object
# haplohh_list <- foreach (
#   chr=chrs, .packages = c("rehh", "glue"), .export = c("vcf_subfile", "chrs"), #.verbose = T,
#   .options.snow = list(progress = progress)) %dopar% {
#     haplohh <- data2haplohh(glue("{vcf_subfile}.vcf.gz"), chr.name = chr, polarize_vcf = F)
#     return(haplohh)
#     }
# 
# names(haplohh_list) <- paste0("chr", chrs)
# 
# # Close progress bar and close cluster
# close(pb)
# stopCluster(clst) 

## -- Without Parallelization -- ##

# Load data into haplohh object
haplohh_list <- lapply(chrs, function(chr) {
  cat(glue("Load data for chromosome {chr}\n\n"))
  cat("-------------------------------\n")
  haplohh <- data2haplohh(glue("{vcf_subfile}.vcf.gz"), chr.name = chr, polarize_vcf = F)
  cat("\n")
  return(haplohh)
}) %>% set_names(paste0("chr", chrs))

# Split chr list by population
haplohh_list <- sapply(target_pops_names, function(pop) {
  lapply(chrs, function(i) {
    haplohh <- haplohh_list[[paste0("chr", i)]]
    haplohh@haplo <- haplohh@haplo[rownames(haplohh@haplo) %in% paste0(target_pops[[pop]], c("_1", "_2")), ]
    return(haplohh)
  }) %>% set_names(paste0("chr", chrs))
}, simplify = F)

# Save haplohh object
dir.create("bin", F)
saveRDS(haplohh_list, "bin/haplohh_list.RDS")

# Make scanhh per population
scanhh_list <- sapply(target_pops_names, function(pop) {
  cat(glue("Calculating EHH for population {pop}...\n\n"))
  
  scanhh_chrlist <- lapply(haplohh_list[[pop]], scan_hh)
  # scanhh_chrlist <- lapply(haplohh_list[[pop]], scan_hh, threads = threads)
  scanhh <- do.call(rbind, unname(scanhh_chrlist))
  return(scanhh)
}, simplify = F)

# Save haplohh object
dir.create("bin", F)
saveRDS(scanhh_list, "bin/scanhh_list.RDS")

## -- Apply Main Statistics -- ##

# Set function for annotating and tidying data
transform_data <- function(df, separate_genes = F) {
  dt <- as.data.table(df, keep.rownames = T) %>%
    setnames("rn", "SNP")
  
  dt[, c("GENE","FUNCTION") := annotate_fromtbl(SNP, ref = gene_dist_ref)]
  
  log_col <- grep("LOGPVALUE", names(dt), value = T)
  dt[, PVALUE := 10^-get(log_col)]
  setcolorder(dt, c("CHR", "POSITION", "SNP", "GENE", "FUNCTION"))
  
  dt[, PADJ := p.adjust(PVALUE, method = "fdr")]
  dt[, LOG_PADJ := -log10(PADJ)]
  
  if(separate_genes) dt <- separate_genes(dt)
  dt <- dt[, .SD, .SDcols = c(1:(length(dt)-2), length(dt), length(dt)-1)][]
  
  return(dt)
}

## -- Applying iHS Statistics -- ##
ihs_amz <- ihh2ihs(scanhh_list[["Amazonia"]], include_freq = T, right = T)
ihs_amz <- transform_data(ihs_amz$ihs)

ihs_and <- ihh2ihs(scanhh_list[["Andes"]], include_freq = T, right = T)
ihs_and <- transform_data(ihs_and$ihs)

## -- Applying XP-EHH Statistics -- ##
xpehh_amz_and <- ies2xpehh (
  scan_pop1 = scanhh_list$Amazonia, scan_pop2 = scanhh_list$Andes,
  include_freq = T, p.side = "right") %>%
  transform_data()

xpehh_amz_mes <- ies2xpehh (
  scan_pop1 = scanhh_list$Amazonia, scan_pop2 = scanhh_list$Mesoamerica,
  include_freq = T, p.side = "right") %>%
  transform_data()

xpehh_amz_eas <- ies2xpehh (
  scan_pop1 = scanhh_list$Amazonia, scan_pop2 = scanhh_list$EastAsia,
  include_freq = T, p.side = "right") %>%
  transform_data()

xpehh_and_amz <- ies2xpehh (
  scan_pop1 = scanhh_list$Andes, scan_pop2 = scanhh_list$Amazonia,
  include_freq = T, p.side = "right") %>%
  transform_data()

xpehh_mes_amz <- ies2xpehh (
  scan_pop1 = scanhh_list$Mesoamerica, scan_pop2 = scanhh_list$Amazonia,
  include_freq = T, p.side = "right") %>%
  transform_data()

# Saving results
dir.create("results", F)
fwrite(ihs_amz, "results/ihs_amz.csv")
fwrite(ihs_and, "results/ihs_and.csv")
fwrite(xpehh_amz_and, "results/xpehh_amz_and.csv")
fwrite(xpehh_amz_mes, "results/xpehh_amz_mes.csv")
fwrite(xpehh_amz_eas, "results/xpehh_amz_eas.csv")
fwrite(xpehh_and_amz, "results/xpehh_and_amz.csv")
fwrite(xpehh_mes_amz, "results/xpehh_mes_amz.csv")

# Delete VCF data
unlink("data", recursive = T)
