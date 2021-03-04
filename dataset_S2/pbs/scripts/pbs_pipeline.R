#### MAIN ANALYSIS ####
####################################################################################################!

pbs_pipeline <- function(plink_file, focal, close, outgroup, focal_name, close_name, outgroup_name,
                         monomorphic_removal = F, maf = 0.05, nchrobs = 2L, filter = "FID", method = "rwc") {
  
  ## ---- Common Arguments ---- ##
  
  # Plink file (dataset)
  plink_file <- plink_file
  
  # Groups
  focal <- focal
  close <- close
  outgroup <- outgroup
  
  # Output names
  if(missing(focal_name)) focal_name <- focal
  if(missing(close_name)) close_name <- close
  if(missing(outgroup_name)) outgroup_name <- outgroup
  comparison <- paste0(focal_name, "_", close_name, "_", outgroup_name)
  
  ## ---- Create result folders ---- ##
  
  dir.create("figures", F)
  dir.create("results", F)
  dir.create("results/RDS", F)
  dir.create("results/XLS", F)
  
  ## ---- Calculate PBS ---- ##
  
  # Calculate PBS with Fst estimates from Reynolds, Weir & Cockerham (1983)
  pbs <- pbs_fst(plink = plink_file, focal = focal, close = close, outgroup = outgroup,
                 monomorphic_removal = monomorphic_removal, maf = maf, nchrobs = nchrobs,
                 filter = filter, method = method)
  
  # Adding p-values
  pbs[, PVALUE := 1 - frankv(PBS, ties = 'max') / (.N + 1)]
  pbs[, LOG_PVALUE := -log(PVALUE, 10)]
  
  # Adding gene annotation
  pbs[, c("GENE", "FUNCTION") := annotate_fromtbl(SNP, ref = "../data/annot/gene_dist.annot")]
  # pbs[GENE=="", .N]
  
  # # Annotate GWAS from top SNPs
  # pbs[LOG_PVALUE >= 2, GWAS := get_gwas_trait(SNP)]
  
  # --- Each SNP --- #
  
  # Getting peaks
  pbs_peaks <- pbs[LOG_PVALUE >= 2]
  qtl <- quantile(pbs[, LOG_PVALUE], c(0.995, 0.999, 0.9999))
  
  # Plotting
  mplot <- manhattan_plot(pbs_data = pbs, col_name = 'LOG_PVALUE') + geom_hline(yintercept = qtl, lty = 2, col = gray(.1))
  mgplot <- mplot.add.genes(data = pbs_peaks, mplot = mplot, by = "CHR", merge_col = "SNP", label = "geom_label_repel",
                            log_pvalue = qtl["99.9%"], size = 3, vjust = 1, nudge_x = 1)
  # mpplot <- mplot.add.points(DT = pbs, rolling_DT = NULL, mplot = mgplot, target_snp = target_snps, merge_col = "SNP", uniq = F)
  
  plt <- mgplot +
    ggtitle(paste(focal_name, "vs", close_name, "vs", outgroup_name)) +
    theme (
      plot.title = element_text(hjust = 0.5, face="bold")
    )
  
  # Save plot
  cowplot::save_plot(paste0("figures/", comparison, ".png"), plt, base_height = 8, base_width = 14)
  
  # --- SNP Mean --- #
  
  # Analyze PBS by sliding window mean
  pbs_mean <- sliding_window(DT = pbs, window_size = 20, step = 5, mode = "exact", col = "PBS", keep.cols = T)
  pbs_mean[, c("PVALUE", "LOG_PVALUE") := NULL]
  pbs_mean[PBS < 0, PBS := 0]
  update_fields(pbs, pbs_mean)
  
  # Calculate p-values
  pbs_mean[, PVALUE := 1 - frankv(PBS, ties = 'max') / (.N + 1)]
  pbs_mean[, LOG_PVALUE := -log(PVALUE, 10)]
  
  # Getting peaks
  pbs_mean_peaks <- pbs_mean[LOG_PVALUE >= 2]
  pbs_mean_peaks_nocons_windows <- exclude_consecutive_windows(pbs_mean_peaks) # exclude consecutive windows
  
  qtl_mean <- quantile(pbs_mean[, LOG_PVALUE], c(0.995, 0.999))
  
  # Rolling Mean
  mplot_mean <- manhattan_plot(pbs_data = pbs_mean, col_name = 'LOG_PVALUE') +
    geom_hline(yintercept = qtl_mean, lty = 2, col = gray(.1))
  
  mgplot_mean <- mplot.add.genes(data = pbs_mean, mplot = mplot_mean, log_pvalue = qtl_mean["99.9%"],
                                 size = 3, vjust = 1, nudge_x = 1)
  # mpplot_mean <- mplot.add.points(DT = pbs, rolling_DT = pbs_mean, mplot = mgplot_mean, target_snp = target_snps, merge_col = "W_ID", uniq = F)
  
  plt_mean <- mgplot_mean +
    ggtitle(paste(focal_name, "vs", close_name, "vs", outgroup_name)) +
    theme (
      plot.title = element_text(hjust = 0.5, face="bold")
    )
  
  # Save plot
  cowplot::save_plot(paste0("figures/", comparison, "_mean.png"), plt_mean, base_height = 8, base_width = 14)
  
  # Save results
  saveRDS(object = list(PBS = pbs, PBS_MEAN = pbs_mean), file = paste0("results/RDS/", comparison, ".Rds"))
  
  # -- Set functions for filtering data -- #
  
  get_top <- function(data, qtl) { 
    data <- copy(data)
    top <- data[, quantile(PBS, qtl)]
    data_top <- data[PBS >= top]
    data_top <- separate_genes(data_top)
    data_top[, N_GENE := .N, by = GENE]
    data_top <- data_top[order(PBS, decreasing = T)]
    data_top <- data_top[data_top[, .I[which.max(PBS)], by = GENE]$V1]
    data_top <- data_top[!grepl("^LINC|^LOC|^MIR\\d", GENE)]
    data_top <- data_top[!grepl("^C.*or.*$", GENE)]
    data_top[, W_ID := NULL]
    
    if("W_GROUP" %in% colnames(data_top))
      data_top[, W_GROUP := NULL]
    
    return(data_top[])
    
  }
  
  get_top_nocons <- function(data, qtl) {
    data <- copy(data)
    top <- data[, quantile(PBS, qtl)]
    data_top <- data[PBS >= top]
    data_top <- exclude_consecutive_windows(data_top)
    data_top <- separate_genes(data_top)
    data_top[, N_GENE := .N, by = GENE]
    data_top <- data_top[order(PBS, decreasing = T)]
    data_top <- data_top[data_top[, .I[which.max(PBS)], by = GENE]$V1]
    data_top <- data_top[!grepl("^LINC|^LOC|^MIR\\d", GENE)]
    data_top <- data_top[!grepl("^C.*or.*$", GENE)]
    data_top[, W_ID := NULL]
    
    if("W_GROUP" %in% colnames(data_top))
      data_top[, W_GROUP := NULL]
    
    return(data_top[])
    
  }
  
  get_top_exons <- function(data, qtl) {
    data <- copy(data)
    top <- data[, quantile(PBS, qtl)]
    data_top <- data[PBS >= top]
    data_top <- separate_genes(data_top)
    data_top <- data_top[FUNCTION == "exonic"]
    data_top <- data_top[order(PBS, decreasing = T)]
    data_top <- data_top[data_top[, .I[which.max(PBS)], by = GENE]$V1]
    data_top[, W_ID := NULL]
    
    if("W_GROUP" %in% colnames(data_top))
      data_top[, W_GROUP := NULL]
    
    return(data_top[])
    
  }
  
  # genes <- get_top(pbs_mean, 0.999)[, GENE]
  # 
  # go.annot <- select(Homo.sapiens, genes, c("ENTREZID", "GO", "TERM"), "SYMBOL") %>% as.data.table()
  # go.annot <- go.annot[, .(ENTREZID, SYMBOL, GO, ONTOLOGY, EVIDENCE, TERM)]
  # go.annot <- go.annot[!(is.na(GO) | is.na(TERM))]
  # go.annot <- go.annot[ONTOLOGY == "BP"]
  # go.annot[, EVIDENCE := NULL]
  # go.annot <- unique(go.annot)
  # 
  # go_all <- go.annot
  # go_by_term <- go.annot[, table(TERM)] %>% as.data.table() %>% .[order(N, decreasing = T)] %>% setnames("N", "FREQ")
  
  xls <- list (
    PBS_MEAN_05 = get_top(pbs_mean, 0.995),
    PBS_MEAN_FILTER_05 = get_top_nocons(pbs_mean, 0.995),
    PBS_MEAN_05_EXONS = get_top_exons(pbs_mean, 0.995),
    PBS_MEAN_01 = get_top(pbs_mean, 0.999),
    PBS_MEAN_FILTER_01 = get_top_nocons(pbs_mean, 0.999),
    PBS_MEAN_01_EXONS = get_top_exons(pbs_mean, 0.999)
    # GO_ALL = go.annot,
    # GO_FREQ_TERMS = go_by_term
  )
  
  # Write tabulated results
  WriteXLS::WriteXLS(xls, paste0("results/XLS/", comparison, ".xls"), row.names = F, AdjWidth = T, BoldHeaderRow = T, AutoFilter = T)
  
}
