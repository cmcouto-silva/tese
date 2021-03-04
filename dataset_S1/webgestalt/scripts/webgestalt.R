#### Settings ####
####################################################################################################!

# WD & libraries
setwd("dataset_S1/webgestalt")
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

## -- Functions -- ##

# Function to get top values
get_top <- function(dt, cutoff, mode = "persnp", pbsw_correction = TRUE) {
  
  dt <- copy(dt)
  
  grp <- c("CHR", "POS", "SNP")
  p_cols <- c("PVALUE", "LOG_PVALUE")
  gene_info_cols <- c("GENE", "GENE_ID", "FUNCTION")
  
  for(method in c("IHS", "XPEHH", "PBS")) {
    if(method %in% names(dt)) setnames(dt, method, "SCORE")
  }
  
  if(mode=="persnp") {
    
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
      mutate(PVALUE = sapply(SCORE, function(value) {1 - sum(value >= .$SCORE) / (nrow(.)+1)})) %>%
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

get_genes <- function(genes) {
  strsplit(genes, ",") %>% unlist %>% unique %>% sort
}

## -- Load Data -- ##

results <- sapply(modes, function(mode) {
  sapply(methods, function(method) {
    read_feather(glue("../results_bin/{mode}/{method}.feather")) %>%
      as.data.table()
  }, simplify = F)
}, simplify = F)

#### ORA ####
####################################################################################################!

for(mode in modes) {
  cat(glue('----  Start Analysis for mode "{mode}"  ----\n\n'))
  cat(paste0(rep('-', 50), collapse = ""), "\n\n")
  for(method in methods) {
    cat(glue('--- Method "{method}" ---\n\n'))
    cat(paste0(rep('-', 50), collapse = ""), "\n\n")
    for(cutoff in cutoff_vals) {
      cat(glue('-- quantile {cutoff} --\n\n'))
      
      outdir <- glue("ORA/{mode}/{method}/pvalue_{cutoff}")
      dir.create(path = outdir, showWarnings = F, recursive = T)
      
      target_genes <- get_top(results[[mode]][[method]], cutoff, mode)[, get_genes(GENE)]
      background_genes <- get_genes(results[[mode]][[method]][["GENE"]])
      
      ora <- tryCatch({
        sapply(gene_sets, function(gs) {
          WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
                      enrichDatabase = gs,
                      interestGene = target_genes,
                      interestGeneType = "genesymbol", 
                      referenceGene = background_genes,
                      referenceGeneType = "genesymbol",
                      minNum = 1,
                      maxNum = 1000, 
                      sigMethod = "top",
                      fdrMethod = "BH",
                      fdrThr = 0.05,
                      topThr = 100,
                      reportNum = 100,
                      perNum = 1000,
                      isOutput = FALSE,
                      outputDirectory = outdir,
                      projectName = gs,
                      dagColor = "continuous",
                      setCoverNum = 10
          ) %>% as.data.table()
        }, simplify = F)
      }, error = function(e) {cat("\n"); NULL})
      saveRDS(ora, glue("{outdir}/pvalue_{cutoff}.RDS"))
    }
    cat("\n")
  }
}

#### GSEA ####
####################################################################################################!

dir.create("GSEA", F)

cat(glue('----------------  GSEA Analysis  ----------------\n\n'))
cat(paste0(rep('-', 50), collapse = ""), "\n\n")

gsea <- sapply(methods, function(method) {

  cat(glue('--- Method "{toupper(method)}" ---\n\n'))
  cat(paste0(rep('-', 50), collapse = ""), "\n\n")

  gsea_input <- results$pergene[[method]][, .(GENE, SCORE)]
  lapply(setdiff(gene_sets, "phenotype_Human_Phenotype_Ontology"), function(gene_set) {

    outdir <- glue("GSEA/{method}/{gene_set}")
    dir.create(outdir, F, T)

    tryCatch({
      WebGestaltR(enrichMethod = "GSEA", organism = "hsapiens",
                  enrichDatabase = gene_set,
                  interestGene = gsea_input,
                  interestGeneType = "genesymbol",
                  collapseMethod = "mean",
                  minNum = 3,
                  maxNum = 1000,
                  sigMethod = "top",
                  topThr = 100,
                  reportNum = 100,
                  perNum = 1000,
                  isOutput = F,
                  outputDirectory = outdir,
                  projectName = gene_set
      ) %>% as.data.table()
    }, error = function(e) NULL)
  }) %>% set_names(setdiff(gene_sets, "phenotype_Human_Phenotype_Ontology"))
}, simplify = F)

# Save GSEA Results
saveRDS(gsea, "GSEA/gsea.RDS")
