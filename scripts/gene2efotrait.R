# GENE Annotation
open_genecards <- function(genes) {
  base_url <- "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
  for(gene in genes) {
    utils::browseURL(paste0(base_url, gene))
  }
}

gene2efotrait <- function(gene) {
  url <- glue::glue (
    "https://www.ebi.ac.uk/gwas/api/search/downloads?q=ensemblMappedGenes:{gene}&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true"
  )
  
  efotraits <- fread(url)[, .(MAPPED_GENE, `REPORTED GENE(S)`, `DISEASE/TRAIT`, MAPPED_TRAIT)]
  if(nrow(efotraits) == 0) return(NA)
  
  check_gene <- function(string, gene) {
    strsplit(string, split = ",|;|-| ") %>%
      sapply(function(genes) gene %in% trimws(genes))
  }  
  
  efotraits <- setDT(tidyr::separate_rows(efotraits[check_gene(MAPPED_GENE, gene)], MAPPED_TRAIT, sep=","))[, .N, MAPPED_TRAIT]
  return(efotraits)
}
