#### SETTINGS ####
####################################################################################################!

# Set working directory
setwd("dataset_S2/fumagwas")
snpsel::pkg.lib(data.table, dplyr, snpsel)

# Set function for reading all sheets from file
read_excel_allsheets <- function(filename, apply_filter = T) {
  sheets <- readxl::excel_sheets(filename) %>%
    sapply(X = ., readxl::read_excel, path = filename, simplify = F) %>%
    lapply(setDT)
  if(apply_filter) {
    sheets <- lapply(sheets, function(dt) dt[Category=="GWAScatalog", .(GeneSet, N_genes, p, adjP)])
  }
  return(sheets)
}

#### MAIN ####
####################################################################################################!

## -- Per SNP -- ##

# PBS
read_excel_allsheets("persnp/pbs_windowed/fumagwas_tbl.xls")
# XP-EHH
read_excel_allsheets("persnp/xpehh/fumagwas_tbl.xls")
# iHS
read_excel_allsheets("persnp/ihs/fumagwas_tbl.xls")

## -- Per Gene -- ##

# PBS
read_excel_allsheets("pergene/pbs_windowed/fumagwas_tbl.xls")
# XP-EHH
read_excel_allsheets("pergene/xpehh/fumagwas_tbl.xls")
# iHS
read_excel_allsheets("pergene/ihs/fumagwas_tbl.xls")

