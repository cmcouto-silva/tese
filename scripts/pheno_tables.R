#### MAIN ####
####################################################################################################!

# Load data
ds1_files <- list.files("dataset_S1/pheno_annot/results/", full.names = T)
ds1_pheno_tables <- lapply(ds1_files, fread) %>%
  set_names(sapply(ds1_files, function(x) rm.extension(basename(x))))

ds2_files <- list.files("dataset_S2/pheno_annot/results/", full.names = T)
ds2_pheno_tables <- lapply(ds2_files, fread) %>%
  set_names(sapply(ds2_files, function(x) rm.extension(basename(x))))

# Define helper functions
make_table <- function(...) {
  dt_list <- list(...)
  max_rows = max(sapply(dt_list, nrow))
  dt_list <- sapply(dt_list, function(dt) {
    dt <- dt[Phenotype != "Annotated by HGMD"][order(-N, Phenotype)]
    dt <- dt[1:max_rows, .(Phenotype, N)]
    setnames(dt, "Phenotype", "Fenótipo")
    dt[]
  }, simplify = F)
  do.call(cbind, dt_list)
}

make_kable <- function(dt) {
  kbl(dt, booktabs = T, format = "latex", linesep = "") %>%
    kable_styling(latex_options = c("striped", "scale_down")) %>%
    add_header_above(c("PBS" = 2, "XP-EHH" = 2, "iHS" = 2))
}

######  DATASET S1  ######

## -- Per SNP -- ##

# Ensembl
ds1_ensembl_persnp <- make_table (
  ds1_pheno_tables$ensembl_persnp_pbs_windowed,
  ds1_pheno_tables$ensembl_persnp_xpehh,
  ds1_pheno_tables$ensembl_persnp_ihs
  ) %>% head(30)

# GWAS
ds1_gwas_persnp <- make_table (
  ds1_pheno_tables$gwas_persnp_pbs_windowed,
  ds1_pheno_tables$gwas_persnp_xpehh,
   ds1_pheno_tables$gwas_persnp_ihs
) %>% head(30)

## -- Per Gene -- ##

# Ensembl
ds1_ensembl_pergene <- make_table (
  ds1_pheno_tables$ensembl_pergene_pbs_windowed,
  ds1_pheno_tables$ensembl_pergene_xpehh,
  ds1_pheno_tables$ensembl_pergene_ihs
) %>% head(30)

# GWAS
ds1_gwas_pergene <- make_table (
  ds1_pheno_tables$gwas_pergene_pbs_windowed,
  ds1_pheno_tables$gwas_pergene_xpehh,
  ds1_pheno_tables$gwas_pergene_ihs
) %>% head(30)

######  DATASET S2  ######

## -- Per SNP -- ##

# Ensembl
ds2_ensembl_persnp <- make_table (
  ds2_pheno_tables$ensembl_persnp_pbs_windowed,
  ds2_pheno_tables$ensembl_persnp_xpehh,
  ds2_pheno_tables$ensembl_persnp_ihs
) %>% head(30)

# GWAS
ds2_gwas_persnp <- make_table (
  ds2_pheno_tables$gwas_persnp_pbs_windowed,
  ds2_pheno_tables$gwas_persnp_xpehh,
  ds2_pheno_tables$gwas_persnp_ihs
) %>% head(30)

## -- Per Gene -- ##

# Ensembl
ds2_ensembl_pergene <- make_table (
  ds2_pheno_tables$ensembl_pergene_pbs_windowed,
  ds2_pheno_tables$ensembl_pergene_xpehh,
  ds2_pheno_tables$ensembl_pergene_ihs
) %>% head(30)

# GWAS
ds2_gwas_pergene <- make_table (
  ds2_pheno_tables$gwas_pergene_pbs_windowed,
  ds2_pheno_tables$gwas_pergene_xpehh,
  ds2_pheno_tables$gwas_pergene_ihs
) %>% head(30)


#### WRITE TABLES ####
####################################################################################################!

# Dataset S1

ds1_ensembl_persnp_tbl <- make_kable(ds1_ensembl_persnp)
save_kable(ds1_ensembl_persnp_tbl, "tables/tex/ds1_ensembl_persnp.tex")

ds1_gwas_persnp_tbl <- make_kable(ds1_gwas_persnp)
save_kable(ds1_gwas_persnp_tbl, "tables/tex/ds1_gwas_persnp.tex")

ds1_ensembl_pergene_tbl <- make_kable(ds1_ensembl_pergene)
save_kable(ds1_ensembl_pergene_tbl, "tables/tex/ds1_ensembl_pergene.tex")

ds1_gwas_pergene_tbl <- make_kable(ds1_gwas_pergene)
save_kable(ds1_gwas_pergene_tbl, "tables/tex/ds1_gwas_pergene.tex")

# Dataset S2

ds2_ensembl_persnp_tbl <- make_kable(ds2_ensembl_persnp)
save_kable(ds2_ensembl_persnp_tbl, "tables/tex/ds2_ensembl_persnp.tex")

ds2_gwas_persnp_tbl <- make_kable(ds2_gwas_persnp)
save_kable(ds2_gwas_persnp_tbl, "tables/tex/ds2_gwas_persnp.tex")

ds2_ensembl_pergene_tbl <- make_kable(ds2_ensembl_pergene)
save_kable(ds2_ensembl_pergene_tbl, "tables/tex/ds2_ensembl_pergene.tex")

ds2_gwas_pergene_tbl <- make_kable(ds2_gwas_pergene)
save_kable(ds2_gwas_pergene_tbl, "tables/tex/ds2_gwas_pergene.tex")

# XLSX

## Save tables as .xlsx
dir.create("tables/xlsx", F, T)

mode <- c("persnp", "pergene")
database <- c("gwas", "ensembl")
method <- c("pbs_windowed", "xpehh", "ihs")

table_order <- apply(expand.grid(database, mode, method), 1, paste, collapse="_")

# Function to write Excel file
write_xlsx <- function(data, out) {
  wb <- createWorkbook()
  
  for(sheet in table_order) {
    addWorksheet(wb, sheet)
    writeData(wb, sheet, data[[sheet]], headerStyle = createStyle(textDecoration = "Bold"))
    addStyle(wb = wb, sheet = sheet, style = createStyle(halign = "right", textDecoration = "Bold"), rows = 1, cols = 3, T)
    setColWidths(wb, sheet, cols = 1:3, widths = c(70, 35, 10))
  }
  saveWorkbook(wb, out, overwrite = TRUE)
}

write_xlsx(ds1_pheno_tables, "tables/xlsx/ds1_pheno_tables.xlsx")
write_xlsx(ds2_pheno_tables, "tables/xlsx/ds2_pheno_tables.xlsx")
