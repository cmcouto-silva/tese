#### Enrichment Analysis ####
####################################################################################################!

##### KEGG #####
####################################################################################################!

## -- KEGG ORA -- ##

ds2_kegg_ora_persnp_top05 <- fread("dataset_S2/enrichr/results/ihs_persnp_pvalue0.995.txt") %>%
  filter(Gene_set == "KEGG_2019_Human" & `Adjusted P-value` <= 0.05) %>%
  select(Term, Overlap, `Combined Score`, `P-value`, `Adjusted P-value`) %>%
  rename(c("Score" = "Combined Score", "-log10(p)" = "P-value", "FDR" = "Adjusted P-value")) %>%
  mutate(`-log10(p)` = -log10(`-log10(p)`)) %>%
  mutate_if(is.numeric, round, 4)

ds2_kegg_ora_persnp_top1 <- fread("dataset_S2/enrichr/results/ihs_persnp_pvalue0.99.txt") %>%
  filter(Gene_set == "KEGG_2019_Human" & `Adjusted P-value` <= 0.05) %>%
  select(Term, Overlap, `Combined Score`, `P-value`, `Adjusted P-value`) %>%
  rename(c("Score" = "Combined Score", "-log10(p)" = "P-value", "FDR" = "Adjusted P-value")) %>%
  mutate(`-log10(p)` = -log10(`-log10(p)`)) %>%
  mutate_if(is.numeric, round, 4)

ds2_kegg_ora_persnp_tbl <- full_join(ds2_kegg_ora_persnp_top05, ds2_kegg_ora_persnp_top1, by = "Term") %>%
  select(-Overlap.x, -Overlap.y) %>%
  mutate_if(is.numeric, ~ replace(., is.na(.), "-"))

#names(ds2_kegg_ora_persnp_tbl) <- gsub(".x$|.y$", "", names(ds2_kegg_ora_persnp_tbl))
names(ds2_kegg_ora_persnp_tbl) <- c("Termo", rep(c("Score", "-log10(p)", "FDR"), 2))

ds2_kegg_ora_persnp_kbl <- kbl(ds2_kegg_ora_persnp_tbl, booktabs = T, format = "latex", linesep = "", align = c("l", rep("r", 6)), escape = F) %>%
  kable_styling(latex_options = c("striped", "scale_down")) %>%
  add_header_above(c(" ", "Top 0,5%" = 3, "Top 1%" = 3)) %>%
  add_header_above(c(" ", "iHS" = 6))

save_kable(ds2_kegg_ora_persnp_kbl, "tables/tex/ds2_kegg_ora_persnp.tex")

## -- KEGG GSEA -- ##

ds1_gsea <- readRDS("dataset_S1/webgestalt/GSEA/gsea.RDS")
ds2_gsea <- readRDS("dataset_S2/webgestalt/GSEA/gsea.RDS")

ds1_gsea_xpehh <- ds1_gsea$xpehh$pathway_KEGG %>%
  filter(FDR <= 0.05) %>%
  select(description, normalizedEnrichmentScore, pValue, FDR) %>%
  mutate_if(is.numeric, round, 4) %>%
  set_names(c("Via de Sinalização", "Enriquecimento", "p-valor", "FDR"))

gsea_kbl <- kbl(ds1_gsea_xpehh, booktabs = T, format = "latex", linesep = "", escape = F) %>%
  kable_classic()

save_kable(gsea_kbl, "tables/tex/gsea.tex")

##### Gene Ontology (GO) #####
####################################################################################################!

# Common function
go_filter <- function(dt) {
  dt %>%
    filter(NS == "BP", depth != 0, p_fdr <= 0.05) %>%
    mutate(p_uncorrected = -log10(p_uncorrected)) %>%
    select(name, p_uncorrected, p_fdr) %>%
    set_names(c("Termo", "-log10(p)", "FDR"))
}

## -- Dataset S1 ─ GO All -- ##

ds1_go_all_persnp_xpehh_0.99 <- fread("dataset_S1/goatools/results/go_all_persnp_xpehh_0.99.csv") %>%
  go_filter() %>% mutate(`Método` = "XP-EHH") %>% relocate(`Método`)

ds1_go_all_persnp_ihs_0.99 <- fread("dataset_S1/goatools/results/go_all_persnp_ihs_0.99.csv") %>%
  go_filter() %>% mutate(`Método` = "iHS") %>% relocate(`Método`)

ds1_go_all_persnp <- rbind(ds1_go_all_persnp_xpehh_0.99, ds1_go_all_persnp_ihs_0.99)

# Removing terms with high hierarchy (non-sense)
ds1_go_all_persnp[!Termo %in% c(
  "regulation of system process", "system process", "multicellular organismal process",
  "regulation of biological quality", "regulation of localization", "multicellular organismal signaling")]

ds1_go_all_persnp_kbl <- kbl(ds1_go_all_persnp, booktabs = T, format = "latex", linesep = "", escape = F) %>%
  kable_styling() %>%
  collapse_rows(1, latex_hline = "major") %>%
  row_spec(1:nrow(ds1_go_all_persnp_xpehh_0.99) - 1, extra_latex_after = "\\rowcolor{gray!6}")

save_kable(ds1_go_all_persnp_kbl, "tables/tex/ds1_go_all_persnp.tex")

## -- Dataset S2 ─ GO All -- ##

ds2_go_all_pergene_xpehh_0.995 <- fread("dataset_S2/goatools/results/go_all_pergene_xpehh_0.995.csv") %>%
  go_filter() %>% mutate(`Método` = "XP-EHH (top 0,5%)") %>% relocate(`Método`)

ds2_go_all_pergene_xpehh_0.99 <- fread("dataset_S2/goatools/results/go_all_pergene_xpehh_0.99.csv") %>%
  go_filter() %>% mutate(`Método` = "XP-EHH (top 1%)") %>% relocate(`Método`)

ds2_go_all_pergene <- rbind(ds2_go_all_pergene_xpehh_0.995, ds2_go_all_pergene_xpehh_0.99)

ds2_go_all_pergene_kbl <- kbl(ds2_go_all_pergene, booktabs = T, format = "latex", linesep = "", escape = F) %>%
  kable_classic() %>%
  collapse_rows(1, latex_hline = "major") %>%
  row_spec(1:nrow(ds2_go_all_pergene_xpehh_0.995) - 1, extra_latex_after = "\\rowcolor{gray!6}")

save_kable(ds2_go_all_pergene_kbl, "tables/tex/ds2_go_all_pergene.tex")

## -- Dataset S2 ─ GO Immune -- ##

ds2_go_immune_pergene_pbs_0.995 <- fread("dataset_S2/goatools/results/go_immune_pergene_pbs_windowed_0.99.csv") %>%
  go_filter()

ds2_go_immune_pergene_xpehh_0.99 <- fread("dataset_S2/goatools/results/go_immune_pergene_xpehh_0.99.csv") %>%
  go_filter()

ds2_go_immune_pergene <- rbind(ds2_go_immune_pergene_pbs_0.995, ds2_go_immune_pergene_xpehh_0.99)

ds2_go_immune_pergene_kbl <- kbl(ds2_go_immune_pergene, booktabs = T, format = "latex", linesep = "", escape = F) %>%
  kable_classic() %>%
  pack_rows("PBS", 1, nrow(ds2_go_immune_pergene_pbs_0.995)) %>%
  pack_rows("XP-EHH", nrow(ds2_go_immune_pergene_pbs_0.995)+1, nrow(ds2_go_immune_pergene_xpehh_0.99))

save_kable(ds2_go_immune_pergene_kbl, "tables/tex/ds2_go_immune_pergene.tex")

##### FUMAGWAS #####
####################################################################################################!

## -- Dataset S1 -- ## 

ds1_fumagwas_ihs_persnp_tbl <- readxl::read_xls("dataset_S1/fumagwas/persnp/ihs/fumagwas_tbl.xls") %>%
  filter(Category=="GWAScatalog") %>%
  mutate(p = -log10(p)) %>%
  mutate_if(is.numeric, round, 4) %>%
  select(GeneSet, p, adjP) %>%
  rename("Fenótipo" = "GeneSet", "-log10(p)" = "p", "FDR" = "adjP")

ds1_fumagwas_ihs_persnp_kbl <- kbl(ds1_fumagwas_ihs_persnp_tbl, booktabs = T, format = "latex", linesep = "") %>%
  kable_classic()

save_kable(ds1_fumagwas_ihs_persnp_kbl, "tables/tex/ds1_fumagwas_ihs_persnp.tex")

## -- Dataset S2 -- ## 

# Per SNP
ds2_fumagwas_pbs_persnp_tbl <- readxl::read_xls("dataset_S2/fumagwas/persnp/pbs_windowed/fumagwas_tbl.xls", "top01") %>%
  mutate(p = -log10(p)) %>%
  mutate_if(is.numeric, round, 4) %>%
  select(GeneSet, p, adjP) %>%
  rename("Fenótipo" = "GeneSet", "-log10(p)" = "p", "FDR" = "adjP")

ds2_fumagwas_pbs_persnp_kbl <- kbl(ds2_fumagwas_pbs_persnp_tbl, booktabs = T, format = "latex", linesep = "") %>%
  kable_classic()

save_kable(ds2_fumagwas_pbs_persnp_kbl, "tables/tex/ds2_fumagwas_pbs_persnp.tex")

# Per gene
ds2_fumagwas_pbs_pergene_tbl <- readxl::read_xls("dataset_S2/fumagwas/pergene/pbs_windowed/fumagwas_tbl.xls", "top01") %>%
  mutate(p = -log10(p)) %>%
  mutate_if(is.numeric, round, 4) %>%
  select(GeneSet, p, adjP) %>%
  rename("Fenótipo" = "GeneSet", "-log10(p)" = "p", "FDR" = "adjP")

ds2_fumagwas_pbs_pergene_kbl <- kbl(ds2_fumagwas_pbs_pergene_tbl, booktabs = T, format = "latex", linesep = "") %>%
  kable_classic()

save_kable(ds2_fumagwas_pbs_pergene_kbl, "tables/tex/ds2_fumagwas_pbs_pergene.tex")

##### Metanalysis #####
####################################################################################################!

metanalysis_table <- fread("metanalysis/study_pop_methods.csv") %>%
  rename(Autores=Authors, Populações=Populations, `Métodos de seleção`=`Selection Methods`) %>%
  mutate(Populações=gsub("others", "outras", Populações)) %>%
  #relocate(PMID, .after = last_col())
  select(-PMID)

metanalysis_table_kbl <- kbl(metanalysis_table, booktabs = T, format = "latex", linesep = "") %>%
  kable_styling(latex_options = c("striped"))

save_kable(metanalysis_table_kbl, "tables/tex/metanalysis.tex")

metanalysis_ensembl <- readxl::read_xlsx("metanalysis/efotable_ensembl.xlsx")

metanalysis_ensembl <- metanalysis_ensembl %>%
  select(Phenotype, n_genes, n_pmid, pbs_s1:ihs_s2) %>%
  arrange(desc(n_genes), Phenotype) %>%
  head(30)

sel_methods <- c("PBS", "XP-EHH", "iHS")
names(metanalysis_ensembl) <- c("Fenótipo", "Genes", "Estudos", rep(sel_methods, 2))

metanalysis_ensembl_kbl <- kbl(metanalysis_ensembl, booktabs = T, format = "latex", linesep = "") %>%
  kable_styling(latex_options = c("striped", "scale_down")) %>%
  add_header_above(c(" " = 3, "Dataset S1" = 3, "Dataset S2" = 3))

save_kable(metanalysis_ensembl_kbl, "tables/tex/metanalysis_ensembl.tex")

#### Spreadsheets ####
####################################################################################################!

names(ds2_kegg_ora_persnp_tbl)[2:length(ds2_kegg_ora_persnp_tbl)] <- paste(names(ds2_kegg_ora_persnp_tbl)[2:4], rep(c("top0.5%", "top1%"), each=3), sep = "_")

enrichment_data <- list (
  ds2_kegg_ora_ihs_persnp = ds2_kegg_ora_persnp_tbl,
  ds1_kegg_gsea_xpehh = ds1_gsea_xpehh,
  ds1_go_all_xpehh_ihs_persnp = ds1_go_all_persnp,
  ds2_go_immune_pbs_xpehh_pergene = ds2_go_immune_pergene,
  ds1_fumagwas_ihs_persnp_tbl = ds1_fumagwas_ihs_persnp_tbl,
  ds2_fumagwas_pbs_persnp_tbl = ds2_fumagwas_pbs_persnp_tbl,
  ds2_fumagwas_pbs_pergene_tbl = ds2_fumagwas_pbs_pergene_tbl
)

WriteXLS::WriteXLS(enrichment_data, "tables/xlsx/enrichment.xls", row.names = F, AdjWidth = T, BoldHeaderRow = T)
