#### iHS ####
####################################################################################################!

#### 1. TABLES ####
####################################################################################################!

######  2.1 DATASET S1  ######

# iHS per SNP
ds1_ihs_persnp_tbl <- ds1_persnp_top$ihs %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds1_ihs_persnp_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

# Translate to portuguese
ds1_ihs_persnp_tbl <- rename(ds1_ihs_persnp_tbl, c("FUNÇÃO"="FUNCTION"))

ds1_ihs_persnp_kbl <- kbl(ds1_ihs_persnp_tbl, booktabs = T, format = "latex", linesep = "") %>%
  kable_styling(latex_options = c("scale_down")) # "striped"

save_kable(ds1_ihs_persnp_kbl, "tables/tex/ds1_ihs_persnp.tex")

wilcox.test(ds1_persnp$ihs$SCORE, mu = min(ds1_persnp_top$ihs$SCORE))

# iHS per Gene
ds1_ihs_pergene_tbl <- ds1_pergene_top$ihs %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds1_ihs_pergene_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

ds1_ihs_pergene_kbl <- kbl(ds1_ihs_pergene_tbl, booktabs = T, format = "latex", linesep = "")

save_kable(ds1_ihs_pergene_kbl, "tables/tex/ds1_ihs_pergene.tex")

wilcox.test(ds1_pergene$ihs$SCORE, mu = min(ds1_pergene_top$ihs$SCORE))

######  2.2 DATASET S2  ######

# iHS per SNP
ds2_ihs_persnp_tbl <- ds2_persnp_top$ihs %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds2_ihs_persnp_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

# Translate to portuguese
ds2_ihs_persnp_tbl <- rename(ds2_ihs_persnp_tbl, c("FUNÇÃO"="FUNCTION"))

ds2_ihs_persnp_kbl <- kbl(ds2_ihs_persnp_tbl, booktabs = T, format = "latex", linesep = "") %>%
  kable_styling(latex_options = c("scale_down")) # "striped"

save_kable(ds2_ihs_persnp_kbl, "tables/tex/ds2_ihs_persnp.tex")

wilcox.test(ds2_persnp$ihs$SCORE, mu = min(ds2_persnp_top$ihs$SCORE))

# iHS per Gene
ds2_ihs_pergene_tbl <- ds2_pergene_top$ihs %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds2_ihs_pergene_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

ds2_ihs_pergene_kbl <- kbl(ds2_ihs_pergene_tbl, booktabs = T, format = "latex", linesep = "")

save_kable(ds2_ihs_pergene_kbl, "tables/tex/ds2_ihs_pergene.tex")

wilcox.test(ds2_pergene$ihs$SCORE, mu = min(ds2_pergene_top$ihs$SCORE))

######  2.3 XLS  ######

ihs_data <- list (
  ds1_top0.01_snp = ds1_ihs_persnp_tbl,
  ds1_top0.01_gene = ds1_ihs_pergene_tbl,
  ds2_top0.01_snp = ds2_ihs_persnp_tbl,
  ds2_top0.01_gene = ds2_ihs_pergene_tbl
)

WriteXLS::WriteXLS(ihs_data, "tables/xlsx/ihs.xls", row.names = F, AdjWidth = T, BoldHeaderRow = T)

#### 2 FIGURES ####
####################################################################################################!

######  2.1 DATASET S1  ######

## -- iHS Per SNP -- ##

ds1_ihs_persnp_plot <- manhattan_plot (
  ds1_persnp$ihs, annot = "GENE", score = "SCORE",
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"), plot.title.size = 18,
  plot.title = "Dataset S1 - iHS por SNP", y.axis.title = "abs(iHS)",
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds1_ihs_persnp_plot <- annotate_plot (
  ds1_ihs_persnp_plot, threshold = .999, nudge_y = 0.01,
  size = 4.2, alpha = .8, segment.alpha = .5
)

ggsave("figures/png/ds1_ihs_persnp.png", ds1_ihs_persnp_plot, height = 6, width = 11)
ggsave("figures/pdf/ds1_ihs_persnp.pdf", ds1_ihs_persnp_plot, height = 6, width = 11)

## -- iHS Per Gene -- ##

ds1_ihs_pergene_plot <- manhattan_plot (
  ds1_pergene$ihs, snp= NULL, annot = "GENE", score = "SCORE",
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"),
  plot.title = "Dataset S1 - iHS por Gene", plot.title.size = 18,
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  y.axis.title = "abs(iHS)", x.axis.title = "cromossomo"
)

ds1_ihs_pergene_plot <- annotate_plot (
  ds1_ihs_pergene_plot, threshold = .999, y.upper.limit = 4.3,
  nudge_y = 0.05, size=4.2, alpha=.8, segment.alpha = .5
)

ggsave("figures/png/ds1_ihs_pergene.png", ds1_ihs_pergene_plot, height = 6, width = 11)
ggsave("figures/pdf/ds1_ihs_pergene.pdf", ds1_ihs_pergene_plot, height = 6, width = 11)

######  2.2 DATASET S2  ######

## -- Plot Per SNP -- ##

ds2_ihs_persnp_plot <- manhattan_plot (
  ds2_persnp$ihs, annot = "GENE", score = "SCORE", 
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"), plot.title.size = 18,
  plot.title = "Dataset S2 - iHS por SNP", y.axis.title = "abs(iHS)",
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds2_ihs_persnp_plot <- annotate_plot (
  ds2_ihs_persnp_plot, threshold = .999, nudge_y = 0.01, y.upper.limit = 7.6,
  size = 4.2, alpha = .8, segment.alpha = .5
)

ggsave("figures/png/ds2_ihs_persnp.png", ds2_ihs_persnp_plot, height = 6, width = 11)
ggsave("figures/pdf/ds2_ihs_persnp.pdf", ds2_ihs_persnp_plot, height = 6, width = 11)

## -- Plot Per Gene -- ##

ds2_ihs_pergene_plot <- manhattan_plot (
  ds2_pergene$ihs, snp= NULL, annot = "GENE", score = "SCORE",
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"),
  plot.title = "Dataset S2 - iHS per Gene", plot.title.size = 18,
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  y.axis.title = "abs(iHS)", x.axis.title = "cromossomo"
)

ds2_ihs_pergene_plot <- annotate_plot (
  ds2_ihs_pergene_plot, threshold = .999, y.upper.limit = 4,
  nudge_y = 0.05, size=4.2, alpha=.8, segment.alpha = .5
)

ggsave("figures/png/ds2_ihs_pergene.png", ds2_ihs_pergene_plot, height = 6, width = 11)
ggsave("figures/pdf/ds2_ihs_pergene.pdf", ds2_ihs_pergene_plot, height = 6, width = 11)

