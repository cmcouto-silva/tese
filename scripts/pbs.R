#### 1 PBS ####
####################################################################################################!

##### 1.1 TABLES #####
####################################################################################################!

###### 1.1.1 DATASET S1 ######

# PBS-windowed per SNP
ds1_pbsw_persnp_tbl <- ds1_persnp_top$pbs_windowed %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

length(ds1_pbsw_persnp_tbl$GENE)
ds1_pbsw_persnp_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

# Translate to portuguese
ds1_pbsw_persnp_tbl <- rename(ds1_pbsw_persnp_tbl, c("FUNÇÃO"="FUNCTION"))

ds1_pbsw_persnp_kbl <- kbl(ds1_pbsw_persnp_tbl, booktabs = T, format = "latex", linesep = "") %>%
  kable_styling(latex_options = c("scale_down")) # "striped"

save_kable(ds1_pbsw_persnp_kbl, "tables/tex/ds1_pbsw_persnp.tex")

wilcox.test(ds1_persnp$pbs_windowed$SCORE, mu = min(ds1_persnp_top$pbs_windowed$SCORE))

# PBS-windowed per Gene
ds1_pbsw_pergene_tbl <- ds1_pergene_top$pbs_windowed %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "POS", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds1_pbsw_pergene_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

ds1_pbsw_pergene_kbl <- kbl(ds1_pbsw_pergene_tbl, booktabs = T, format = "latex", linesep = "")
save_kable(ds1_pbsw_pergene_kbl, "tables/tex/ds1_pbsw_pergene.tex")

wilcox.test(ds1_pergene$pbs_windowed$SCORE, mu = min(ds1_pergene_top$pbs_windowed$SCORE))

###### 1.1.2 DATASET S2 ######

# PBS-windowed per SNP
ds2_pbsw_persnp_tbl <- ds2_persnp_top$pbs_windowed %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds2_pbsw_persnp_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

# Translate to portuguese
ds2_pbsw_persnp_tbl <- rename(ds2_pbsw_persnp_tbl, c("FUNÇÃO"="FUNCTION"))

ds2_pbsw_persnp_kbl <- kbl(ds2_pbsw_persnp_tbl, booktabs = T, format = "latex", linesep = "") %>%
  kable_styling(latex_options = c("scale_down")) # "striped"

save_kable(ds2_pbsw_persnp_kbl, "tables/tex/ds2_pbsw_persnp.tex")

wilcox.test(ds2_persnp$pbs_windowed$SCORE, mu = min(ds2_persnp_top$pbs_windowed$SCORE))

# PBS-windowed per Gene
ds2_pbsw_pergene_tbl <- ds2_pergene_top$pbs_windowed %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "POS", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds2_pbsw_pergene_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

ds2_pbsw_pergene_kbl <- kbl(ds2_pbsw_pergene_tbl, booktabs = T, format = "latex", linesep = "")

save_kable(ds2_pbsw_pergene_kbl, "tables/tex/ds2_pbsw_pergene.tex")

wilcox.test(ds2_pergene$pbs_windowed$SCORE, mu = min(ds2_pergene_top$pbs_windowed$SCORE))

###### 1.1.3 AMBOS DATASETS ######

# PBS Windowed
pbsw_pergene_tbl <- rbind (
  data.table(Dataset = "Dataset S1", ds1_pbsw_pergene_tbl),
  data.table(Dataset = "Dataset S2", ds2_pbsw_pergene_tbl)
)

pbsw_pergene_kbl <- kbl(pbsw_pergene_tbl, booktabs = T, format = "latex") %>%
  kable_styling() %>%
  collapse_rows(1, latex_hline = "major") %>%
  row_spec(1:nrow(ds1_pbsw_pergene_tbl) - 1, extra_latex_after = "\\rowcolor{gray!6}")

save_kable(pbsw_pergene_kbl, "tables/tex/pbsw_pergene.tex") 

# Write XLS

pbsw_data <- list (
  ds1_top0.01_snp = ds1_pbsw_persnp_tbl,
  ds1_top0.01_gene = ds1_pbsw_pergene_tbl,
  ds2_top0.01_snp = ds2_pbsw_persnp_tbl,
  ds2_top0.01_gene = ds2_pbsw_pergene_tbl
)

WriteXLS::WriteXLS(pbsw_data, "tables/xlsx/pbs.xls", row.names = F, AdjWidth = T, BoldHeaderRow = T)

##### 1.2 FIGURES #####
####################################################################################################!

######  1.2.1 DATASET S1  ######

## -- PBS Per SNP -- ##

ds1_pbs_persnp_plot <- manhattan_plot (
  ds1_persnp$pbs_windowed, annot = "GENE", score = "LOG_PVALUE", 
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"),
  plot.title = "Dataset S1 - PBS por SNP", plot.title.size = 18,
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds1_pbs_persnp_plot <- annotate_plot (
  ds1_pbs_persnp_plot, y.upper.limit = 5.25, threshold = .999,
  nudge_y = 0.5, size=4.2, alpha=.8, segment.alpha = .5
)

ggsave("figures/png/ds1_pbsw_persnp.png", ds1_pbs_persnp_plot, height = 6, width = 11)
ggsave("figures/pdf/ds1_pbsw_persnp.pdf", ds1_pbs_persnp_plot, height = 6, width = 11)

## -- Plot Per Gene -- ##

ds1_pbs_pergene_plot <- manhattan_plot (
  ds1_pergene$pbs_windowed, snp= NULL, annot = "GENE", score = "LOG_PVALUE",
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"),
  plot.title = "Dataset S1 - PBS por Gene", plot.title.size = 18,
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds1_pbs_pergene_plot <- annotate_plot (
  ds1_pbs_pergene_plot, threshold = .999, y.upper.limit = 4.3,
  nudge_y = 0.3, size=4.2, alpha=.8, segment.alpha = .5
)

ggsave("figures/png/ds1_pbsw_pergene.png", ds1_pbs_pergene_plot, height = 6, width = 11)
ggsave("figures/pdf/ds1_pbsw_pergene.pdf", ds1_pbs_pergene_plot, height = 6, width = 11)

## -- Plot SNP & Gene -- ##

# Dataset S1 - PBS - snp vs gene
ds1_pbs_title <- ggdraw() + draw_label("Dataset S1 - PBS", fontface='bold', size = 21)

ds1_pbs_grid <- plot_grid (
  ds1_pbs_persnp_plot + labs(x='cromossomo', title=""), ds1_pbs_pergene_plot + labs(x='cromossomo', title=''), 
  labels = "AUTO", label_size = 18, nrow = 2, ncol = 1, rel_heights = c(.55, .45)
)

ds1_pbs_grid <- plot_grid(ds1_pbs_title, ds1_pbs_grid, ncol=1, rel_heights=c(0.08, 1))

ggsave("figures/png/ds1_pbsw_both.png", height = 11, width = 11, dpi=300)
ggsave("figures/pdf/ds1_pbsw_both.pdf", height = 11, width = 11, dpi=300)


######  1.2.2 DATASET S2  ######

## -- Plot Per SNP -- ##

ds2_pbs_persnp_plot <- manhattan_plot (
  ds2_persnp$pbs_windowed, annot = "GENE", score = "LOG_PVALUE", 
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"),
  plot.title = "Dataset S2 - PBS por SNP", plot.title.size = 18,
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds2_pbs_persnp_plot <- annotate_plot (
  ds2_pbs_persnp_plot, threshold = .999, y.upper.limit = 5.25,
  nudge_y = 0.5, size=4.2, alpha=.8, segment.alpha = .5,
  include_genes = "ACACA",
)

ggsave("figures/png/ds2_pbsw_persnp.png", ds2_pbs_persnp_plot, height = 6, width = 11)
ggsave("figures/pdf/ds2_pbsw_persnp.pdf", ds2_pbs_persnp_plot, height = 6, width = 11)

## -- Plot Per Gene -- ##

ds2_pbs_pergene_plot <- manhattan_plot (
  ds2_pergene$pbs_windowed, snp= NULL, annot = "GENE", score = "LOG_PVALUE",
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"),
  plot.title = "Dataset S2 - PBS por Gene", plot.title.size = 18,
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds2_pbs_pergene_plot <- annotate_plot (
  ds2_pbs_pergene_plot, threshold = .999, y.upper.limit = 4.4,
  nudge_y=0.3, size=4.2, alpha=.8, segment.alpha=.5
)

ggsave("figures/png/ds2_pbsw_pergene.png", ds2_pbs_pergene_plot, height = 6, width = 11)
ggsave("figures/pdf/ds2_pbsw_pergene.pdf", ds2_pbs_pergene_plot, height = 6, width = 11)

## -- Plot SNP & Gene -- ##

# Dataset S2 - PBS - snp vs gene
ds2_pbs_title <- ggdraw() + draw_label("Dataset S2 - PBS", fontface='bold', size = 21)

ds2_pbs_grid <- plot_grid (
  ds2_pbs_persnp_plot + labs(x='cromossomo', title=""), ds2_pbs_pergene_plot + labs(x='cromossomo', title=''), 
  labels = "AUTO", label_size = 18, nrow = 2, ncol = 1, rel_heights = c(.55, .45)
)

ds2_pbs_grid <- plot_grid(ds2_pbs_title, ds2_pbs_grid, ncol=1, rel_heights=c(0.08, 1))

ggsave("figures/png/ds2_pbsw_both.png", ds2_pbs_grid, height = 11, width = 11, dpi=300)
ggsave("figures/pdf/ds2_pbsw_both.pdf", ds2_pbs_grid, height = 11, width = 11, dpi=300)

###########################################  BREAK  ################################################!

#### 2 XP-EHH ####
####################################################################################################!

##### 2.1 TABLES #####
####################################################################################################!

######  2.1.1 DATASET S1  ######

######  2.1.2 DATASET S2  ######

##### 2.2 FIGURES #####
####################################################################################################!

######  2.2.1 DATASET S1  ######

######  2.2.2 DATASET S2  ######

############################################  STOP  ################################################!


#### 3 iHS ####
####################################################################################################!

##### 3.1 TABLES #####
####################################################################################################!

######  3.1.1 DATASET S1  ######

######  3.1.2 DATASET S2  ######

##### 3.2 FIGURES #####
####################################################################################################!

######  3.2.1 DATASET S1  ######

######  3.2.2 DATASET S2  ######
