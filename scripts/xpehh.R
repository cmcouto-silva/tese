#### XP-EHH ####
####################################################################################################!

#### 1. TABLES ####
####################################################################################################!

######  2.1 DATASET S1  ######

# XP-EHH per SNP
ds1_xpehh_persnp_tbl <- ds1_persnp_top$xpehh %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds1_xpehh_persnp_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% sort %>% length

# Translate to portuguese
ds1_xpehh_persnp_tbl <- rename(ds1_xpehh_persnp_tbl, c("FUNÇÃO"="FUNCTION"))

ds1_xpehh_persnp_kbl <- kbl(ds1_xpehh_persnp_tbl, booktabs = T, format = "latex", linesep = "") %>%
  kable_styling(latex_options = c("scale_down")) # "striped"

save_kable(ds1_xpehh_persnp_kbl, "tables/tex/ds1_xpehh_persnp.tex")

wilcox.test(ds1_persnp$xpehh$SCORE, mu = min(ds1_persnp_top$xpehh$SCORE))

# XP-EHH per Gene
ds1_xpehh_pergene_tbl <- ds1_pergene_top$xpehh %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds1_xpehh_pergene_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

ds1_xpehh_pergene_kbl <- kbl(ds1_xpehh_pergene_tbl, booktabs = T, format = "latex", linesep = "")

save_kable(ds1_xpehh_pergene_kbl, "tables/tex/ds1_xpehh_pergene.tex")

wilcox.test(ds1_pergene$xpehh$SCORE, mu = min(ds1_pergene_top$xpehh$SCORE))

######  2.2 DATASET S2  ######

# XP-EHH per SNP
ds2_xpehh_persnp_tbl <- ds2_persnp_top$xpehh %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds2_xpehh_persnp_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

# Translate to portuguese
ds2_xpehh_persnp_tbl <- rename(ds2_xpehh_persnp_tbl, c("FUNÇÃO"="FUNCTION"))

ds2_xpehh_persnp_kbl <- kbl(ds2_xpehh_persnp_tbl, booktabs = T, format = "latex", linesep = "") %>%
  kable_styling(latex_options = c("scale_down")) # "striped"

save_kable(ds2_xpehh_persnp_kbl, "tables/tex/ds2_xpehh_persnp.tex")

wilcox.test(ds2_persnp$xpehh$SCORE, mu = min(ds2_persnp_top$xpehh$SCORE))

# XP-EHH per Gene
ds2_xpehh_pergene_tbl <- ds2_pergene_top$xpehh %>%
  filter(!grepl(exclude_gene_pattern, GENE)) %>%
  mutate(PBS = round(SCORE, 4), `-log10(p)` = round(LOG_PVALUE, 4)) %>%
  select(-all_of(c("SCORE", "PVALUE", "LOG_PVALUE", "GENE_ID")))

ds2_xpehh_pergene_tbl$GENE %>% strsplit(",") %>% unlist %>% unique %>% length

ds2_xpehh_pergene_kbl <- kbl(ds2_xpehh_pergene_tbl, booktabs = T, format = "latex", linesep = "")

save_kable(ds2_xpehh_pergene_kbl, "tables/tex/ds2_xpehh_pergene.tex")

wilcox.test(ds2_pergene$xpehh$SCORE, mu = min(ds2_pergene_top$xpehh$SCORE))

######  2.3 XLS  ######

xpehh_data <- list (
  ds1_top0.01_snp = ds1_xpehh_persnp_tbl,
  ds1_top0.01_gene = ds1_xpehh_pergene_tbl,
  ds2_top0.01_snp = ds2_xpehh_persnp_tbl,
  ds2_top0.01_gene = ds2_xpehh_pergene_tbl
)

WriteXLS::WriteXLS(xpehh_data, "tables/xlsx/xpehh.xls", row.names = F, AdjWidth = T, BoldHeaderRow = T)

#### 2 FIGURES ####
####################################################################################################!

######  2.1 DATASET S1  ######

## -- XP-EHH Per SNP -- ##

ds1_xpehh_persnp_plot <- manhattan_plot (
  ds1_persnp$xpehh, annot = "GENE", score = "LOG_PVALUE", 
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"), plot.title.size = 18,
  plot.title = "Dataset S1 - XPEHH por SNP",
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds1_xpehh_persnp_plot <- annotate_plot (
  ds1_xpehh_persnp_plot, threshold = .999, y.upper.limit = 6,
  nudge_y = 0.5, size = 4.2, alpha = .8, segment.alpha = .5
)

ggsave("figures/png/ds1_xpehh_persnp.png", ds1_xpehh_persnp_plot, height = 6, width = 11)
ggsave("figures/pdf/ds1_xpehh_persnp.pdf", ds1_xpehh_persnp_plot, height = 6, width = 11)

## -- XP-EHH Per Gene -- ##

ds1_xpehh_pergene_plot <- manhattan_plot (
  ds1_pergene$xpehh, snp = NULL, annot = "GENE", score = "LOG_PVALUE",
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"),
  plot.title = "Dataset S1 - XP-EHH por Gene", plot.title.size = 18,
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds1_xpehh_pergene_plot <- annotate_plot (
  ds1_xpehh_pergene_plot, threshold = .999, y.upper.limit = 4.8,
  nudge_y = 0.1, size=4.2, alpha=.8, segment.alpha = .5
)

ggsave("figures/png/ds1_xpehh_pergene.png", ds1_xpehh_pergene_plot, height = 6, width = 11)
ggsave("figures/pdf/ds1_xpehh_pergene.pdf", ds1_xpehh_pergene_plot, height = 6, width = 11)

## -- XP-EHH SNP & Gene -- ##

ds1_xpehh_title <- ggdraw() + draw_label("Dataset S1 - XP-EHH", fontface='bold', size = 21)

ds1_xpehh_grid <- plot_grid (
  ds1_xpehh_persnp_plot + labs(x='cromossomo', title=""), ds1_xpehh_pergene_plot + labs(x='cromossomo', title=''), 
  labels = "AUTO", label_size = 18, nrow = 2, ncol = 1, rel_heights = c(.55, .45)
)

ds1_xpehh_grid <- plot_grid(ds1_xpehh_title, ds1_xpehh_grid, ncol=1, rel_heights=c(0.08, 1))

ggsave("figures/png/ds1_xpehh_both.png", ds1_xpehh_grid, height = 11, width = 11, dpi=300)
ggsave("figures/pdf/ds1_xpehh_both.pdf", ds1_xpehh_grid, height = 11, width = 11, dpi=300)

######  2.2 DATASET S2  ######

## -- Plot Per SNP -- ##

ds2_xpehh_persnp_plot <- manhattan_plot (
  ds2_persnp$xpehh, annot = "GENE", score = "LOG_PVALUE", 
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"),
  plot.title.size = 18, plot.title = "Dataset S2 - XPEHH por SNP",
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds2_xpehh_persnp_plot <- annotate_plot (
  ds2_xpehh_persnp_plot, threshold = .999,
  nudge_y = 0.5, size = 4.2, alpha = .8, segment.alpha = .5
)

ggsave("figures/png/ds2_xpehh_persnp.png", ds2_xpehh_persnp_plot, height = 6, width = 11)
ggsave("figures/pdf/ds2_xpehh_persnp.pdf", ds2_xpehh_persnp_plot, height = 6, width = 11)

## -- Plot Per Gene -- ##

ds2_xpehh_pergene_plot <- manhattan_plot (
  ds2_pergene$xpehh, snp = NULL, annot = "GENE", score = "LOG_PVALUE",
  point.alpha = 1, point.color = c("#0086AB", "#9A9A9A"),
  plot.title = "Dataset S2 - XP-EHH por Gene", plot.title.size = 18,
  threshold_lines = c(.995, .999), threshold_colors = c("orange", "red"),
  x.axis.title = "cromossomo"
)

ds2_xpehh_pergene_plot <- annotate_plot (
  ds2_xpehh_pergene_plot, threshold = .999, 
  nudge_y = 0.3, size=4.2, alpha=.8, segment.alpha = .5
)

ggsave("figures/png/ds2_xpehh_pergene.png", ds2_xpehh_pergene_plot, height = 6, width = 11)
ggsave("figures/pdf/ds2_xpehh_pergene.pdf", ds2_xpehh_pergene_plot, height = 6, width = 11)

## -- XP-EHH SNP & Gene -- ##

# Dataset S2 - XP-EHH - snp vs gene
ds2_xpehh_title <- ggdraw() + draw_label("Dataset S2 - XP-EHH", fontface='bold', size = 21)

ds2_xpehh_grid <- plot_grid (
  ds2_xpehh_persnp_plot + labs(x='cromossomo', title=""), ds2_xpehh_pergene_plot + labs(x='cromossomo', title=''), 
  labels = "AUTO", label_size = 18, nrow = 2, ncol = 1, rel_heights = c(.55, .45)
)

ds2_xpehh_grid <- plot_grid(ds2_xpehh_title, ds2_xpehh_grid, ncol=1, rel_heights=c(0.08, 1))

ggsave("figures/png/ds2_xpehh_both.png", ds2_xpehh_grid, height = 11, width = 11, dpi=300)
ggsave("figures/pdf/ds2_xpehh_both.pdf", ds2_xpehh_grid, height = 11, width = 11, dpi=300)

