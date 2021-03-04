# Make density dataframe for Dataset S1

# Set PVALUE & LOG_PVALUE
for(method in names(ds1_persnp)) {
  col <- toupper(method)
  if(method=="pbs_windowed") col <- "PBS"
  
  ds1_persnp[[method]][, PVALUE := 1 - frank(SCORE) / (.N+1)]
  ds1_persnp[[method]][, LOG_PVALUE := -log10(PVALUE)]
  
  ds2_persnp[[method]][, PVALUE := 1 - frank(SCORE) / (.N+1)]
  ds2_persnp[[method]][, LOG_PVALUE := -log10(PVALUE)]
}

for(method in names(ds1_pergene)) {
  col <- toupper(method)
  if(method=="pbs_windowed") col <- "PBS"
  
  ds1_pergene[[method]][, PVALUE := 1 - frank(SCORE) / (.N+1)]
  ds1_pergene[[method]][, LOG_PVALUE := -log10(PVALUE)]
  
  ds2_pergene[[method]][, PVALUE := 1 - frank(SCORE) / (.N+1)]
  ds2_pergene[[method]][, LOG_PVALUE := -log10(PVALUE)]
}

## -- Dataset S1 -- ##

# Get tables
ds1_persnp_density_dt <- rbindlist(list(
  ds1_persnp$pbs_windowed[, .(SCORE, Method="PBS")],
  ds1_persnp$xpehh[, .(SCORE, Method="XP-EHH")],
  ds1_persnp$ihs[, .(SCORE, Method="iHS")]
))

ds1_pergene_density_dt <- rbindlist(list(
  ds1_pergene$pbs_windowed[, .(SCORE, Method="PBS")],
  ds1_pergene$xpehh[, .(SCORE, Method="XP-EHH")],
  ds1_pergene$ihs[, .(SCORE, Method="iHS")]
))

ds1_persnp_density_dt[, Method := factor(Method, unique(Method))]
ds1_pergene_density_dt[, Method := factor(Method, unique(Method))]

# Get quantiles
ds1_persnp_density_ref <- cbind (
  ds1_persnp_density_dt[, .(top01 = quantile(SCORE, .999)), by = Method],
  ds1_persnp_density_dt[, .(top05 = quantile(SCORE, .995)), by = Method][, .(top05)]
)

ds1_pergene_density_ref <- cbind (
  ds1_pergene_density_dt[, .(top01 = quantile(SCORE, .999)), by = Method],
  ds1_pergene_density_dt[, .(top05 = quantile(SCORE, .995)), by = Method][, .(top05)]
)

# ds1_persnp_density_dt <- ds1_persnp_density_dt[SCORE>=0]
# ds1_pergene_density_dt <- ds1_pergene_density_dt[SCORE>=0]

# Make density plot ─ Dataset S1
ds1_density_dt <- list(SNP = ds1_persnp_density_dt, Gene = ds1_pergene_density_dt)
ds1_density_ref <- list(SNP = ds1_persnp_density_ref, Gene = ds1_pergene_density_ref)

ds1_density_plt <- sapply(1:2, function(i) {
  dt <- ds1_density_dt[[i]]
  method <- names(ds1_density_dt)[i]
  method <- ifelse(method=="Gene", "gene", method)
  
  ggplot(dt, aes(x=SCORE)) +
    geom_density() +
    facet_wrap(~Method, scales = "free") +
    geom_vline(data = ds1_density_ref[[i]], aes(xintercept = top05), linetype = "dashed", color = "orange") +
    geom_vline(data = ds1_density_ref[[i]], aes(xintercept = top01), linetype = "dashed", color = "red") +
    labs(x = "valores", y = "densidade", title = glue::glue("Dataset S1 - Distribuição dos resultados por {method}")) +
    theme_bw() +
    theme (
      plot.title = element_text(size = 16, face = "bold", margin = margin(t = 10, b = 15), hjust = 0.5),
      strip.text = element_text(size = 14, face = "bold"),
      text = element_text(size = 12)
    )}, simplify = F) %>%
  set_names(c("persnp", "pergene"))

ggsave("figures/png/ds1_persnp_density.png", ds1_density_plt$persnp, width = 8, height = 4)
ggsave("figures/pdf/ds1_persnp_density.pdf", ds1_density_plt$persnp, width = 8, height = 4)
ggsave("figures/png/ds1_pergene_density.png", ds1_density_plt$pergene, width = 8, height = 4)
ggsave("figures/pdf/ds1_pergene_density.pdf", ds1_density_plt$pergene, width = 8, height = 4)

## -- Dataset S2 -- ##

# Get tables
ds2_persnp_density_dt <- rbindlist(list(
  ds2_persnp$pbs_windowed[, .(SCORE, Method="PBS")],
  ds2_persnp$xpehh[, .(SCORE, Method="XP-EHH")],
  ds2_persnp$ihs[, .(SCORE, Method="iHS")]
))

ds2_pergene_density_dt <- rbindlist(list(
  ds2_pergene$pbs_windowed[, .(SCORE, Method="PBS")],
  ds2_pergene$xpehh[, .(SCORE, Method="XP-EHH")],
  ds2_pergene$ihs[, .(SCORE, Method="iHS")]
))

ds2_persnp_density_dt[, Method := factor(Method, unique(Method))]
ds2_pergene_density_dt[, Method := factor(Method, unique(Method))]

# Get quantiles
ds2_persnp_density_ref <- cbind (
  ds2_persnp_density_dt[, .(top01 = quantile(SCORE, .999)), by = Method],
  ds2_persnp_density_dt[, .(top05 = quantile(SCORE, .995)), by = Method][, .(top05)]
)

ds2_pergene_density_ref <- cbind (
  ds2_pergene_density_dt[, .(top01 = quantile(SCORE, .999)), by = Method],
  ds2_pergene_density_dt[, .(top05 = quantile(SCORE, .995)), by = Method][, .(top05)]
)

# ds2_persnp_density_dt <- ds2_persnp_density_dt[SCORE>=0]
# ds2_pergene_density_dt <- ds2_pergene_density_dt[SCORE>=0]

# Make density plot ─ Dataset S2
ds2_density_dt <- list(SNP = ds2_persnp_density_dt, Gene = ds2_pergene_density_dt)
ds2_density_ref <- list(SNP = ds2_persnp_density_ref, Gene = ds2_pergene_density_ref)

ds2_density_plt <- sapply(1:2, function(i) {
  dt <- ds2_density_dt[[i]]
  method <- names(ds2_density_dt)[i]
  method <- ifelse(method=="Gene", "gene", method)
  
  ggplot(dt, aes(x=SCORE)) +
    geom_density() +
    facet_wrap(~Method, scales = "free") +
    geom_vline(data = ds2_density_ref[[i]], aes(xintercept = top05), linetype = "dashed", color = "orange") +
    geom_vline(data = ds2_density_ref[[i]], aes(xintercept = top01), linetype = "dashed", color = "red") +
    labs(x = "valores", y = "densidade", title = glue::glue("Dataset S2 - Distribuição dos resultados por {method}")) +
    theme_bw() +
    theme (
      plot.title = element_text(size = 16, face = "bold", margin = margin(t = 10, b = 15), hjust = 0.5),
      strip.text = element_text(size = 14, face = "bold"),
      text = element_text(size = 12)
    )}, simplify = F) %>%
  set_names(c("persnp", "pergene"))

ggsave("figures/png/ds2_persnp_density.png", ds2_density_plt$persnp, width = 8, height = 4)
ggsave("figures/pdf/ds2_persnp_density.pdf", ds2_density_plt$persnp, width = 8, height = 4)
ggsave("figures/png/ds2_pergene_density.png", ds2_density_plt$pergene, width = 8, height = 4)
ggsave("figures/pdf/ds2_pergene_density.pdf", ds2_density_plt$pergene, width = 8, height = 4)

# Make Grid Density Plot ─ Dataset S1

t <- theme (
  plot.title = element_text(size = 14, face = "bold.italic", margin = margin(t = 10, b = 15), hjust = 0.5),
  strip.text = element_text(size = 10, face = "bold"),
  text = element_text(size = 12)
)

ds1_density_title <- ggdraw() + draw_label("       Dataset S1", fontface='bold', size = 18)

ds1_density_plt_grid <- plot_grid (
  ds1_density_plt$persnp + labs(title = "Distribuição dos resultados por SNP", x = "") + t,
  ds1_density_plt$pergene + labs(title = "Distribuição dos resultados por gene") + t,
  nrow = 2, labels = "AUTO", label_y = 0.95
  )

ds1_density_plt_grid <- plot_grid (
  ds1_density_title,
  ds1_density_plt_grid, nrow = 2, rel_heights = c(0.05, 1)
)

ggsave("figures/png/ds1_density.png", ds1_density_plt_grid, width = 8, height = 6)
ggsave("figures/pdf/ds1_density.pdf", ds1_density_plt_grid, width = 8, height = 6)

ds2_density_title <- ggdraw() + draw_label("       Dataset S2", fontface='bold', size = 18)

ds2_density_plt_grid <- plot_grid (
  ds2_density_plt$persnp + labs(title = "Distribuição dos resultados por SNP", x = "") + t,
  ds2_density_plt$pergene + labs(title = "Distribuição dos resultados por gene") + t,
  nrow = 2, labels = "AUTO", label_y = 0.95
)

ds2_density_plt_grid <- plot_grid (
  ds2_density_title,
  ds2_density_plt_grid, nrow = 2, rel_heights = c(0.05, 1)
)

ggsave("figures/png/ds2_density.png", ds2_density_plt_grid, width = 8, height = 6)
ggsave("figures/pdf/ds2_density.pdf", ds2_density_plt_grid, width = 8, height = 6)
