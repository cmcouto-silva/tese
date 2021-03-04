#### SETTTINGS ####
####################################################################################################!

get_n <- function(col) {
  strsplit(col, ",") %>%
    unlist %>%
    table %>%
    as.data.table() %>%
    set_names(c("Tissue", "N")) %>%
    arrange(-N, Tissue) %>%
    mutate(Tissue = factor(Tissue, unique(Tissue)))
}

#### MAIN ####
####################################################################################################!

datasets <- c("dataset_S1", "dataset_S2")
methods <- c("pbsw", "xpehh")

eqtl_plots <- sapply(datasets, function(dataset) {
  ds <- ifelse(dataset=="dataset_S1", "ds1", "ds2")
  ds_title <- stringr::str_to_title(gsub("_", " ", dataset))
  
  sapply(methods, function(method) {
    
    dt <- fread(glue("{dataset}/gtex_annot/eqtl_persnp_{method}.csv"))
    method_title <- ifelse(method=="pbsw", "PBS", "XP-EHH")
    
    dt_up <- dt[, get_n(UP)]
    dt_down <- dt[, get_n(DOWN)]
    
    dt <- merge(dt_up, dt_down, by = "Tissue", all = T)
    setnames(dt, c("N.x", "N.y"), c("UP", "DOWN"))
    
    dt[, c("UP", "DOWN") := lapply(.SD, function(x) as.numeric(nafill(x, fill=0))), .SDcols = c("UP", "DOWN")]
    dt[, SCORE := UP - DOWN]
    dt <- dt[order(-SCORE)]
    dt[, Tissue := factor(Tissue, unique(Tissue))]
    
    dt_melt <- melt(dt, "Tissue", c("UP", "DOWN"))
    dt_melt[, Tissue := factor(Tissue, levels(dt$Tissue))]
    dt_melt[value==0, value := 0.05]
    
    eqtl_score_plt <- ggplot(dt[SCORE!=0], aes(x=Tissue, y=SCORE, fill=SCORE)) +
      scale_fill_gradientn(colours = c("darkred", "black", "darkgreen"),
                           values = scales::rescale(c(-8, 0, 2))) +
      geom_col() +
      geom_hline(yintercept = 0) +
      theme_bw() +
      theme (
        plot.title = element_text(size=18, face = "bold", margin = margin(t = 10, b = 15)),
        axis.title.y.left = element_text(size=11, margin = margin(l = 10, r = 5)),
        axis.text.x.bottom = element_text(angle = 60, hjust = 1), 
        axis.line.x.bottom = element_blank(),
        text = element_text(size = 12),
        legend.position = "none"
      ) +
      labs(title=glue("{ds_title} - {method_title} eQTL Score"), x="")
    
    eqtl_plt <- ggplot(dt_melt, aes(x=Tissue, y=value, fill=variable)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("darkgreen", "darkred")) +
      theme_bw() +
      theme (
        legend.position = "top",
        legend.justification = "left",
        text = element_text(size=12),
        plot.title = element_text(size=18, face = "bold", hjust = .5, margin = margin(t = 10, b = 15)),
        axis.title.y.left = element_text(size=11, margin = margin(l = 10, r = 5)),
        axis.text.x.bottom = element_text(angle = 60, hjust = 1), 
        panel.grid = element_blank()
      ) +
      labs(title=glue("{ds_title} - {method_title} eQTL"), y="count", x="", fill="Expression") +
      geom_vline(xintercept = seq(0.5, length(unique(dt_melt$Tissue)), by = 1), color="gray", size=.5, alpha=.5)
    
    plots <- list(eqtl=eqtl_plt, eqtl_score=eqtl_score_plt)
    
    return(plots)
    
  }, simplify = F)
}, simplify = F)

# Datset S1 ─ PBS
eqtl_plots$dataset_S1$pbsw$eqtl <- eqtl_plots$dataset_S1$pbsw$eqtl + theme (
  axis.title.y.left = element_text(margin = margin(l = 10, r = 10))
)

ggsave("figures/png/ds1_pbs_eqtl_score.png", eqtl_plots$dataset_S1$pbsw$eqtl_score, width = 9, height = 6)
ggsave("figures/png/ds1_pbs_eqtl.png", eqtl_plots$dataset_S1$pbsw$eqtl, width = 11, height = 6)

eqtl_plots$dataset_S1$pbsw$eqtl_score <- eqtl_plots$dataset_S1$pbsw$eqtl_score + theme (
  axis.title.y.left = element_text(margin = margin(l = 10, r = 10))
)

ggsave("figures/pdf/ds1_pbs_eqtl_score.pdf", eqtl_plots$dataset_S1$pbsw$eqtl_score, width = 9, height = 6)
ggsave("figures/pdf/ds1_pbs_eqtl.pdf", eqtl_plots$dataset_S1$pbsw$eqtl, width = 11, height = 6)

# Datset S1 ─ XP-EHH
eqtl_plots$dataset_S1$xpehh$eqtl <- eqtl_plots$dataset_S1$xpehh$eqtl + theme (
  axis.title.y.left = element_text(margin = margin(l = 5, r = 5))
)

ggsave("figures/png/ds1_xpehh_eqtl_score.png", eqtl_plots$dataset_S1$xpehh$eqtl_score, width = 10, height = 6)
ggsave("figures/png/ds1_xpehh_eqtl.png", eqtl_plots$dataset_S1$xpehh$eqtl, width = 11, height = 6)

ggsave("figures/pdf/ds1_xpehh_eqtl_score.pdf", eqtl_plots$dataset_S1$xpehh$eqtl_score, width = 10, height = 6)
ggsave("figures/pdf/ds1_xpehh_eqtl.pdf", eqtl_plots$dataset_S1$xpehh$eqtl, width = 11, height = 6)

# Datset S2 ─ PBS
eqtl_plots$dataset_S2$pbsw$eqtl <- eqtl_plots$dataset_S2$pbsw$eqtl + theme (
  axis.title.y.left = element_text(margin = margin(l = 10, r = 10))
)

eqtl_plots$dataset_S2$pbsw$eqtl_score <- eqtl_plots$dataset_S2$pbsw$eqtl_score + theme (
  axis.title.y.left = element_text(margin = margin(l = 10, r = 10))
)

ggsave("figures/png/ds2_pbs_eqtl_score.png", eqtl_plots$dataset_S2$pbsw$eqtl_score, width = 9, height = 6)
ggsave("figures/png/ds2_pbs_eqtl.png", eqtl_plots$dataset_S2$pbsw$eqtl, width = 11, height = 6)

ggsave("figures/pdf/ds2_pbs_eqtl_score.pdf", eqtl_plots$dataset_S2$pbsw$eqtl_score, width = 9, height = 6)
ggsave("figures/pdf/ds2_pbs_eqtl.pdf", eqtl_plots$dataset_S2$pbsw$eqtl, width = 11, height = 6)

# Datset S2 ─ XP-EHH
eqtl_plots$dataset_S2$xpehh$eqtl <- eqtl_plots$dataset_S2$xpehh$eqtl + theme (
  axis.title.y.left = element_text(margin = margin(l = 5, r = 5))
)

ggsave("figures/png/ds2_xpehh_eqtl_score.png", eqtl_plots$dataset_S2$xpehh$eqtl_score, width = 10, height = 6)
ggsave("figures/png/ds2_xpehh_eqtl.png", eqtl_plots$dataset_S2$xpehh$eqtl, width = 11, height = 6)

ggsave("figures/pdf/ds2_xpehh_eqtl_score.pdf", eqtl_plots$dataset_S2$xpehh$eqtl_score, width = 10, height = 6)
ggsave("figures/pdf/ds2_xpehh_eqtl.pdf", eqtl_plots$dataset_S2$xpehh$eqtl, width = 11, height = 6)

## -- Comparação Total -- ##

# eQTL score

eqtl_score <- rbindlist(list(
  eqtl_plots$dataset_S1$pbsw$eqtl_score$data,
  eqtl_plots$dataset_S2$pbsw$eqtl_score$data,
  eqtl_plots$dataset_S1$xpehh$eqtl_score$data,
  eqtl_plots$dataset_S2$xpehh$eqtl_score$data
))

eqtl_score[, SCORE := UP + DOWN]
eqtl_score <- eqtl_score[, .(SCORE=sum(SCORE)), Tissue][order(SCORE, -Tissue)]
eqtl_score[, Tissue := factor(Tissue, unique(Tissue))]

eqtl_score_plt <- ggplot(eqtl_score[SCORE!=0], aes(y=Tissue, x=SCORE, fill=SCORE)) +
  scale_fill_gradientn(colours = c("gray", "#1FA075")) +
  geom_col() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme (
    plot.title = element_text(size=18, face = "bold", margin = margin(t = 10, b = 20)),
    axis.text.y.left = element_text(size=12),
    text = element_text(size = 12)
  ) +
  labs(title="eQTL entre métodos e datasets", x="score", y="", fill="Score")

ggsave("figures/png/overall_eqtl_score.png", eqtl_score_plt, width = 10, height = 10)
ggsave("figures/pdf/overall_eqtl_score.pdf", eqtl_score_plt, width = 10, height = 10)

# eQTL count
eqtl <- rbindlist(list(
  eqtl_plots$dataset_S1$pbsw$eqtl$data,
  eqtl_plots$dataset_S2$pbsw$eqtl$data,
  eqtl_plots$dataset_S1$xpehh$eqtl$data,
  eqtl_plots$dataset_S2$xpehh$eqtl$data  
))

eqtl[value==0.05, value := 0]
eqtl <- eqtl[, .(value=sum(value)), c("Tissue", "variable")]

eqtl[, Tissue := factor(Tissue, sort(unique(as.character(Tissue))))]
eqtl[, sort(unique(Tissue))]

eqtl_count_plt <- ggplot(eqtl, aes(x=Tissue, y=value, fill=variable)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#2385A8", "#d1495b")) +
  theme_bw() +
  theme (
    legend.position = "top",
    legend.justification = "left",
    text = element_text(size=12),
    plot.title = element_text(size=18, face = "bold", hjust = .5, margin = margin(t = 15, b = 10)),
    axis.title.y.left = element_text(size=11, margin = margin(l = 10, r = 10)),
    axis.text.x.bottom = element_text(angle = 60, hjust = 1), 
    panel.grid = element_blank()
  ) +
  labs(title="eQTL significativos entre métodos e datasets", x="", y="count", fill="Expression") +
  geom_vline(xintercept = seq(0.5, length(unique(eqtl$Tissue)), by = 1), color="gray", size=.5, alpha=.5)

ggsave("figures/png/overall_eqtl_count.png", eqtl_count_plt, width = 12, height = 7)
ggsave("figures/pdf/overall_eqtl_count.pdf", eqtl_count_plt, width = 12, height = 7)

eqtl_split <- copy(eqtl)
eqtl_split[variable=="DOWN", value := value*-1]
eqtl_split[, Tissue := factor(Tissue, sort(unique(Tissue), decreasing = T))]

eqtl_split_plt <- ggplot(eqtl_split, aes(x=value, y=Tissue, fill=variable)) +
  geom_col() +
  scale_fill_manual(values = c("#d1495b", "#2385A8")) +
  theme_bw() +
  theme (
    legend.position = "top",
    legend.justification = "left",
    text = element_text(size=12),
    plot.title = element_text(size=18, face = "bold", hjust = 0, margin = margin(t = 15, b = 10)),
    axis.text.y.left = element_text(size=12),
    axis.text.x.bottom = element_text(hjust = 1)
  ) +
  labs(title="eQTL significativos entre métodos e datasets", x="score", y="", fill="Expression")

ggsave("figures/png/overall_eqtl_count_split.png", eqtl_split_plt, width = 10, height = 11)
ggsave("figures/pdf/overall_eqtl_count_split.pdf", eqtl_split_plt, width = 10, height = 11)

## -- Spreadsheet -- ##

eqtl_data <- list (
  ds1_pbs_eqtl = eqtl_plots$dataset_S1$pbsw$eqtl_score$data,
  ds1_xpehh_eqtl = eqtl_plots$dataset_S1$xpehh$eqtl_score$data,
  ds2_pbs_eqtl = eqtl_plots$dataset_S2$pbsw$eqtl_score$data,
  ds2_xpehh_eqtl = eqtl_plots$dataset_S2$xpehh$eqtl_score$data
)

WriteXLS::WriteXLS(eqtl_data, "tables/xlsx/eqtl_score.xls", row.names = F, AdjWidth = T, BoldHeaderRow = T)

