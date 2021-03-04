## -- Dataset S1 vs Dataset S2 -- ##

# PBS per SNP
pbs_persnp_intersection.999 <- list (
  `Dataset S1` = get_top_intersection(ds1_persnp$pbs_windowed, .999, "snp"),
  `Dataset S2` = get_top_intersection(ds2_persnp$pbs_windowed, .999, "snp")
)

pbs_persnp_intersection.995 <- list (
  `Dataset S1` = get_top_intersection(ds1_persnp$pbs_windowed, .995, "snp"),
  `Dataset S2` = get_top_intersection(ds2_persnp$pbs_windowed, .995, "snp")
)

pbs_persnp_intersection.999_plt <- ggven_plot(pbs_persnp_intersection.999, n_cat = 2, internal_annot_size = 5, 
                                              external_annot_size = 6, external_annot_fontface = "", size=.75)

pbs_persnp_intersection.995_plt <- ggven_plot(pbs_persnp_intersection.995, n_cat = 2, internal_annot_size = 5,
                                              external_annot_size = 6, external_annot_fontface = "", size=.75)

ggsave("figures/png/pbs_persnp_intersection.999.png", pbs_persnp_intersection.999_plt, height = 3, width = 7)
ggsave("figures/pdf/pbs_persnp_intersection.999.pdf", pbs_persnp_intersection.999_plt, height = 3, width = 7)

ggsave("figures/png/pbs_persnp_intersection.995.png", pbs_persnp_intersection.995_plt, height = 3, width = 7)
ggsave("figures/pdf/pbs_persnp_intersection.995.pdf", pbs_persnp_intersection.995_plt, height = 3, width = 7)

# PBS per Gene
pbs_pergene_intersection.999 <- list (
  `Dataset S1` = get_top_intersection(ds1_pergene$pbs_windowed, .999, "gene"),
  `Dataset S2` = get_top_intersection(ds2_pergene$pbs_windowed, .999, "gene")
)

pbs_pergene_intersection.995 <- list (
  `Dataset S1` = get_top_intersection(ds1_pergene$pbs_windowed, .995, "gene"),
  `Dataset S2` = get_top_intersection(ds2_pergene$pbs_windowed, .995, "gene")
)

pbs_pergene_intersection.999_plt <- ggven_plot(pbs_pergene_intersection.999, n_cat = 2, internal_annot_size = 5,
                                               external_annot_size = 6, external_annot_fontface = "", size=.75)

pbs_pergene_intersection.995_plt <- ggven_plot(pbs_pergene_intersection.995, n_cat = 2, internal_annot_size = 5, 
                                               external_annot_fontface = "", external_annot_size = 6, size=.75)

ggsave("figures/png/pbs_pergene_intersection.999.png", pbs_pergene_intersection.999_plt, height = 3, width = 7)
ggsave("figures/pdf/pbs_pergene_intersection.999.pdf", pbs_pergene_intersection.999_plt, height = 3, width = 7)

ggsave("figures/png/pbs_pergene_intersection.995.png", pbs_pergene_intersection.995_plt, height = 3, width = 7)
ggsave("figures/pdf/pbs_pergene_intersection.995.pdf", pbs_pergene_intersection.995_plt, height = 3, width = 7)

# PBS Both
pbs_intersection_title <- ggdraw() + draw_label("PBS - Intersecção de genes entre os conjuntos de dados", fontface='bold', size = 21)

pbs_intersection_subtitle01 <- ggdraw() + draw_label("p-valor < 0,01", fontface='bold', size = 18)
pbs_intersection_subtitle02 <- ggdraw() + draw_label("p-valor < 0,05", fontface='bold', size = 18)
pbs_intersection_subtitle <- plot_grid(pbs_intersection_subtitle01, pbs_intersection_subtitle02)

pbs_persnp_intersection <- plot_grid (
  pbs_persnp_intersection.999_plt, pbs_persnp_intersection.995_plt,
  label_x = 0.25, label_y = 0.95, rel_heights = 0.5, label_size = 16
  )

pbs_pergene_intersection <- plot_grid (
  pbs_pergene_intersection.999_plt, pbs_pergene_intersection.995_plt, 
  label_x = 0.25, label_y = 0.95, rel_heights = 0.5, label_size = 16
)

pbs_intersection <- plot_grid(pbs_persnp_intersection, pbs_pergene_intersection, labels = "AUTO", nrow = 2, label_size = 18)
pbs_intersection <- plot_grid(pbs_intersection_title, pbs_intersection_subtitle, pbs_intersection, ncol=1, rel_heights=c(0.1, 0.1, 1))

ggsave("figures/png/pbs_intersection.png", pbs_intersection, height = 7.5, width = 10)
ggsave("figures/pdf/pbs_intersection.pdf", pbs_intersection, height = 7.5, width = 10)


## -- Methods' Intersection in dataset S1 -- ##

get_top_intersection <- function(dt, cutoff, mode = "snp", pbsw_correction = TRUE) {
  dt <- get_top(dt, cutoff, mode = mode, pbsw_correction = pbsw_correction)
  strsplit(dt$GENE, ",") %>% unlist %>% unique %>% sort
}

ds1_persnp_top_intersection.999 <- lapply(ds1_persnp[-1], get_top_intersection, .999, pbsw_correction = T)
ds1_persnp_top_intersection.995 <- lapply(ds1_persnp[-1], get_top_intersection, .995, pbsw_correction = T)

ds1_persnp_top_intersection.999_plt <- ggven_plot (internal_label = "count",
  ds1_persnp_top_intersection.999, internal_annot_size = 5, external_annot_size = 5,
  size=.75, external_label = c("PBS", "iHS", "XP-EHH"), external_annot_fontface = ""
  )

ds1_persnp_top_intersection.995_plt <- ggven_plot (internal_label = "count",
  ds1_persnp_top_intersection.995, internal_annot_size = 5, external_annot_size = 5,
  size=.75, external_label = c("PBS", "iHS", "XP-EHH"), external_annot_fontface = ""
)

ds1_pergene_top_intersection.999 <- lapply(ds1_pergene[-1], get_top_intersection, .999, mode = "gene", pbsw_correction = T)
ds1_pergene_top_intersection.995 <- lapply(ds1_pergene[-1], get_top_intersection, .995, mode = "gene", pbsw_correction = T)

ds1_pergene_top_intersection.999_plt <- ggven_plot (internal_label = "count",
  ds1_pergene_top_intersection.999, internal_annot_size = 5, external_annot_size = 5,
  size=.75, external_label = c("PBS", "iHS", "XP-EHH"), external_annot_fontface = ""
)

ds1_pergene_top_intersection.995_plt <- ggven_plot (internal_label = "count",
  ds1_pergene_top_intersection.995, internal_annot_size = 5, external_annot_size = 5,
  size=.75, external_label = c("PBS", "iHS", "XP-EHH"), external_annot_fontface = ""
)

ds1_intersection_plt_list <- list (
  ds1_persnp_top_intersection.999_plt,
  ds1_persnp_top_intersection.995_plt,
  ds1_pergene_top_intersection.999_plt,
  ds1_pergene_top_intersection.995_plt
)

ds1_intersection_title01 <- ggdraw() +
  draw_label("Dataset S1", fontface='bold', size = 20)

ds1_intersection_title02 <- ggdraw() +
  draw_label("Intersecção de genes candidatos entre métodos", fontface='bold', size = 18)

ds1_intersection_title <- plot_grid(ds1_intersection_title01, ds1_intersection_title02, nrow = 2, rel_heights = c(1.2, 0.8))

ds1_intersection_subtitle01 <- ggdraw() + draw_label("Top 0,1%", fontface='bold', size = 16)
ds1_intersection_subtitle02 <- ggdraw() + draw_label("Top 0,5%", fontface='bold', size = 16)
ds1_intersection_subtitle <- plot_grid(ds1_intersection_subtitle01, ds1_intersection_subtitle02)

ds1_intersection_persnp_plt <- plot_grid(ds1_persnp_top_intersection.999_plt, ds1_persnp_top_intersection.995_plt)
ds1_intersection_pergene_plt <- plot_grid(ds1_pergene_top_intersection.999_plt, ds1_pergene_top_intersection.995_plt)

ds1_intersection_plt <- plot_grid(ds1_intersection_persnp_plt, ds1_intersection_pergene_plt, nrow = 2, labels = "AUTO")
ds1_intersection_plt <- plot_grid(ds1_intersection_title, ds1_intersection_subtitle, ds1_intersection_plt, ncol=1, rel_heights=c(0.13, 0.08, .9))

ggsave("figures/png/ds1_intersection.png", ds1_intersection_plt, height = 7, width = 9)
ggsave("figures/pdf/ds1_intersection.pdf", ds1_intersection_plt, height = 7, width = 9)

## -- Methods' Intersection in dataset S2 -- ##

get_top_intersection <- function(dt, cutoff, mode = "snp", pbsw_correction = TRUE) {
  dt <- get_top(dt, cutoff, mode = mode, pbsw_correction = pbsw_correction)
  strsplit(dt$GENE, ",") %>% unlist %>% unique %>% sort
}

ds2_persnp_top_intersection.999 <- lapply(ds2_persnp[-1], get_top_intersection, .999, pbsw_correction = T)
ds2_persnp_top_intersection.995 <- lapply(ds2_persnp[-1], get_top_intersection, .995, pbsw_correction = T)

ds2_persnp_top_intersection.999_plt <- ggven_plot (internal_label = "count",
                                                   ds2_persnp_top_intersection.999, internal_annot_size = 5, external_annot_size = 5,
                                                   size=.75, external_label = c("PBS", "iHS", "XP-EHH"), external_annot_fontface = ""
)

ds2_persnp_top_intersection.995_plt <- ggven_plot (internal_label = "count",
                                                   ds2_persnp_top_intersection.995, internal_annot_size = 5, external_annot_size = 5,
                                                   size=.75, external_label = c("PBS", "iHS", "XP-EHH"), external_annot_fontface = ""
)

ds2_pergene_top_intersection.999 <- lapply(ds2_pergene[-1], get_top_intersection, .999, mode = "gene", pbsw_correction = T)
ds2_pergene_top_intersection.995 <- lapply(ds2_pergene[-1], get_top_intersection, .995, mode = "gene", pbsw_correction = T)

ds2_pergene_top_intersection.999_plt <- ggven_plot (internal_label = "count",
                                                    ds2_pergene_top_intersection.999, internal_annot_size = 5, external_annot_size = 5,
                                                    size=.75, external_label = c("PBS", "iHS", "XP-EHH"), external_annot_fontface = ""
)

ds2_pergene_top_intersection.995_plt <- ggven_plot (internal_label = "count",
                                                    ds2_pergene_top_intersection.995, internal_annot_size = 5, external_annot_size = 5,
                                                    size=.75, external_label = c("PBS", "iHS", "XP-EHH"), external_annot_fontface = ""
)

ds2_intersection_plt_list <- list (
  ds2_persnp_top_intersection.999_plt,
  ds2_persnp_top_intersection.995_plt,
  ds2_pergene_top_intersection.999_plt,
  ds2_pergene_top_intersection.995_plt
)

ds2_intersection_title01 <- ggdraw() +
  draw_label("Dataset S2", fontface='bold', size = 20)

ds2_intersection_title02 <- ggdraw() +
  draw_label("Intersecção de genes candidatos entre métodos", fontface='bold', size = 18)

ds2_intersection_title <- plot_grid(ds2_intersection_title01, ds2_intersection_title02, nrow = 2, rel_heights = c(1.2, 0.8))

ds2_intersection_subtitle01 <- ggdraw() + draw_label("Top 0,1%", fontface='bold', size = 16)
ds2_intersection_subtitle02 <- ggdraw() + draw_label("Top 0,5%", fontface='bold', size = 16)
ds2_intersection_subtitle <- plot_grid(ds2_intersection_subtitle01, ds2_intersection_subtitle02)

ds2_intersection_persnp_plt <- plot_grid(ds2_persnp_top_intersection.999_plt, ds2_persnp_top_intersection.995_plt)
ds2_intersection_pergene_plt <- plot_grid(ds2_pergene_top_intersection.999_plt, ds2_pergene_top_intersection.995_plt)

ds2_intersection_plt <- plot_grid(ds2_intersection_persnp_plt, ds2_intersection_pergene_plt, nrow = 2, labels = "AUTO")
ds2_intersection_plt <- plot_grid(ds2_intersection_title, ds2_intersection_subtitle, ds2_intersection_plt, ncol=1, rel_heights=c(0.13, 0.08, .9))

ggsave("figures/png/ds2_intersection.png", ds2_intersection_plt, height = 7, width = 9)
ggsave("figures/pdf/ds2_intersection.pdf", ds2_intersection_plt, height = 7, width = 9)
