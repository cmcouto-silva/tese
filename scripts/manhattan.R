manhattan_plot <- function(data, chr = "CHR", pos = "POS", snp = "SNP", gene = "GENE", annot = FALSE, score = "PVALUE", 
                           plot.title = "", y.axis.title = expression(-log[10]~(italic(p))), x.axis.title = "chromosome",
                           point.color = c("gray10", "gray60"), point.alpha = 0.8, point.size = 2,
                           text.size = 18, plot.title.size = 32,
                           threshold_lines = NULL, threshold_colors = NULL, threshold_mode = "quantile") {
  
  # Select target columns
  target_cols <- c(chr, pos)
  new_col_names <- c("CHR", "POS")
  
  # Keep SNP column
  if(!is.null(snp)) {
    target_cols <- append(target_cols, snp)
    new_col_names <- append(new_col_names, "SNP")
  }
  # Keep GENE column
  if(!is.null(gene)) {
    target_cols <- append(target_cols, gene)
    new_col_names <- append(new_col_names, "GENE")
  }
  # Keep Score column
  target_cols <- c(target_cols, score)
  new_col_names <- c(new_col_names, "SCORE")
  
  # Filter and Create Dataframe
  df <- data %>% select(all_of(target_cols)) %>%
    set_colnames(new_col_names) %>%
    mutate(CHR=as.integer(CHR),
           POS=as.numeric(POS),
           SCORE=as.numeric(SCORE))
  
  # Get Annotation column
  if(is.character(annot)) {
    df <- df %>%
      mutate(ANNOT=data[[annot]])
  }
  
  ## -- Set SNP Position -- ##
  df <- df %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(POS), .groups = 'drop') %>% 
    # Calculate cumulative position of each chromosome
    mutate(chr_tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, POS) %>%
    mutate(pos_cum=POS+chr_tot)
  
  x.axis.label <- df %>% group_by(CHR) %>% 
    summarize(center=(max(pos_cum) + min(pos_cum)) / 2, .groups = 'drop')
  x.axis.label$CHR[seq(19,22, 2)] <- ""
  
  ## -- Plot -- ##
  if(!is.null(threshold_colors)) {
    if(length(threshold_lines) != length(threshold_colors)) {
      stop("Length of threshold_colors must be equal to threshold_lines")
    }
    plt <- ggplot(df[SCORE < quantile(SCORE, min(threshold_lines))], aes(x=pos_cum, y=SCORE)) +
      
      # Show all points
      geom_point(aes(color=as.factor(CHR)), alpha=point.alpha, size=point.size) +
      scale_color_manual(values = rep(point.color, length.out=22)) +
      
      # custom X axis:
      scale_x_continuous(label = x.axis.label$CHR, breaks = x.axis.label$center)
    
    for(i in seq_along(threshold_lines)) {
      if(threshold_lines[i]!=threshold_lines[length(threshold_lines)]){
        df_tmp <- df[SCORE >= quantile(SCORE, threshold_lines[i]) & SCORE < quantile(SCORE, threshold_lines[i+1])]
      } else {
        df_tmp <- df[SCORE >= quantile(SCORE, threshold_lines[i])]
      }
      plt <- plt +
        geom_point(data = df_tmp, color = threshold_colors[i])
    }
  } else {
    plt <- ggplot(df, aes(x=pos_cum, y=SCORE)) +
      
      # Show all points
      geom_point(aes(color=as.factor(CHR)), alpha=point.alpha, size=point.size) +
      scale_color_manual(values = rep(point.color, length.out=22)) +
      
      # custom X axis:
      scale_x_continuous(label = x.axis.label$CHR, breaks = x.axis.label$center)
  }
  
  # Set custom theme
  plt <- plt +
    theme_bw() +
    theme ( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      # panel.grid.minor.x = element_blank(),
      axis.line.y.left = element_line(colour = "black"),
      axis.line.x.bottom = element_line(colour = "black"),
      text = element_text(size = text.size), 
      plot.title = element_text(face = "bold", size = plot.title.size, hjust = .5, 
                                margin = margin(b = 30, t = 10)),
      axis.title.y.left = element_text(margin = margin(r = 10)),
      axis.title.x.bottom = element_text(margin = margin(t = 10)),
      axis.text.x.bottom = element_text(margin = margin(t = 5)),
      axis.ticks.length.x.bottom = unit(.25, "cm")
    )
  
  # Set threshold lines
  if(!is.null(threshold_lines)) {
    plt <- plt + 
      geom_hline(yintercept = df[, quantile(SCORE, threshold_lines)], linetype="dashed")
  }
  
  # Set axis label names
  plt <- plt + 
    labs(title = plot.title, x = x.axis.title, y = y.axis.title)
  
  return(plt)
}

annotate_plot <- function(plt, threshold = 0, ..., groupby = "CHR", 
                          include_genes = NULL, exclude_gene_pattern = "^LOC|^LINC|^$",
                          y.upper.limit, nudge_y = .05, nudge_x = 0, direction = "y") {
  
  idx <- sapply(plt$layers, function(layer) "GeomPoint" %in% class(layer$geom))
  idx[1] <- FALSE
  
  geom_list <- lapply(plt$layers[idx], function(layer) {
    layer$data
  })
  
  geom_list <- append(list(plt$data), geom_list)
  df <- rbindlist(geom_list)
  
  if(groupby != "CHR" & threshold < 0.9) {
    stop("When grouping by other column than 'CHR', please set threshold >= 0.9.")
  }
  
  annot_df <- df %>% snpsel::separate_genes(setdt = F) %>%
    filter(SCORE >= quantile(SCORE, threshold)) %>%
    filter(!grepl(exclude_gene_pattern, ANNOT)) %>%
    group_by_at(groupby) %>%
    slice(which.max(SCORE)) %>%
    mutate(GENE=paste(GENE, collapse = ",")) %>%
    unique()
  
  if(!is.null(include_genes)) {
    annot_genes <- df[GENE %in% include_genes] %>% 
      group_by_at("GENE") %>%
      slice(which.max(SCORE)) %>%
      ungroup %>% filter(!GENE %in% annot_df$GENE)
    annot_df <- rbind.data.frame(annot_df, annot_genes)
    }
  
  annot_plt <- plt + 
    ggrepel::geom_label_repel (
      aes(pos_cum, SCORE, label=ANNOT), annot_df, ...,
      nudge_y = nudge_y, nudge_x = nudge_x, direction = direction
      )
  
  if(!missing(y.upper.limit)) {
    annot_plt <- annot_plt + ylim(0, y.upper.limit)
  }
  
  annot_plt$custom_attr <- list(annot=annot_df$ANNOT)
  return(annot_plt)  
  
}
