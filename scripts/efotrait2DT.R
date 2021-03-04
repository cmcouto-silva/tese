efolist2DT <- function(efolist) {
  # Remove genes without annotation
  efolist <- efolist[sapply(efolist, function(dt) nrow(dt) > 0)]
  # Put gene name on DT
  invisible(lapply(seq_along(efolist), function(i) {
    efolist[[i]][, GENE := names(efolist)[i]]
  }))
  # Concatenate list into a data.table
  dt <- rbindlist(efolist)
  # Make unique records per trait
  dt <- dt[, .(
    N_STUDIES = sum(N),
    N_GENE = length(unique(GENE)),
    GENES = paste(unique(GENE), collapse = ",")),
    by = MAPPED_TRAIT]
  dt <- dt[MAPPED_TRAIT!=""] # remove empty trait
  dt <- dt[order(-N_GENE, -N_STUDIES, MAPPED_TRAIT)] # order data.table
  return(dt)
}
