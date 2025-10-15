# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# filtering functions


#' Filter gwas results
#' @param gwas [data.frame] output of the gwas function
#' @param filter_pAdj [numeric] threshold to remove points
#' with pAdj < filter_pAdj from the plot (default no filtering)
#' @param filter_nPoints [numeric] threshold to keep only the filter_nPoints
#' with the lowest p-values for the plot (default no filtering)
#' @param filter_quant [numeric] threshold to keep only the filter_quant*100 %
#' of the points with the lowest p-values for the plot (default no filtering)
filterGWAS <- function(gwas,
                       filter_pAdj = NULL,
                       filter_nPoints = NULL,
                       filter_quant = NULL) {
  logger <- Logger$new("r-filterGWAS()")

  # filter results ----
  logger$log("Filter points ...")

  # filter according to quantile
  # NOTE: We should filter on quantile first because it depends on the number of
  # rows in gwas
  if (!is.null(filter_quant)) {
    logger$log("Filtering on quantile...")
    gwas <- filterGWAS_quant(gwas, filter_quant)
  } else {
    logger$log("Skip filtering on quantile...")
  }

  # filter according to a threshold on pAdj
  if (!is.null(filter_pAdj)) {
    logger$log("Filtering on pAdj...")
    gwas <- filterGWAS_pAdj(gwas, filter_pAdj)
  } else {
    logger$log("Skip filtering on pAdj...")
  }

  # filter according to a fixed number of point
  if (!is.null(filter_nPoints)) {
    logger$log("Filtering on number of point...")
    gwas <- filterGWAS_nPoints(gwas, filter_nPoints)
  } else {
    logger$log("Skip filtering on number of point...")
  }

  logger$log("Filter points ...")
  return(gwas)
}


filterGWAS_pAdj <- function(gwas, filter_pAdj = NULL) {
  if (is.null(filter_pAdj) || filter_pAdj == 1 || nrow(gwas) == 0) {
    return(gwas)
  }
  if (is.null(gwas$p_adj)) {
    warning(
      "gwas's p-values haven't been adjusted. Filtering according to",
      "p_ajd is not possible, no filtering applied"
    )
    return(gwas)
  }
  if (filter_pAdj < 0 || filter_pAdj > 1) {
    bad_argument("filter_pAdj", must_be = "between 0  and 1", not = filter_pAdj)
  }
  gwas <- gwas[ifelse(is.na(gwas$p_adj), FALSE, gwas$p_adj <= filter_pAdj), ]
  if (nrow(gwas) == 0) {
    warning("filter_pAdj removed all the points of the graph")
  }
  return(gwas)
}


filterGWAS_nPoints <- function(gwas, filter_nPoints = NULL) {
  if (is.null(filter_nPoints) || nrow(gwas) == 0) {
    return(gwas)
  }
  if (filter_nPoints < 0) {
    bad_argument("filter_nPoints", must_be = "a positive number", not = filter_nPoints)
  }
  if (filter_nPoints >= nrow(gwas)) {
    return(gwas)
  }
  gwas <- gwas[order(gwas$p), , drop = FALSE]
  gwas <- gwas[seq_len(filter_nPoints), , drop = FALSE]
  if (nrow(gwas) == 0) {
    warning("filter_pAdj removed all the points of the graph")
  }
  return(gwas)
}


filterGWAS_quant <- function(gwas, filter_quant = NULL) {
  if (is.null(filter_quant) || filter_quant == 1 || nrow(gwas) == 0) {
    return(gwas)
  }
  if (filter_quant < 0 || filter_quant > 1) {
    bad_argument("filter_quant", must_be = "between 0  and 1", not = filter_quant)
  }
  gwas <- gwas[order(gwas$p), , drop = FALSE]
  nPoints <- floor(nrow(gwas) * filter_quant)
  gwas <- gwas[seq_len(nPoints), , drop = FALSE]
  if (nrow(gwas) == 0) {
    warning("filter_pAdj removed all the points of the graph")
  }
  return(gwas)
}
