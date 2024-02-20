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
                       filter_pAdj = 1,
                       filter_nPoints = Inf,
                       filter_quant = 1) {

  logger <- Logger$new("r-filterGWAS()")


  # filter results ----
  logger$log("Filter points ...")
  nTotalSnp <- nrow(gwas)
  # filter according to a threshold on pAdj
  if (nrow(gwas) != 0 && filter_pAdj != 1) {
    if (is.null(gwas$p_adj)) {
      warning("gwas's p-values haven't been adjusted. Filtering according to",
              "p_ajd is not possible")
    } else if (filter_pAdj < 0 || filter_pAdj > 1) {
      bad_argument("filter_pAdj", expected = "between 0  and 1", not = filter_pAdj, "value")
    } else {
      gwas <- gwas[gwas$p_adj <= filter_pAdj,]
      if (nrow(gwas) == 0) {
        warning('filter_pAdj removed all the points of the graph')
      }
    }
  } else {
    logger$log('skip filter_pAdj')
  }

  # filter according to quantile
  if (nrow(gwas) != 0 && filter_quant != 1) {
    if (filter_quant < 0 || filter_quant > 1) {
      bad_argument('filter_quant', expected = "between 0  and 1", not = filter_quant, "value")
    }
    gwas <- gwas[order(gwas$p),]
    gwas <- gwas[seq_len(min(nrow(gwas),
                             floor(nTotalSnp*filter_quant))),]
    if (nrow(gwas) == 0) {
      warning('filter_quant removed all the points of the graph')
    }
  } else {
    logger$log('skip filter_quant')
  }

  # filter according to a fixed number of point
  if (nrow(gwas) != 0 && filter_nPoints < nrow(gwas)) {
    if (filter_nPoints < 0) {
      bad_argument("filter_nPoints", expected = "a positive number", not = filter_nPoints, "value")
    }
    gwas <- gwas[order(gwas$p),]
    gwas <- gwas[seq_len(min(nrow(gwas), filter_nPoints)),]
    if (nrow(gwas) == 0) {
      warning('filter_nPoints removed all the points of the graph')
    }
  } else {
    logger$log('skip filter_nPoints')
  }

  gwas
}
