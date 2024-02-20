# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Functions for running GWAS analysis


#' perform GWAS analysis
#'
#' @param data List return by `prepareDta` function
#' @param trait Chraracter of length 1, name of the trait to analyze. Could be a
#'   column name of the phenotypic file
#' @param test Which test to use. Either `"score"`,  `"wald"` or `"lrt"`. For
#'   binary phenotypes, test = `"score"` is mandatory. For more information
#'   about this parameters see: `??gaston::association.test`
#' @param fixed Number of Principal Components to include in the model with
#'   fixed effect (for test = `"wald"` or `"lrt"`). Default value is 0. For more
#'   information about this parameters see: `??gaston::association.test`
#' @param response Character of length 1, Either "quantitative" or "binary". Is
#'   the trait a quantitative or a binary phenotype? Default value is
#'   "quantitative"
#' @param thresh_maf Threshold for filtering markers. Only markers with minor
#'   allele frequency > `thresh_maf` will be kept.
#' @param thresh_callrate Threshold for filtering markers. Only markers with a
#'   callrate > `thresh_callrate` will be kept.
#'
#' @details For the calculation, the genetic relationship matrix need to be calculated. This is done based on the genetic matrix standardized by the the genetic mean "mu" and the genetic variance "sigma", after having filtering according to the `thresh_maf` and `thresh_callrate`.
#' @return `data.frame`
gwas <- function(data,
                 trait,
                 test,
                 fixed = 0,
                 response = "quantitative",
                 thresh_maf,
                 thresh_callrate) {
  logger <- Logger$new("r-gwas()")

  ### Aggregate data in bed matrix ----
  logger$log("aggregate data in bed matrix ...")
  bm <- data$genoData
  bm@ped$pheno <- data$phenoData[, trait]
  logger$log("aggregate data in bed matrix DONE")

  ### FILTER INDS  ----
  empty <- FALSE

  # remove individuals with missing phenotypic values
  logger$log("remove individuals with missing phenotypic values ...")
  bm <- gaston::select.inds(bm, !is.na(pheno))
  if (nrow(bm) == 0) {
    warning('Filtering process removed all the individuals.')
    empty <- TRUE
  }
  logger$log("remove samples with missing phenotypic values DONE")


  ### FILTER SNPs ----
  # keep marker with a large enough MAF (>0.05)
  # and low missing rate (callrate>0.9)
  logger$log("filter SNPs ...")
  bm <- gaston::select.snps(bm, maf > thresh_maf)
  bm <- gaston::select.snps(bm, callrate > thresh_callrate)
  if (ncol(bm) == 0) {
    warning('Filtering process removed all the SNPs.')
    empty <- TRUE
  }
  logger$log("filter SNPs DONE")

  ### FAST RETURN ----
  if (empty) {
    resCols <- list("score" = c("chr", "pos", "id", "A1", "A2",
                                "freqA2", "score", "p"),
                    "wald" = c("chr", "pos", "id", "A1", "A2",
                               "freqA2", "h2", "beta", "sd","p"),
                    "lrt" = c("chr", "pos", "id", "A1", "A2",
                              "freqA2", "h2", "LRT", "p"))
    logger$log("DONE, return output.")
    gwa <- data.frame(matrix(ncol = length(resCols[[test]]), nrow = 0))
    colnames(gwa) <- resCols[[test]]
    # need to add one line with NA, in order to let `saveGWAS` and `readGWAS`
    # save and read the object as a data.frame. If not, readGWAS will consider
    # an empty list instead
    gwa[1,] <- NA
    return(gwa)
  }


  # calculate genetic relational matrix
  logger$log("calculate genetic relatinoal matrix ...")
  gaston::standardize(bm) <- "mu_sigma"
  K <- gaston::GRM(bm, autosome.only = FALSE)
  logger$log("calculate genetic relatinoal matrix DONE")


  ### FIT MODEL ----
  logger$log("fit model ...")
  if (test != "score") {
    gwa <- gaston::association.test(
      bm,
      method = "lmm",
      response = response,
      test = test,
      eigenK = eigen(K),
      p = fixed
    )
  } else {
    gwa <- gaston::association.test(
      bm,
      method = "lmm",
      response = response,
      test = test,
      K = K
    )
  }
  logger$log("fit model DONE")
  logger$log("DONE, return output.")

  return(gwa)

}



#' Adjust P-values for Multiple Comparisons
#'
#' @param vector of p-values
#' @param adj_method  correction method: "holm", "hochberg",
#' "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
#' @param thresh_p optional value of the p value significant threshold (default NULL)
#'
#' @details The method "hommel" is not implemented because it is too long to calculate.
#' @return list of two elements: "p_adj" vector of adjusted p values, "thresh_adj" the adjusted threshold (if thresh_p is preovided, NULL if not)
adjustPval <- function(p, adj_method, thresh_p = NULL){
  logger <- Logger$new("r-adjustPval()")

  # check adjMethod
  logger$log("Check adj_method ...")
  logger$log("Check adj_method DONE")

  # check p
  logger$log("Check p values ...")
  if (any(p < 0) || any(p > 1)) {
    # The `p.adjust` function will proceed to the calculation even if the p-values
    # are invalid. To avoid unexpected behaviour, I'd prefer this function
    # to crash in such case.
    bad_argument('p', must_be = "between 0 and 1", class = NULL)
  }
  logger$log("Check p values DONE")

  # P-Values adjustment
  logger$log("Adjust p-values ...")
  p_adj <- p.adjust(p,
                    method = adj_method,
                    n = length(p))
  logger$log("Adjust p-values DONE")

  if (!is.null(thresh_p)) {
    if (!adj_method == "none") {
    # get adjusted threshold
    logger$log("Adjust threshold ...")
    thresh_adj <- uniroot(
      function(log10p){
        # accuracy is bad if we don't transform with log10
        x <- 10^(log10p)
        p.adjust(c(x, p),
                 method = adj_method,
                 n = length(p) + 1)[1] - thresh_p
      },
      c(log10(2e-12), log10(1))
    )
    thresh_adj <- 10^(thresh_adj$root)
    logger$log("Adjust threshold DONE")
    } else {
      thresh_adj <- thresh_p
    }
  } else {
    thresh_adj <- NULL
  }
  logger$log("DONE, return output")

  list("p_adj" = p_adj,
       "thresh_adj" = thresh_adj)
}
