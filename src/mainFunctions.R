# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# main function of the GWAS engnine

#' Title
#'
#' @param genoFile path of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param phenoFile path of the phenotypic data file (`csv` file)
#' @param genoUrl url of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param phenoUrl url of the phenotypic data file (`csv` file)
#' @param trait Chraracter of length 1, name of the trait to analyze. Could be a
#'   column name of the phenotypic file
#' @param test Which test to use. Either `"score"`,  `"wald"` or `"lrt"`. For
#'   binary phenotypes, test = `"score"` is mandatory. For more information
#'   about this parameters see: `??gaston::association.test`
#' @param  fixed Number of Principal Components to include in the model with
#'   fixed effect (for test = `"wald"` or `"lrt"`). Default value is 0. For more
#'   information about this parameters see: `??gaston::association.test`
#' @param response Character of length 1, Either "quantitative" or "binary". Is
#'   the trait a quantitative or a binary phenotype? Default value is
#'   "quantitative"
#' @param thresh_maf Threshold for filtering markers. Only markers with minor
#'   allele frequency > `thresh_maf` will be kept.
#' @param thresh_callrate Threshold for filtering markers. Only markers with a
#'   callrate > `thresh_callrate` will be kept.
#' @param  dir directory where to save the data,
#' by default it is a temporary directory
#'
#' @return list with 2 elements `gwasRes` for the results of the gwas analysis in json and `file` path of the json file containing the results (id `dir` is not `NULL`)
run_gwas <- function(genoFile = NULL,
                     phenoFile = NULL,
                     genoUrl = NULL,
                     phenoUrl = NULL,
                     trait,
                     test,
                     fixed = 0,
                     response = "quantitative",
                     thresh_maf,
                     thresh_callrate,
                     dir = tempdir()){

  logger <- logger$new("r-run_gwas()")

  logger$log("Get data ...")
  if (!is.null(genoFile) && !is.null(phenoFile)
      && is.null(genoUrl) && is.null(phenoUrl)) {
    data <- readData(genoFile = genoFile,
                     phenoFile = phenoFile)
  } else if (!is.null(genoUrl) && !is.null(phenoUrl)
            && is.null(genoFile) && is.null(phenoFile)) {
    data <- downloadData(genoUrl = genoUrl,
                         phenoUrl = phenoUrl)
  } else {
    stop("Error: either genoFile and phenoFile or genoUrl and phenoUrl should be provided")
  }
  logger$log("Get data DONE")


  logger$log("GWAS analysis ...")
  gwasRes <- gwas(data = data,
                  trait = trait,
                  test = test,
                  fixed = fixed,
                  response = response,
                  thresh_maf = thresh_maf,
                  thresh_callrate = thresh_callrate)
  logger$log("GWAS analysis DONE")

  if (!is.null(dir)) {
  logger$log("Save results ...")
    file <- saveGWAS(gwas = gwasRes, dir = dir)
  logger$log("Save results DONE")
  } else {file <- NULL}

  return(list(
    "gwasRes" = jsonlite::toJSON(gwasRes, dataframe = "rows", pretty = T),
    "file" = file
  ))
}



#' Draw an interactive Manhattan Plot
#'
#' @param gwasFile path of the gwas result data file (json file)
#' @param gwasUrl url of the gwas result data file (json file)
#' @param adj_method correction method: "holm", "hochberg",
#' "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
#' @param thresh_p p value significant threshold (default 0.05)
#' @param chr name of the chromosome to show (show all if NA)
#' @param title Title of the plot. Default is "Manhattan Plot"
#'
#' @return plotly graph
draw_manhattanPlot <- function(gwasFile = NULL,
                               gwasUrl = NULL,
                               adj_method = "bonferroni",
                               thresh_p = 0.05,
                               chr = NA,
                               title = "Manhattan Plot") {
  logger <- logger$new("r-draw_manhattanPlot()")

  logger$log("Get data ...")
  if (!is.null(gwasFile) &&  is.null(gwasUrl)) {
    gwas <- readGWAS(gwasFile)
  } else if (!is.null(gwasUrl) && is.null(gwasFile)) {
    gwas <- downloadGWAS(gwasUrl)
  } else {
    stop("Error: either gwasFile or gwasUrl should be provided")
  }
  logger$log("Get data DONE")

  logger$log("Draw Manhattan Plot ...")
  p <- manPlot(gwas = gwas,
               adj_method = adj_method,
               thresh_p = thresh_p,
               chr = chr,
               title = title)
  logger$log("Draw Manhattan Plot DONE")

  p
}




#' Draw an LD Plot
#'
#' @param genoFile path of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param phenoFile path of the phenotypic data file (`csv` file)
#' @param from lower bound of the range of SNPs for which the LD is computed
#' @param to upper bound of the range of SNPs for which the LD is computed
#' @param dir path of the directory where to save the file, by default dir is set to a temporary directory, if null, the plot will be displayed.
#'
#' @return path of the created file (or NULL if `dir` is NULL)
draw_ldPlot <- function(genoFile = NULL,
                        genoUrl = NULL,
                        from,
                        to,
                        dir = tempdir()) {
  logger <- logger$new("r-draw_ldPlot()")

  logger$log("Get data ...")
  if (!is.null(genoFile) &&  is.null(genoUrl)) {
    geno <- readGenoData(genoFile)
  } else if (!is.null(genoUrl) && is.null(genoFile)) {
    geno <- downloadGenoData(genoUrl)
  } else {
    stop("Error: either genoFile or genoUrl should be provided")
  }
  logger$log("Get data DONE")

  logger$log("Draw LD Plot ...")
  imgFile <- LDplot(geno = geno,
                    from = from,
                    to = to,
                    dir = dir)
  logger$log("Draw LD Plot DONE")

  imgFile
}



#' Adjust GWAS p-values
#'
#' @param gwasFile path of the gwas result data file (json file)
#' @param gwasUrl url of the gwas result data file (json file)
#' @param adj_method correction method: "holm", "hochberg",
#' "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
#' @param  dir directory where to save the data,
#' by default it is a temporary directory
#'
#' @return list with 2 elements `gwasAdjusted` for the results of the gwas analysis in json with adjusted p-values and `file` path of the json file containing the results (if `dir` is not `NULL`)
run_resAdjustment <- function(gwasFile = NULL,
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              dir = tempdir()) {
  logger <- logger$new("r-run_resAdjustment()")

  logger$log("Get data ...")
  if (!is.null(gwasFile) &&  is.null(gwasUrl)) {
    gwas <- readGWAS(gwasFile)
  } else if (!is.null(gwasUrl) && is.null(gwasFile)) {
    gwas <- downloadGWAS(gwasUrl)
  } else {
    stop("Error: either gwasFile or gwasUrl should be provided")
  }
  logger$log("Get data DONE")

  # P-Values adjustment
  logger$log("Adjust p-values ...")
  adj <- adjustPval(gwas$p, adj_method)
  gwas$p_adj <- adj$p_adj
  logger$log("Adjust p-values DONE")

  if (!is.null(dir)) {
    logger$log("Save results ...")
    file <- saveGWAS(gwas = gwas, dir = dir)
    logger$log("Save results DONE")
  } else {file <- NULL}

  return(list(
    "gwasAdjusted" = jsonlite::toJSON(gwas, dataframe = "rows", pretty = T),
    "file" = file
  ))
}
