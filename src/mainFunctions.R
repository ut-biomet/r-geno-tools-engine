# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# main function of the GWAS engnine

#' Run GWAS analysis
#'
#' @param genoFile path of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param phenoFile path of the phenotypic data file (`csv` file)
#' @param genoUrl url of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param phenoUrl url of the phenotypic data file (`csv` file)
#' @param trait Chraracter of length 1, name of the trait to analyze. Must be a
#'   column name of the phenotypic file.
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
#' @param outFile path of the output file. If `NULL`, the output will not be
#' written in any file. By default write in an tempoary `.json` file.
#'
#' @return list with 3 elements `gwasRes` for the results of the gwas analysis in json, `metadata` a list of metadata of these analysis and `file` path of the json file containing the results (id `dir` is not `NULL`)
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
                     outFile = tempfile(fileext = ".json")){

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

  logger$log("Save metadata ...")
  metadata <- list(genoFP = digest::digest(data$genoData),
                   phenoFP = digest::digest(data$phenoData),
                   trait = trait,
                   test = test,
                   fixed = fixed,
                   response = response,
                   thresh_maf = thresh_maf,
                   thresh_callrate = thresh_callrate,
                   date = Sys.time()
  )
  logger$log("Save metadata DONE")


  if (!is.null(outFile)) {
    logger$log("Save results ...")
    file <- saveGWAS(gwasRes = gwasRes, metadata = metadata, file = outFile)
    logger$log("Save results DONE")
  }

  return(list(
    "gwasRes" = jsonlite::toJSON(gwasRes, complex = "list", pretty = T),
    "metadata" = metadata,
    "file" = file
  ))
}



#' Draw a Manhattan Plot
#'
#' @param gwasFile path of the gwas result data file (json file)
#' @param gwasUrl url of the gwas result data file (json file)
#' @param adj_method correction method: "holm", "hochberg",
#' "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
#' @param thresh_p p value significant threshold (default 0.05)
#' @param chr name of the chromosome to show (show all if NA)
#' @param interactive [bool] should the plot be interactive (the default)
#' @param filter_pAdj [numeric] threshold to remove points
#' with pAdj < filter_pAdj from the plot (default no filtering)
#' @param filter_nPoints [numeric] threshold to keep only the filter_nPoints
#' with the lowest p-values for the plot (default no filtering)
#' @param filter_quant [numeric] threshold to keep only the filter_quant*100 %
#' of the points with the lowest p-values for the plot (default no filtering)
#' @param outFile path of the file containing the plot. If `NULL`, the
#' output will not be written in any file. By default write in an tempoary
#' file.
#'
#' @details If several filtering rules are given, the filtering process apply
#' the filtering process sequentially (this lead to having the same result
#' that if only the strongest rules were given).
#' Moreover, the number of points kept for the plot will be display
#' in the plot title.
#' @return plotly graph if interactive is TRUE, or NULL if not.
draw_manhattanPlot <- function(gwasFile = NULL,
                               gwasUrl = NULL,
                               adj_method = "bonferroni",
                               thresh_p = 0.05,
                               chr = NA,
                               interactive = TRUE,
                               filter_pAdj = 1,
                               filter_nPoints = Inf,
                               filter_quant = 1,
                               outFile = tempfile()) {
  logger <- logger$new("r-draw_manhattanPlot()")

  logger$log("Check outFile ...")
  if (!is.null(outFile)) {
    if (length(outFile) != 1) {
      stop("Error: `outFile` must be of length 1.")
    }
    if (tools::file_ext(outFile) == '') {
      if (interactive) {
        outFile <- paste0(outFile, '.html')
      } else {
        outFile <- paste0(outFile, '.png')
      }
    }
  }
  logger$log("Check outFile DONE")

  logger$log("Get data ...")
  if (!is.null(gwasFile) &&  is.null(gwasUrl)) {
    gwas <- readGWAS(gwasFile)
  } else if (!is.null(gwasUrl) && is.null(gwasFile)) {
    gwas <- downloadGWAS(gwasUrl)
  } else {
    stop("Error: either gwasFile or gwasUrl should be provided")
  }
  logger$log("Get data DONE")

  if (!interactive && !is.null(outFile)) {
    logger$log("Open connexion to draw the png plot ...")
    png(filename = outFile,
        width = 40,
        height = 30,
        units = 'cm',
        res = 177,
        pointsize = 18,
        type = "cairo")
    logger$log("Open connexion to draw the png plot DONE")
  }

  logger$log("Draw Manhattan Plot ...")
  p <- manPlot(gwas = gwas$gwas,
               adj_method = adj_method,
               thresh_p = thresh_p,
               chr = chr,
               interactive = interactive,
               title = paste(gwas$metadata$trait,
                             adj_method,
                             thresh_p, sep = " - "),
               filter_pAdj = filter_pAdj,
               filter_quant = filter_quant,
               filter_nPoints = filter_nPoints)
  logger$log("Draw Manhattan Plot DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
    if (interactive) {
      htmlwidgets::saveWidget(plotly::partial_bundle(p),
                              outFile, selfcontained = TRUE)
    } else {
      dev.off()
    }
    logger$log("Save results DONE")
  }
  p
}




#' Adjust GWAS p-values
#'
#' @param gwasFile path of the gwas result data file (json file)
#' @param gwasUrl url of the gwas result data file (json file)
#' @param adj_method correction method: "holm", "hochberg",
#' "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
#' @param filter_pAdj [numeric] threshold to remove points
#' with pAdj < filter_pAdj from the plot (default no filtering)
#' @param filter_nPoints [numeric] threshold to keep only the filter_nPoints
#' with the lowest p-values for the plot (default no filtering)
#' @param filter_quant [numeric] threshold to keep only the filter_quant*100 %
#' of the points with the lowest p-values for the plot (default no filtering)
#' @param outFile path of the output file. If `NULL`, the output will not be
#' written in any file. By default write in an tempoary `.json` file.
#'
#' @return list with 3 elements `gwasAdjusted` for the results of the gwas analysis in json with adjusted p-values, `metadata` a list of metadata of the gwas analysis in json with adjusted p-values, and `file` path of the json file containing the results (if `dir` is not `NULL`)
run_resAdjustment <- function(gwasFile = NULL,
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              filter_pAdj = 1,
                              filter_nPoints = Inf,
                              filter_quant = 1,
                              outFile = tempfile(fileext = ".json")){
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
  adj <- adjustPval(gwas$gwas$p, adj_method)
  gwas$gwas$p_adj <- adj$p_adj
  logger$log("Adjust p-values DONE")

  # filter results
  gwas$gwas <- filterGWAS(gwas = gwas$gwas,
                          filter_pAdj = filter_pAdj,
                          filter_nPoints = filter_nPoints,
                          filter_quant = filter_quant)

  metadata <- gwas$metadata
  metadata$adj_method <- adj_method

  if (!is.null(outFile)) {
    logger$log("Save results ...")
    file <- saveGWAS(gwasRes = gwas$gwas, metadata = metadata, file = outFile)
    logger$log("Save results DONE")
  }
  return(list(
    "gwasAdjusted" = jsonlite::toJSON(gwas$gwas, dataframe = "rows", pretty = T),
    "metadata" = metadata,
    "file" = file
  ))
}






#' Draw an LD Plot
#'
#' @param genoFile path of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param phenoFile path of the phenotypic data file (`csv` file)
#' @param from lower bound of the range of SNPs for which the LD is computed
#' @param to upper bound of the range of SNPs for which the LD is computed
#' @param outFile path of the png file to save the plot. If `NULL`, the image file will not be
#' created. By default write in an tempoary `.png` file.
#'
#' @return path of the created file (or NULL if `file` is NULL)
draw_ldPlot <- function(genoFile = NULL,
                        genoUrl = NULL,
                        from,
                        to,
                        outFile = tempfile(fileext = ".png")) {
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
                    file = outFile)
  logger$log("Draw LD Plot DONE")

  imgFile
}
