# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# main function of the GWAS engnine

#' Run GWAS analysis
#'
#' @param genoFile path of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param phenoFile path of the phenotypic data file (`csv` file). Individuals'
#' name should be the first column of the file and no duplication is allowed.
#' @param genoUrl url of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param phenoUrl url of the phenotypic data file (`csv` file) Individuals'
#' name should be the first column of the file and no duplication is allowed.
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
#' @return list with 3 elements `gwasRes` for the results of the gwas analysis in json, `metadata` a list of metadata of these analysis and `file` path of the json file containing the results
run_gwas <- function(genoFile = NULL,
                     phenoFile = NULL,
                     genoUrl = NULL,
                     phenoUrl = NULL,
                     trait,
                     test,
                     fixed = NULL,
                     response = "quantitative",
                     thresh_maf,
                     thresh_callrate,
                     outFile = tempfile(fileext = ".json")){

  logger <- Logger$new("r-run_gwas()")

  check_test(test)
  check_fixed(fixed, test)
  if (test %in% c("wald", "lrt") && is.null(fixed)) {
    fixed <- 0
  }
  check_response(response)
  check_outFile(outFile, accept_null = TRUE)

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
    engineError("Either `genoFile` and `phenoFile` or `genoUrl` and `phenoUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "genotype"))
  }
  logger$log("Get data DONE")

  check_trait(trait, colnames(data$phenoData))
  check_thresh_maf(thresh_maf, geno = data$genoData)
  check_thresh_callrate(thresh_callrate, geno = data$genoData)



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
  } else {
    file <- NULL
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
  logger <- Logger$new("r-draw_manhattanPlot()")

  check_adj_method(adj_method)
  check_interactive(interactive)

  thresh_p <- as.numeric(thresh_p)
  check_thresh_p(thresh_p)

  filter_pAdj <- as.numeric(filter_pAdj)
  check_filter_pAdj(filter_pAdj)

  filter_nPoints <- as.numeric(filter_nPoints)
  check_filter_nPoints(filter_nPoints)

  filter_quant <- as.numeric(filter_quant)
  check_filter_quant(filter_quant)


  check_outFile(outFile, accept_null = TRUE)
  if (!is.null(outFile)) {
    ext <- tools::file_ext(outFile)
    if (ext == "") {
      if (interactive) {
        outFile <- paste0(outFile, '.html')
      } else {
        outFile <- paste0(outFile, '.png')
      }
    }
  }


  logger$log("Get data ...")
  if (!is.null(gwasFile) &&  is.null(gwasUrl)) {
    gwas <- readGWAS(gwasFile)
  } else if (!is.null(gwasUrl) && is.null(gwasFile)) {
    gwas <- downloadGWAS(gwasUrl)
  } else {
    engineError("Either `gwasFile` or `gwasUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "GWAS"))
  }
  logger$log("Get data DONE")

  check_chr(chr, unique(gwas$chr))

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
      save_plotly(p, outFile)
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
  logger <- Logger$new("r-run_resAdjustment()")

  check_adj_method(adj_method)

  filter_pAdj <- as.numeric(filter_pAdj)
  check_filter_pAdj(filter_pAdj)

  filter_nPoints <- as.numeric(filter_nPoints)
  check_filter_nPoints(filter_nPoints)

  filter_quant <- as.numeric(filter_quant)
  check_filter_quant(filter_quant)

  check_outFile(outFile, accept_null = TRUE)

  logger$log("Get data ...")
  if (!is.null(gwasFile) &&  is.null(gwasUrl)) {
    gwas <- readGWAS(gwasFile)
  } else if (!is.null(gwasUrl) && is.null(gwasFile)) {
    gwas <- downloadGWAS(gwasUrl)
  } else {
    engineError("Either `gwasFile` or `gwasUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "GWAS"))
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
  } else {
    file <- NULL
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
#' @param from lower bound of the range of SNPs for which the LD is computed
#' @param to upper bound of the range of SNPs for which the LD is computed
#' @param outFile path of the png file to save the plot. If `NULL`, the image file will not be
#' @param n_max maximum number of marker to show (to avoid unreadable plot)
#' created. By default write in an tempoary `.png` file.
#'
#' @return path of the created file (or NULL if `file` is NULL)
draw_ldPlot <- function(genoFile = NULL,
                        genoUrl = NULL,
                        from,
                        to,
                        n_max = 50,
                        outFile = tempfile(fileext = ".png")) {
  logger <- Logger$new("r-draw_ldPlot()")

  check_from_to(from, to, n_max)

  check_outFile(outFile, accept_null = TRUE)
  if (!is.null(outFile)) {
    ext <- tools::file_ext(outFile)
    if (ext == "") {
        outFile <- paste0(outFile, '.png')
    }
  }

  logger$log("Get data ...")
  if (!is.null(genoFile) &&  is.null(genoUrl)) {
    geno <- readGenoData(genoFile)
  } else if (!is.null(genoUrl) && is.null(genoFile)) {
    geno <- downloadGenoData(genoUrl)
  } else {
    engineError("Either `genoFile` or `genoUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "genotype"))
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


#' Calculate pedigree relationship matrix
#'
#' @param pedFile path of the pedigree data file (`csv` file).
#' @param pedUrl url of the pedigree data file (`csv` file).
#' @param unknown_string  [default: ""] a character vector of strings which are
#' to be interpreted as "unknown parent". By default: missing value in the file.
#' @param header  [default: TRUE] a logical value indicating whether the file
#' contains the names of the variables as its first line. The default value is
#' TRUE. In any cases, the column 1 will be interpreted as the individual id,
#' column 2 as the first parent, column 3 as the second parent.
#' @param outFile path of the output file. If `NULL`, the output will not be
#' written in any file. By default write in an tempoary `.json` file.
#' @param outFormat Format of the output file, either `csv` or `json`.
#' by default it will use the file extension of `outfile`.
#'
#' @details
#' For `csv` output, the file will include some metadata lines (starting by a `#`
#' symbol), a header and the row ids in its first column.
#'
#' @return list with 3 elements `relMat` the relationship matrix, `metadata` a
#' list of metadata of these analysis (pedigree fingerprint,
#' number of individuals, creation time) and `file` path
#' of the file containing the results.
calc_pedRelMat <- function(pedFile = NULL,
                           pedUrl = NULL,
                           unknown_string = "",
                           header = TRUE,
                           outFile = tempfile(fileext = ".csv"),
                           outFormat = tools::file_ext(outFile)) {
  logger <- Logger$new("r-calc_pedRelMAt()")

  check_outFile(outFile, accept_null = TRUE)

  logger$log("Get data ...")
  if (!is.null(pedFile) &&  is.null(pedUrl)) {
    ped <- readPedData(file = pedFile,
                       unknown_string = unknown_string,
                       header = header)
  } else if (!is.null(pedUrl) && is.null(pedFile)) {
    ped <- downloadPedData(url = pedUrl,
                           unknown_string = unknown_string,
                           header = header)
  } else {
    engineError("Either `pedFile` or `pedUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "pedigree"))
  }
  logger$log("Get data DONE")


  logger$log("Calcualte pedigree relationship matrix ...")
  relMat <- pedRelMat(ped = ped)
  logger$log("Calcualte pedigree relationship matrix DONE")
  logger$log("Get metadata ...")
  metadata <- list(info = "R-geno-engine, Pedigree relationship matrix",
                   date = Sys.time(),
                   nInds = ncol(relMat),
                   pedFP = digest::digest(ped$data))
  logger$log("Get metadata DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
    file <- saveRelMat(relMat = relMat, metadata = metadata, file = outFile, format = outFormat)
    logger$log("Save results DONE")
  } else {
    file <- NULL
  }

  return(list(
    "relMat" = relMat,
    "metadata" = metadata,
    "file" = file
  ))
}





#' Calculate genomic relationship matrix
#'
#' @param genoFile path of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param genoUrl url of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param outFile path of the output file. If `NULL`, the output will not be
#' written in any file. By default write in an tempoary `.json` file.
#' @param outFormat Format of the output file, either `csv` or `json`.
#' by default it will use the file extension of `outfile`.
#'
#' @details
#' For `csv` output, the file will include some metadata lines (starting by a `#`
#' symbol), a header and the row ids in its first column.
#'
#' @return list with 3 elements `relMat` the relationship matrix, `metadata` a
#' list of metadata of these analysis (pedigree fingerprint,
#' number of individuals, creation time) and `file` path
#' of the file containing the results.
calc_genoRelMat <- function(genoFile = NULL,
                            genoUrl = NULL,
                            outFile = tempfile(fileext = ".csv"),
                            outFormat = tools::file_ext(outFile)) {
  logger <- Logger$new("r-calc_genoRelMAt()")

  check_outFile(outFile, accept_null = TRUE)

  logger$log("Get data ...")
  if (!is.null(genoFile) &&  is.null(genoUrl)) {
    geno <- readGenoData(genoFile)
  } else if (!is.null(genoUrl) && is.null(genoFile)) {
    geno <- downloadGenoData(genoUrl)
  } else {
    engineError("Either `genoFile` or `genoUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "genotype"))
  }

  logger$log("Calcualte genomic relationship matrix ...")
  relMat <- genoRelMat(geno = geno)
  logger$log("Calcualte genomic relationship matrix DONE")
  logger$log("Get metadata ...")
  metadata <- list(info = "R-geno-engine, genomic relationship matrix",
                   date = Sys.time(),
                   nInds = ncol(relMat),
                   genoFP = digest::digest(geno))
  logger$log("Get metadata DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
    file <- saveRelMat(relMat = relMat, metadata = metadata, file = outFile, format = outFormat)
    logger$log("Save results DONE")
  } else {
    file <- NULL
  }

  return(list(
    "relMat" = relMat,
    "metadata" = metadata,
    "file" = file
  ))
}



#' Combined (pedigree + genomic) Relationship Matrix
#'
#' Correct a pedigree relationship matrix by using genomic relationship matrix.
#'
#' @param pedRelMatFile path of a pedigree relationship matrix generated by the
#' the engine.
#' @param pedRelMatUrl url of a pedigree relationship matrix generated by the
#' the engine.
#' @param genoRelMatFile path of a genomic relationship matrix generated by the
#' the engine.
#' @param genoRelMatUrl url of a genomic relationship matrix generated by the
#' the engine.
#' @param method method to use, either "Legarra" or "Martini"
#' @param tau tau parameter of the Martini's method
#' @param omega omega parameter of the Martini's method
#' @param outFile path of the output file. If `NULL`, the output will not be
#' written in any file. By default write in an tempoary `.json` file.
#' @param outFormat Format of the output file, either `csv` or `json`.
#' by default it will use the file extension of `outfile`.
#'
#' @details
#' This method correct the pedigree matrix with the genomic relationship matrix.
#' Therefore, individuals in the genomic relationship matrix not in the pedigree
#' relationship matrix will be ignored.
#'
#' Using the Martini's method with `tau=1`, and `omega=1` is equivalent of
#' Legarra's method.
#'
#' For `csv` output, the file will include some metadata lines (starting by a `#`
#' symbol), a header and the row ids in its first column.
#'
#' @references
#' *Martini, JW, et al. 2018 The effect of the H-1 scaling factors tau and omega on the structure of H in the single-step procedure. Genetics Selection Evolution 50(1), 16*
#'
#' *Legarra, A, et al. 2009 A relationship matrix including full pedigree and genomic information. Journal of Dairy Science 92, 4656–4663*
#'
#' @return list with 3 elements `relMat` the relationship matrix, `metadata` a
#' list of metadata of these analysis (pedigree fingerprint,
#' number of individuals, creation time) and `file` path
#' of the file containing the results.
calc_combinedRelMat <- function(pedRelMatFile = NULL,
                                pedRelMatUrl = NULL,
                                genoRelMatFile = NULL,
                                genoRelMatUrl  = NULL,
                                method = "Legarra",
                                tau = NULL,
                                omega = NULL,
                                outFile = tempfile(fileext = ".csv"),
                                outFormat = tools::file_ext(outFile)) {

  logger <- Logger$new("r-calc_combinedRelMat()")

  check_method(method, tau, omega)
  if (method == "Martini") {check_tau_omega(tau, omega)}
  check_outFile(outFile, accept_null = TRUE)

  logger$log("Get data ...")
  if (!is.null(pedRelMatFile) &&  is.null(pedRelMatUrl)) {
    ped_rm <- readRelMat(pedRelMatFile)
  } else if (!is.null(pedRelMatUrl) && is.null(pedRelMatFile)) {
    ped_rm <- downloadRelMat(pedRelMatUrl)
  } else {
    engineError("Either `pedRelMatFile` or `pedRelMatUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "pedigree relationship matrix"))
  }
  if (!is.null(genoRelMatFile) &&  is.null(genoRelMatUrl)) {
    geno_rm <- readRelMat(genoRelMatFile)
  } else if (!is.null(genoRelMatUrl) && is.null(genoRelMatFile)) {
    geno_rm <- downloadRelMat(genoRelMatUrl)
  } else {
    engineError("Either `genoRelMatFile` or `genoRelMatUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "genomic relationship matrix"))
  }
  logger$log("Get data DONE")

  if (!any(rownames(geno_rm) %in% rownames(ped_rm))) {
    engineError(paste('No common individuals between genomic',
                      'and pedigree relationship matrices.'),
                extra = list(code = errorCode("INCONSISTENT_RELATIONSHIP_MATRICES")))
  }

  logger$log("Calcualte combined relationship matrix ...")
  relMat <- combinedRelMat(ped_rm = ped_rm,
                           geno_rm = geno_rm,
                           method = method,
                           tau = tau,
                           omega = omega)
  logger$log("Calcualte combined relationship matrix DONE")
  logger$log("Get metadata ...")
  metadata <- list(info = "R-geno-engine, combined relationship matrix",
                   date = Sys.time(),
                   nInds = ncol(relMat),
                   geno_relMatFP = digest::digest(geno_rm),
                   ped_relMatFP = digest::digest(ped_rm))
  logger$log("Get metadata DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
    file <- saveRelMat(relMat = relMat, metadata = metadata, file = outFile, format = outFormat)
    logger$log("Save results DONE")
  } else {
    file <- NULL
  }

  return(list(
    "relMat" = relMat,
    "metadata" = metadata,
    "file" = file
  ))
}





#' Draw a heatmap of a relationship matrix
#'
#' @param relMatFile path of a file generated by the function `saveRelMat`
#' @param relMatUrl url of a file generated by the function `saveRelMat`
#' @param interactive [bool] should the plot be interactive (the default) or not
#' @param outFile path of the file containing the plot. If `NULL`, the
#' output will not be written in any file. By default write in an tempoary
#' file.
#'
#' @return plotly graph if interactive is TRUE, or NULL if not.
draw_relHeatmap <- function(relMatFile = NULL,
                            relMatUrl = NULL,
                            format = NULL,
                            interactive = TRUE,
                            outFile = tempfile()) {
  logger <- Logger$new("r-draw_relHeatmap()")

  check_interactive(interactive)

  check_outFile(outFile, accept_null = TRUE)
  if (!is.null(outFile)) {
    ext <- tools::file_ext(outFile)
    if (ext == "") {
      if (interactive) {
        outFile <- paste0(outFile, '.html')
      } else {
        outFile <- paste0(outFile, '.png')
      }
    }
  }


  logger$log("Get data ...")
  if (!is.null(relMatFile) &&  is.null(relMatUrl)) {
    if (is.null(format)) {
      format <- tools::file_ext(relMatFile)
    }
    relMat <- readRelMat(file = relMatFile, format = format)
  } else if (!is.null(relMatUrl) && is.null(relMatFile)) {
    if (is.null(format)) {
      format <- tools::file_ext(relMatUrl)
    }
    relMat <- downloadRelMat(url = relMatUrl, format = format)
  } else {
    engineError("Either `relMatFile` or `relMatFile` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "relationship matrix"))
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

  logger$log("Draw relationship heatmap ...")
  p <- relMatHeatmap(relMat = relMat,
                     interactive = interactive)
  logger$log("Draw relationship heatmap DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
    if (interactive) {
      save_plotly(p, outFile)
    } else {
      dev.off()
    }
    logger$log("Save results DONE")
  }
  p
}







#' Draw interactive pedigree network
#'
#' @param pedFile path of the pedigree data file (`csv` file).
#' @param pedUrl url of the pedigree data file (`csv` file).
#' @param unknown_string  [default: ""] a character vector of strings which are
#' to be interpreted as "unknown parent". By default: missing value in the file.
#' @param header  [default: TRUE] a logical value indicating whether the file
#' contains the names of the variables as its first line. The default value is
#' TRUE. In any cases, the column 1 will be interpreted as the individual id,
#' column 2 as the first parent, column 3 as the second parent.
#' @param outFile path of the file containing the plot. If `NULL`, the
#' output will not be written in any file. By default write in an temporary
#' file.
#'
#' @return plotly graph if interactive is TRUE, or NULL if not.
draw_pedNetwork <- function(pedFile = NULL,
                            pedUrl = NULL,
                            unknown_string = "",
                            header = TRUE,
                            outFile = tempfile(fileext = ".html")) {
  logger <- Logger$new("r-draw_pedNetwork()")

  logger$log("Get data ...")
  if (!is.null(pedFile) &&  is.null(pedUrl)) {
    ped <- readPedData(file = pedFile,
                       unknown_string = unknown_string,
                       header = header)
  } else if (!is.null(pedUrl) && is.null(pedFile)) {
    ped <- downloadPedData(url = pedUrl,
                           unknown_string = unknown_string,
                           header = header)
  } else {
    engineError("Either `pedFile` or `pedUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "pedigree"))
  }
  logger$log("Get data DONE")

  logger$log("Draw pedigree interactive network ...")
  p <- pedNetwork(ped)
  logger$log("Draw pedigree interactive network DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
      save_plotly(p, outFile)
    logger$log("Save results DONE")
  }
  p
}





#' Simulate the genotypes of offspring given the parent genotypes.
#'
#' @param genoFile phased VCF file path (ext `.vcf` or `.vcf.gz`)
#' @param genoUrl url of a phased VCF file path (ext `.vcf` or `.vcf.gz`)
#' @param crossTableFile path of the crossing table data file
#' (`csv` file of 2 or 3 columns). It must contain the names of the variables
#' as its first line. The column 1 and 2 will be interpreted as the parents
#' ids. The optional third column will be interpreted as the offspring base
#' name.
#' @param crossTableUrl URL of a crossing table file
#' @param SNPcoordFile path of the SNPs coordinates file
#' (`csv` file). This `.csv` file should have 4 named columns:
#' - `chr`: Chromosome holding the SNP
#' - `physPos`: SNP physical position on the chromosome
#' - `linkMapPos`: SNP linkage map position on the chromosome in Morgan
#' - `SNPid`: SNP's IDs
#' @param SNPcoordUrl URL of a SNP coordinate file
#' @param nCross number of cross to simulate for each parent pair defined
#' in the crossing table.
#' @param outFile path of the `.vcf.gz` file containing the simulated genotypes
#' of the offspring. It must end by `.vcf.gz`. By default write in an temporary
#' file.
#'
#' @return path of the `.vcf.gz` file containing the simulated genotypes
#' of the offspring.
crossingSimulation <- function(genoFile = NULL,
                               genoUrl = NULL,
                               crossTableFile = NULL,
                               crossTableUrl = NULL,
                               SNPcoordFile = NULL,
                               SNPcoordUrl = NULL,
                               nCross = 30,
                               outFile = tempfile(fileext = ".vcf.gz")) {
  logger <- Logger$new("r-crossingSimulation()")

  check_outFile(outFile, expected_ext = "vcf.gz")
  check_nCross(nCross)


  # load input data
  logger$log("Get data ...")
  if (!is.null(genoFile) &&  is.null(genoUrl)) {
    g <- readPhasedGeno(genoFile)
  } else if (is.null(genoFile) && !is.null(genoUrl)) {
    g <- downloadPhasedGeno(genoUrl)
  } else {
    engineError("Either `genoFile` or `genoUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "genotype"))

  }

  if (!is.null(SNPcoordFile) &&  is.null(SNPcoordUrl)) {
    SNPcoord <- readSNPcoord(SNPcoordFile)
  } else if (is.null(SNPcoordFile) && !is.null(SNPcoordUrl)) {
    SNPcoord <- downloadSNPcoord(SNPcoordUrl)
  } else {
    engineError("Either `SNPcoordFile` or `SNPcoordUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "SNP coordinates"))
  }

  if (!is.null(crossTableFile) &&  is.null(crossTableUrl)) {
    crossTable <- readCrossTable(crossTableFile)
  } else if (is.null(crossTableFile) && !is.null(crossTableUrl)) {
    crossTable <- downloadCrossTable(crossTableUrl)
  } else {
    engineError("Either `crossTableFile` or `crossTableUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "crossing table"))
  }
  crossTable$n <- nCross
  logger$log("Get data DONE")

  # check input data
  logger$log("Check SNP's coordinates consistency between",
             "`.vcf` and SNPcoordinate file ...")
  SNPcoord <- checkAndFilterSNPcoord(SNPcoord, g$SNPcoord)
  g$SNPcoord <- NULL # release some memory
  logger$log("Check SNP's coordinates consistency between",
             " `.vcf` and SNPcoordinate file DONE")

  logger$log("Check individuals' names consistency between",
             " `.vcf` and `.csv` file ...")
  checkIndNamesConsistency(crossTable, g$haplotypes)
  logger$log("Check individuals' names consistency between",
             " `.vcf` and `.csv` file DONE")



  # Initialise sumulation population
  logger$log("Initialise simulation ...")
  parentPopulation <- initializeSimulation(haplotypes = g$haplotypes,
                                           SNPcoord = SNPcoord)
  logger$log("Initialise simulation DONE")

  # Simulate crosses
  logger$log("Crossing simulation ...")
  simulatedIndividuals <- breedSimulatR::makeCrosses(crosses = crossTable,
                                                     pop = parentPopulation)
  logger$log("Crossing simulation DONE")

  # write results
  logger$log("Write output file ...")
  simulatedPop <- breedSimulatR::population$new(inds = simulatedIndividuals,
                                                verbose = FALSE)

  if (file.exists(outFile)) {
    file.remove(outFile)
  }
  saveVcf(outFile, simulatedPop, SNPcoord)
  logger$log("Write output file DONE")


  return(outFile)
}





#' Calculate progenies BLUPs variance and expected values based on parents'
#' genotype and markers effects
#'
#' @param genoFile phased VCF file path (ext `.vcf` or `.vcf.gz`)
#' @param genoUrl url of a phased VCF file path (ext `.vcf` or `.vcf.gz`)
#' @param crossTableFile path of the crossing table data file
#' (`csv` file of 2 or 3 columns). It must contain the names of the variables
#' as its first line. The column 1 and 2 will be interpreted as the parents
#' ids. The optional third column will be interpreted as the offspring base
#' name.
#' @param crossTableUrl URL of a crossing table file
#' @param SNPcoordFile path of the SNPs coordinates file
#' (`csv` file). This `.csv` file should have 4 named columns:
#' - `chr`: Chromosome holding the SNP
#' - `physPos`: SNP physical position on the chromosome
#' - `linkMapPos`: SNP linkage map position on the chromosome in Morgan
#' - `SNPid`: SNP's IDs
#' @param SNPcoordUrl URL of a SNP coordinate file
#' @param markerEffectsFile path of the marker effects file (`csv` or `json` file).
#' @param markerEffectsUrl URL of a marker effect file
#' @param outFile `.json` file path where to save the data. If the file already exists,
#' it will be overwritten.
#'
#' @return data.frame containing the calculations results
calc_progenyBlupEstimation <- function(genoFile = NULL,
                                       genoUrl = NULL,
                                       crossTableFile = NULL,
                                       crossTableUrl = NULL,
                                       SNPcoordFile = NULL,
                                       SNPcoordUrl = NULL,
                                       markerEffectsFile = NULL,
                                       markerEffectsUrl = NULL,
                                       outFile = tempfile(fileext = ".json")) {
  logger <- Logger$new("r-progenyBlupVarExp()")

  check_outFile(outFile, accept_null = TRUE)

  # load input data
  logger$log("Get data ...")
  if (!is.null(genoFile) &&  is.null(genoUrl)) {
    g <- readPhasedGeno(genoFile)
  } else if (is.null(genoFile) && !is.null(genoUrl)) {
    g <- downloadPhasedGeno(genoUrl)
  } else {
    engineError("Either `genoFile` or `genoUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "genotype"))
  }

  if (!is.null(SNPcoordFile) &&  is.null(SNPcoordUrl)) {
    SNPcoord <- readSNPcoord(SNPcoordFile)
  } else if (is.null(SNPcoordFile) && !is.null(SNPcoordUrl)) {
    SNPcoord <- downloadSNPcoord(SNPcoordUrl)
  } else {
    engineError("Either `SNPcoordFile` or `SNPcoordUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "SNP coordinates"))

  }

  if (!is.null(crossTableFile) &&  is.null(crossTableUrl)) {
    crossTable <- readCrossTable(crossTableFile)
  } else if (is.null(crossTableFile) && !is.null(crossTableUrl)) {
    crossTable <- downloadCrossTable(crossTableUrl)
  } else {
    engineError("Either `crossTableFile` or `crossTableUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "crossing table"))
  }

  if (!is.null(markerEffectsFile) &&  is.null(markerEffectsUrl)) {
    markerEffects <- readMarkerEffects(markerEffectsFile)
  } else if (is.null(markerEffectsFile) && !is.null(markerEffectsUrl)) {
    markerEffects <- downloadMarkerEffects(markerEffectsUrl)
  } else {
    engineError("Either `markerEffectsFile` or `markerEffectsUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "markers effects"))
  }
  logger$log("Get data DONE")

  # check input data
  logger$log("Check individuals' names consistency between",
             " `.vcf` and `.csv` file ...")
  checkIndNamesConsistency(crossTable, g$haplotypes)
  logger$log("Check individuals' names consistency between",
             " `.vcf` and `.csv` file DONE")

  logger$log("Check SNP's coordinates consistency between",
             "`.vcf` and SNPcoordinate file ...")
  SNPcoord <- checkAndFilterSNPcoord(SNPcoord, g$SNPcoord)
  g$SNPcoord <- NULL # release some memory
  logger$log("Check SNP's coordinates consistency between",
             " `.vcf` and SNPcoordinate file DONE")

  logger$log("Check SNPs' ids consistency between",
             "SNPcoordinate and markerEffects file ...")
  if (!all(SNPcoord$SNPid %in% row.names(markerEffects$SNPeffects))) {
    missingSNP <- which(!SNPcoord$SNPid %in% row.names(markerEffects$SNPeffects))
    msg <- paste('The SNPs coordinate file miss', length(missingSNP),
                 'SNPs when compared with the provided marker effects.',)
    engineError(msg,
      extra = list(
        "code" = errorCode("SNP_COORD_MISSING_SNP"),
        "n_missing_SNP" = length(missingSNP),
        "n_expected_SNP" = nrow(markerEffects$SNPeffects),
        "missing_SNP" = missingSNP,
        "reference" = "marker effects"
    ))
  }
  logger$log("Check SNPs' ids consistency between",
             "SNPcoordinate and markerEffects file DONE")


  # calculate recombination rate
  r <- calcRecombRate(SNPcoord)

  # initialize results
  nCrosses <- nrow(crossTable)
  blupVarExp <- vector(mode = 'list', length = nCrosses)
  names(blupVarExp) <- crossTable$names

  # calculation for each couple
  logger$log("BLUP variance and expected value calculation for each crosses ...")
  for (i in seq(nCrosses)) {
    logger$log(paste0("Calculating cross: ", i, "/", nCrosses))
    p1.id <- which(grepl(crossTable$ind1[i], colnames(g$haplotypes)))
    p2.id <- which(grepl(crossTable$ind2[i], colnames(g$haplotypes)))

    blupExp <- calcProgenyBlupExpected(SNPcoord,
                                       g$haplotypes,
                                       p1.id,
                                       p2.id,
                                       markerEffects)
    blupCovar <- calcProgenyBlupCovariance(SNPcoord = SNPcoord,
                                           r = r,
                                           haplo = g$haplotypes,
                                           p1.id = p1.id,
                                           p2.id = p2.id,
                                           markerEffects = markerEffects,
                                           blupExpectedValues = blupExp)
    blupVar <- diag(blupCovar)

    blupVarExp[[i]]$ind1 <- crossTable$ind1[i]
    blupVarExp[[i]]$ind2 <- crossTable$ind2[i]
    for (trait in names(blupExp)) {
      blupVarExp[[i]]$blup_exp[[trait]] <- blupExp[[trait]]$sum
      blupVarExp[[i]]$blup_var[[trait]] <- blupVar[trait]
    }
    blupVarExp[[i]]$cov <- as.data.frame(blupCovar)
  }

  logger$log("BLUP variance and expected value calculation for each crosses DONE")
  # save and return results
  if (!is.null(outFile)) {
    logger$log("Save results ...")
    file <- save_blupVarExp_as_json(blupVarExp, file = outFile)
    logger$log("Save results DONE")
  } else {
    file <- NULL
  }

  return(blupVarExp)

}





#' Draw a plot of the progenies BLUPs' expected values with error bars
#'
#' X axis is the crosses, and Y axis the blups. The points are located at the
#' expected value and the error bar length is the standard deviation.
#'
#' @param progEstimFile path of the progeny BLUP estimation file generated by
#' r-geno-tools-engine containing the blup estimations of the progenies of some
#' crosses (`json` file).
#' @param progEstimUrl URL of a progeny BLUP estimation file
#' @param sorting method to sort the individuals (X axis) can be:
#'   - "asc": sort the BLUP expected value in ascending order (from left to right)
#'   - "dec": sort the BLUP expected value in decreasing order (from left to right)
#'   - any other value will sort the individuals in alphabetical order (from left to right)
#' @param outFile outFile path of the file containing the plot. If `NULL`, the
#' output will not be written in any file. By default write in an temporary
#' file.
#' @param y_axisName The Y axis name (default = "genetic values")
#' @param errorBarInterval length of XX\% interval of interest represented by the error bars (default=0.95)
#'
#' @return plotly graph
draw_progBlupsPlot <- function(progEstimFile = NULL,
                               progEstimUrl = NULL,
                               errorBarInterval= 0.95,
                               y_axisName = "Genetic values",
                               sorting = 'alpha',
                               trait = NULL,
                               outFile = tempfile(fileext = ".html")) {
  logger <- Logger$new("r-draw_progBlupsPlot()")

  check_sorting(sorting)

  check_outFile(outFile, accept_null = TRUE)
  ext <- tools::file_ext(outFile)
  if (ext == "") {
    if (interactive) {
      outFile <- paste0(outFile, '.html')
    } else {
      outFile <- paste0(outFile, '.png')
    }
  }


  logger$log("Get data ...")
  if (!is.null(progEstimFile) &&  is.null(progEstimUrl)) {
    progBlup <- readProgBlupEstim(progEstimFile)
  } else if (!is.null(progEstimUrl) && is.null(progEstimFile)) {
    progBlup <- downloadProgBlupEstim(progEstimUrl)
  } else {
    engineError("Either `progEstimFile` or `progEstimUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "progeny BLUP estimations"))
  }
  logger$log("Get data DONE")


  available_traits <- sapply(progBlup, function(blupDta_cross) {
    names(blupDta_cross$blup_exp)
  })
  if (!is.null(trait)) {
    check_trait(trait, available_traits)
  } else {
    trait <- available_traits[1]
  }


  logger$log("Draw progenies' blup plot ...")
  p <- plotBlup_1trait(progBlup, sorting = sorting,
                       y_axisName = y_axisName,
                       errorBarInterval = errorBarInterval,
                       trait = trait)
  logger$log("Draw progenies' blup plot DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
      save_plotly(p, outFile)
    logger$log("Save results DONE")
  }
  p
}

#' Draw a plotly graph of blups data for 2 traits
#'
#' The points are located at the expected value and the ellipses
#' size represent the `confidenceLevel` prediction interval.
#'
#' @param progEstimFile path of the progeny BLUP estimation file generated by
#' r-geno-tools-engine containing the blup estimations of the progenies of some
#' crosses (`json` file).
#' @param progEstimUrl URL of a progeny BLUP estimation file
#' @param x_trait name of the trait to show on the x axis
#' @param y_trait name of the trait to show on the y axis
#' @param confidenceLevel level of the prediction ellipses (default 0.95, ie 95%
#' ellypses)
#' @param x_suffix suffix to add to the x axis's name
#' @param y_suffix suffix to add to the y axis's name
#' @param ellipses_npoints number of points used to draw the ellipses (default 100)
#' @param outFile outFile path of the file containing the plot. If `NULL`, the
#' output will not be written in any file. By default write in an temporary
#' file.
#' @return plotly graph
draw_progBlupsPlot_2traits <- function(progEstimFile = NULL,
                                       progEstimUrl = NULL,
                                       x_trait,
                                       y_trait,
                                       confidenceLevel = 0.95,
                                       x_suffix = "",
                                       y_suffix = "",
                                       ellipses_npoints = 100,
                                       outFile = tempfile(fileext = ".html")) {
  logger <- Logger$new("r-draw_progBlupsPlot_2traits()")

  check_confidenceLevel(confidenceLevel)
  check_suffix(x_suffix)
  check_suffix(y_suffix)

  check_outFile(outFile, accept_null = TRUE)
  if (!is.null(outFile)) {
    ext <- tools::file_ext(outFile)
    if (ext == "") {
      if (interactive) {
        outFile <- paste0(outFile, '.html')
      } else {
        outFile <- paste0(outFile, '.png')
      }
    }
  }

  logger$log("Get data ...")
  if (!is.null(progEstimFile) && is.null(progEstimUrl)) {
    progBlup <- readProgBlupEstim(progEstimFile)
  } else if (!is.null(progEstimUrl) && is.null(progEstimFile)) {
    progBlup <- downloadProgBlupEstim(progEstimUrl)
  } else {
    engineError("Either `progEstimFile` or `progEstimUrl` should be provided",
                extra = list(code = errorCode("INPUT_FILE_NOT_PROVIDED"),
                             input_file = "progeny BLUP estimations"))
  }
  logger$log("Get data DONE")

  available_traits <- sapply(progBlup, function(blupDta_cross) {
    names(blupDta_cross$blup_exp)
  })
  check_trait(x_trait, available_traits)
  check_trait(y_trait, available_traits)

  logger$log("Draw progenies' blup plot ...")
  p <- plotBlup_2traits(blupDta = progBlup,
                        x_trait = x_trait,
                        y_trait = y_trait,
                        confidenceLevel = confidenceLevel,
                        x_suffix = x_suffix,
                        y_suffix = y_suffix,
                        ellipses_npoints = ellipses_npoints)
  logger$log("Draw progenies' blup plot DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
      save_plotly(p, outFile)
    logger$log("Save results DONE")
  }
  p
}
