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
  logger <- logger$new("r-draw_manhattanPlot()")

  logger$log("Check outFile ...")
  if (!is.null(outFile)) {
    if (length(outFile) != 1) {
      stop("Error: `outFile` must be of length 1.")
    }
    if (missing(outFile)) {
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
calc_pedRelMAt <- function(pedFile = NULL,
                           pedUrl = NULL,
                           unknown_string = "",
                           header = TRUE,
                           outFile = tempfile(fileext = ".csv"),
                           outFormat = tools::file_ext(outFile)) {
  logger <- logger$new("r-calc_pedRelMAt()")

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
    stop("Error: either `pedFile` or `pedUrl` should be provided")
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
calc_genoRelMAt <- function(genoFile = NULL,
                            genoUrl = NULL,
                            outFile = tempfile(fileext = ".csv"),
                            outFormat = tools::file_ext(outFile)) {
  logger <- logger$new("r-calc_genoRelMAt()")

  logger$log("Get data ...")
  if (!is.null(genoFile) &&  is.null(genoUrl)) {
    geno <- readGenoData(genoFile)
  } else if (!is.null(genoUrl) && is.null(genoFile)) {
    geno <- downloadGenoData(genoUrl)
  } else {
    stop("Error: either genoFile or genoUrl should be provided")
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

  logger <- logger$new("r-calc_combinedRelMat()")

  logger$log("Get data ...")
  if (!is.null(pedRelMatFile) &&  is.null(pedRelMatUrl)) {
    ped_rm <- readRelMat(pedRelMatFile)
  } else if (!is.null(pedRelMatUrl) && is.null(pedRelMatFile)) {
    ped_rm <- downloadRelMat(pedRelMatUrl)
  } else {
    stop("Error: either pedRelMatFile or pedRelMatUrl should be provided")
  }
  if (!is.null(genoRelMatFile) &&  is.null(genoRelMatUrl)) {
    geno_rm <- readRelMat(genoRelMatFile)
  } else if (!is.null(genoRelMatUrl) && is.null(genoRelMatFile)) {
    geno_rm <- downloadRelMat(genoRelMatUrl)
  } else {
    stop("Error: either genoRelMatFile or genoRelMatUrl should be provided")
  }
  logger$log("Get data DONE")

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
  logger <- logger$new("r-draw_relHeatmap()")

  logger$log("Check outFile ...")
  if (!is.null(outFile)) {
    if (length(outFile) != 1) {
      stop("Error: `outFile` must be of length 1.")
    }
    if (missing(outFile)) {
      if (interactive) {
        outFile <- paste0(outFile, '.html')
      } else {
        outFile <- paste0(outFile, '.png')
      }
    }
  }
  logger$log("Check outFile DONE")

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
    stop("Error: either `relMatFile` or `relMatFile` should be provided")
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
      htmlwidgets::saveWidget(plotly::partial_bundle(p),
                              outFile, selfcontained = TRUE)
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
  logger <- logger$new("r-draw_pedNetwork()")

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
    stop("Error: either `pedFile` or `pedUrl` should be provided")
  }
  logger$log("Get data DONE")

  logger$log("Draw pedigree interactive network ...")
  p <- pedNetwork(ped)
  logger$log("Draw pedigree interactive network DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
    htmlwidgets::saveWidget(p,
                            outFile, selfcontained = TRUE)
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
  logger <- logger$new("r-crossingSimulation()")

  # load input data
  logger$log("Get data ...")
  if (!is.null(genoFile) &&  is.null(genoUrl)) {
    g <- readPhasedGeno(genoFile)
  } else if (is.null(genoFile) && !is.null(genoUrl)) {
    g <- downloadPhasedGeno(genoUrl)
  } else {
    stop("Error: either `genoFile` or `genoUrl` should be provided")
  }

  if (!is.null(SNPcoordFile) &&  is.null(SNPcoordUrl)) {
    SNPcoord <- readSNPcoord(SNPcoordFile)
  } else if (is.null(SNPcoordFile) && !is.null(SNPcoordUrl)) {
    SNPcoord <- downloadSNPcoord(SNPcoordUrl)
  } else {
    stop("Error: either `SNPcoordFile` or `SNPcoordUrl` should be provided")
  }

  if (!is.null(crossTableFile) &&  is.null(crossTableUrl)) {
    crossTable <- readCrossTable(crossTableFile)
  } else if (is.null(crossTableFile) && !is.null(crossTableUrl)) {
    crossTable <- downloadCrossTable(crossTableUrl)
  } else {
    stop("Error: either `crossTableFile` or `crossTableUrl` should be provided")
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

  logger$log("Check output file extention ...")
  ext1 <- tools::file_ext(outFile)
  ext2 <- tools::file_ext(tools::file_path_sans_ext(outFile))
  if (ext1 != "gz" || ext2 != "vcf") {
    stop('The output file must end by `.vcf.gz`')
  }
  logger$log("Check output file extention DONE")



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
#' @param markerEffectsFile path of the marker effects file (`csv` file). This `.csv` file should
#' have 2 named columns:
#' - `SNPid`: Marker id
#' - `effects`: effect of the corresponding marker
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
  logger <- logger$new("r-progenyBlupVarExp()")

  # load input data
  logger$log("Get data ...")
  if (!is.null(genoFile) &&  is.null(genoUrl)) {
    g <- readPhasedGeno(genoFile)
  } else if (is.null(genoFile) && !is.null(genoUrl)) {
    g <- downloadPhasedGeno(genoUrl)
  } else {
    stop("Error: either `genoFile` or `genoUrl` should be provided")
  }

  if (!is.null(SNPcoordFile) &&  is.null(SNPcoordUrl)) {
    SNPcoord <- readSNPcoord(SNPcoordFile)
  } else if (is.null(SNPcoordFile) && !is.null(SNPcoordUrl)) {
    SNPcoord <- downloadSNPcoord(SNPcoordUrl)
  } else {
    stop("Error: either `SNPcoordFile` or `SNPcoordUrl` should be provided")
  }

  if (!is.null(crossTableFile) &&  is.null(crossTableUrl)) {
    crossTable <- readCrossTable(crossTableFile)
  } else if (is.null(crossTableFile) && !is.null(crossTableUrl)) {
    crossTable <- downloadCrossTable(crossTableUrl)
  } else {
    stop("Error: either `crossTableFile` or `crossTableUrl` should be provided")
  }

  if (!is.null(markerEffectsFile) &&  is.null(markerEffectsUrl)) {
    markerEffects <- readMarkerEffects(markerEffectsFile)
  } else if (is.null(markerEffectsFile) && !is.null(markerEffectsUrl)) {
    markerEffects <- downloadMarkerEffects(markerEffectsUrl)
  } else {
    stop("Error: either `markerEffectsFile` or `markerEffectsUrl` should be provided")
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
  if (!all(SNPcoord$SNPid %in% row.names(markerEffects))) {
    stop("Missing marker effects for some SNPs of the genetic data.")
  }
  logger$log("Check SNPs' ids consistency between",
             "SNPcoordinate and markerEffects file DONE")


  if (!is.null(outFile)) {
    logger$log("Check output file extention ...")
    ext <- tools::file_ext(outFile)
    if (ext != "json") {
      stop('The output file must end by `.json`')
    }
    logger$log("Check output file extention DONE")
  }


  # calculate recombination rate
  r <- calcRecombRate(SNPcoord)

  # initialize results data.frame
  blupVarExp <- crossTable[, c('ind1', 'ind2')]
  blupVarExp$blup_var <- NA
  blupVarExp$blup_exp <- NA

  # calculation for each couple
  logger$log("BLUP variance and expected value calculation for each crosses ...")
  iter <- 1
  nCrosses <- nrow(crossTable)
  for (i in seq(nCrosses)) {
    logger$log(paste0("Calculating cross: ", i, "/", nCrosses))
    p1.id <- which(grepl(crossTable$ind1[i], colnames(g$haplotypes)))
    p2.id <- which(grepl(crossTable$ind2[i], colnames(g$haplotypes)))

    geneticCovar <- calcProgenyGenetCovar(SNPcoord, r, g$haplotypes, p1.id, p2.id)
    blupVar <- calcProgenyBlupVariance(SNPcoord, markerEffects, geneticCovar)
    blupExp <- calcProgenyBlupExpected(SNPcoord, g$haplotypes,
                                p1.id, p2.id, markerEffects)

    blupVarExp$blup_var[i] <- blupVar
    blupVarExp$blup_exp[i] <- blupExp
  }
  logger$log("BLUP variance and expected value calculation for each crosses DONE")

  # save and return results
  if (!is.null(outFile)) {
    logger$log("Save results ...")
    file <- save_dataFrame_as_json(blupVarExp, file = outFile)
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
#'
#' @return plotly graph
draw_progBlupsPlot <- function(progEstimFile = NULL,
                               progEstimUrl = NULL,
                               sorting = 'alpha',
                               outFile = tempfile(fileext = ".html")) {
  logger <- logger$new("r-draw_progBlupsPlot()")

  logger$log("Check outFile ...")
  if (!is.null(outFile)) {
    if (length(outFile) != 1) {
      stop("Error: `outFile` must be of length 1.")
    }
    if (missing(outFile)) {
      if (interactive) {
        outFile <- paste0(outFile, '.html')
      } else {
        outFile <- paste0(outFile, '.png')
      }
    }
  }
  logger$log("Check outFile DONE")

  logger$log("Get data ...")
  if (!is.null(progEstimFile) &&  is.null(progEstimUrl)) {
    progBlup <- readProgBlupEstim(progEstimFile)
  } else if (!is.null(progEstimUrl) && is.null(progEstimFile)) {
    progBlup <- downloadProgBlupEstim(progEstimUrl)
  } else {
    stop("Error: either progEstimFile or progEstimUrl should be provided")
  }
  logger$log("Get data DONE")

  logger$log("Draw progenies' blup plot ...")
  p <- plotBlup(progBlup, sorting = sorting)
  logger$log("Draw progenies' blup plot DONE")

  if (!is.null(outFile)) {
    logger$log("Save results ...")
    htmlwidgets::saveWidget(p,
                            outFile, selfcontained = TRUE)
    logger$log("Save results DONE")
  }
  p

}
