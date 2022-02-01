# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 The University of Tokyo
#
# Description:
# Functions of the GWAS API related to data loading


#' Download geno data
#'
#' @param url url of the geno data file (.vcf.gz file)
#'
#' @return `gaston::bed.matrix`
downloadGenoData <- function(url) {
  logger <- logger$new("r-downloadGenoData()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "geno",
                        tmpdir = tempdir(),
                        fileext = ".vcf.gz")
  logger$log("Download geno file ...")

  # download data file:
  logger$log("Download genotypic file ... ")
  download.file(url, localFile, quiet = TRUE)
  logger$log("Download genotypic file DONE")
  logger$log("Read genotypic file ...")
  dta <- readGenoData(localFile)
  logger$log("Read genotypic file DONE")

  logger$log("DONE, return output.")
  dta
}

#' Download phenotypic data
#'
#' @param url url of the phenotypic data file (csv file)
#'
#' @details The individuals' names must be on the first column. No duplication
#' is allowed.
#' @return `data.frame`
downloadPhenoData <- function(url){
  logger <- logger$new("r-downloadPhenoData()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "pheno",
                        tmpdir = tempdir(),
                        fileext = ".csv")

  # download data file:
  logger$log("Download phenotypic file ... ")
  download.file(url, localFile, quiet = TRUE)
  logger$log("Download phenotypic file DONE")
  logger$log("Read phenotypic file ...")
  dta <- readPhenoData(localFile)
  logger$log("Read phenotypic file DONE")

  logger$log("DONE, return output.")
  dta
}



#' Download and prepare data for GWAS analysis
#'
#' @param genoUrl url of the geno data file (.vcf.gz file)
#' @param phenoUrl url of the phenotypic data file (csv file)
#'
#' @return List
downloadData <- function(genoUrl, phenoUrl){
  logger <- logger$new("r-downloadData()")
  logger$log("get geno data ...")
  mDta <- downloadGenoData(genoUrl)
  logger$log("get geno data DONE")

  logger$log("get pheno data ...")
  pDta <- downloadPhenoData(phenoUrl)
  logger$log("get pheno data DONE")

  logger$log("prepare data ...")
  dta <- prepareData(mDta, pDta)
  logger$log("prepare data DONE")

  logger$log("DONE, return output.")

  dta
}


#' Download a gwas reults
#'
#' @param url url of the result data file (json file)
#'
#' @return `data.frame`
downloadGWAS <- function(url){
  logger <- logger$new("r-downloadGWAS()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "downloadedResult",
                        tmpdir = tempdir(),
                        fileext = ".json")
  logger$log("Download result file ... ")
  download.file(url, localFile, quiet = TRUE)

  logger$log("Read result file ... ")
  gwasRes <- readGWAS(localFile)
  logger$log("Read result file DONE ")

  logger$log("DONE, return output.")
  gwasRes
}


#' Download pedigree data
#'
#' @param url url of the result data file (csv file)
#'
#' @resturn List of 2: `data` pedigree data, `graph` "igraph" object of the pedigree graph.
downloadPedData <- function(url) {
  logger <- logger$new("r-downloadPedData()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "downloadedResult",
                        tmpdir = tempdir(),
                        fileext = ".csv")
  logger$log("Download result file ... ")
  download.file(url, localFile, quiet = TRUE)

  logger$log("Read result file ... ")
  ped <- readPedData(localFile)
  logger$log("Read result file DONE ")

  logger$log("DONE, return output.")
  ped
}


#' Read geno data from a file
#'
#' @param file VCF file path (ext `.vcf` or `.vcf.gz`)
#'
#' @return `gaston::bed.matrix`
readGenoData <- function(file) {
  logger <- logger$new("r-readGenoData()")
  logger$log("Check file extention ... ")

  if (!file.exists(file)) {
    stop("Genotypic file do not exists")
  }

  ext <- tools::file_ext(file)
  if (identical(ext, "gz")) {
    ext <- paste0(
      tools::file_ext(tools::file_path_sans_ext(file, compression = FALSE)),
      ".", ext)
  }

  if (!(identical(ext, "vcf.gz") || identical(ext, "vcf"))) {
    stop("Genotypic file should be in VCF format. `.vcf`, or `.vcf.gz`. Provided file extention is ", ext)
  }

  logger$log("Read geno file ... ")

  dta <- gaston::read.vcf(file,
                          verbose = TRUE,
                          convert.chr = FALSE)
  logger$log("Read geno file DONE ")
  logger$log("DONE, return output.")

  dta
}


#' Read phenotypic data file
#'
#' @param file file path
#' @param ind.names [default 1] a single number giving the column of the table
#'  which contains the individuals' names.
#' @param ... Further arguments to be passed to `read.csv`
#'
#' @details Any duplication in the phenotypic file is forbidden.
#' @return `data.frame`
readPhenoData <- function(file, ind.names = 1, ...) {
  logger <- logger$new("r-readPhenoData()")

  # read data file:
  if (!file.exists(file)) {
    stop("Phenotypic file do not exists")
  }
  logger$log("Read phenotypic file ... ")
  dta <- read.csv(file, ...)
  logger$log("Read phenotypic file DONE ")

  logger$log("Check individuals unicity ...")
  inds <- dta[, ind.names]
  if (any(duplicated(inds))) {
    stop("Duplicated individuals found in the phenotypic file. ",
         "Individuals must apprear only once in the phenotypic file.")
  }
  logger$log("Check individuals unicity DONE")

  logger$log("Set pheno data's row names ...")
  row.names(dta) <- dta[, ind.names]
  dta <- dta[, -ind.names]
  logger$log("Set pheno data's row names DONE")

  logger$log("DONE, return output.")
  dta
}


#' Read and prepare data for GWAS result
#'
#' @param genoFile path of the geno data file (`.vcf` or `.vcf.gz` file)
#' @param phenoFile path of the phenotypic data file (`csv` file)
#'
#' @return List
readData <- function(genoFile, phenoFile){
  logger <- logger$new("r-readData()")
  logger$log("get geno data ...")
  mDta <- readGenoData(genoFile)
  logger$log("get geno data DONE")

  logger$log("get pheno data ...")
  pDta <- readPhenoData(phenoFile)
  logger$log("get pheno data DONE")

  logger$log("prepare data ...")
  dta <- prepareData(mDta, pDta)
  logger$log("prepare data DONE")

  logger$log("DONE, return output.")

  dta
}


#' Read and prepare pedigree data
#'
#' @param file path of the pedigree data file (`csv` file).
#' @param unknown_string [default: ""] a character vector of strings which are to
#' be interpreted as "unknown parent". By default: missing value in the file.
#' If this information is wrong,
#' @param header [default: TRUE] a logical value indicating whether the file
#' contains the names of the variables as its first line. The default value is
#' TRUE. In any cases, the column 1 will be interpreted as the individual id,
#' column 2 as the first parent, column 3 as the second parent.
#'
#' @details
#' We consider here only allo-fecundation or auto-fecundation. For
#' auto-fecundation, use the parental individual id in both column 2 and 3.
#' Doubles haploids can not be interpreted, please avoid them in the file.
#'
#' Please be sure that all individuals id in columns 2 and 3 are defined in the
#' column 1. If columns 2 and/or 3 contain id of individuals
#' that are not in the first column, a warning will be raised and these
#' individuals will be added to the pedigree with unknown
#' parents as founder individuals.
#'
#'
#' @resturn List of 2: `data` pedigree data, `graph` "igraph" object of the pedigree graph.
readPedData <- function(file, unknown_string = "", header = TRUE) {
  logger <- logger$new('r-readPedData')

  logger$log('Read pedigree file ...')
  if (!file.exists(file)) {
    stop("pedigree file do not exists")
  }
  if (!identical(tools::file_ext(file), 'csv')) {
    stop('Pedigree file shoucd be a `.csv` file.')
  }
  ped <- read.csv(file,
                  na.strings = unknown_string,
                  header = header,
                  stringsAsFactors = FALSE,
                  comment.char = '#')
  logger$log('Read pedigree file DONE')


  logger$log('Check pedigree file ...')
  # file dimention
  if (ncol(ped) != 3) {
    stop('Pedigree file should have exactly 3 columns', ncol(ped), 'detected.')
  }
  if (nrow(ped) == 0) {
    stop('Pedigree file should have at least one row', ncol(ped), 'detected.')
  }

  # NA in first colunm
  if (any(is.na(ped[,1]))) {
    stop('The first colunm of the pedigree file should not have any unknown ',
         'individual.')
  }

  # duplicated rows
  dupRows <- duplicated(ped)
  if (any(dupRows)) {
    warning(sum(dupRows),
            ' duplicated line(s) found in the pedigree file.',
            ' They will be removed.')
    ped <- ped[!dupRows,]
  }

  # inconsistent entries
  incEnt <- duplicated(ped[, 1]) | duplicated(ped[, 1], fromLast = TRUE)
  if (any(incEnt)) {
    stop(length(unique(ped[incEnt, 1])),
         ' inconsistent pedigree entrie(s) found at lines: ',
         paste0('"', which(incEnt), '"',
                collapse = ", ")
         )
  }

  # find missing individual in column 1
  missInd_2 <- na.omit(ped[!ped[, 2] %in% ped[, 1], 2])
  missInd_3 <- na.omit(ped[!ped[, 3] %in% ped[, 1], 3])
  missInd <- unique(c(missInd_2, missInd_3))
  if (length(missInd) != 0) {
    warning(length(missInd), ' unspecified individual(s) found',
            ' in column 2 and/or 3: ',
            ifelse(length(missInd) >= 5,
                   paste(paste0('"', missInd[1:5], '"', collapse = " ,"), "..."),
                   paste(paste0('"', missInd, '"', collapse = " ,"), ".")),
            " These individuals are added in the pedigree with unknown parents."
    )

    missInd_df <- data.frame(missInd, NA, NA,
                             stringsAsFactors = FALSE)
    ped <- rbind(setNames(missInd_df, colnames(ped)),
                 ped,
                 stringsAsFactors = FALSE)

    # Unique founder individual.
    # If there is only one founder corresponding to the only missing
    # individual, `unknown_string`, might have been miss-specified...(It will be
    # the case only if `unknown_string` is miss-specified and all the
    # the founder individuals are defined in the pedigree. This explains the
    # recommendation about specifying all the founders in the pedigree.)
    founders <- ped[is.na(ped[,2]) & is.na(ped[,3]), 1]
    if (length(founders) == 1 && length(missInd) == 1) {
      if (identical(founders, missInd)) {
        warning('This individual is the only founder found in your pedigree. ',
                'Please be sure you have specified the id of unknown parrent ',
                'correctly.')
      }
    }
  }

  # inconsistent genealogy (may be long for big pedigree file)
  # check that individual are not a parent of their parents.
  edges <- c(as.vector(t(na.omit(ped[,c(2,1)]))),
             as.vector(t(na.omit(ped[,c(3,1)]))))
  g <- igraph::make_graph(edges, directed = TRUE)

  cycle <- list()
  for (i in seq_along(igraph::V(g))) {
    i <- igraph::V(g)[i]
    parents <- igraph::neighbors(g, i, mode = "in")
    for (p in parents) {
      paths <- igraph::all_simple_paths(g, i, p, 'out')
      if (length(paths) != 0) {
        cycle <- append(cycle, lapply(paths, function(p){names(sort(p))}))
      }
    }
  }
  cycle <- cycle[!duplicated(cycle)]
  if (length(cycle) != 0) {
    stop(length(cycle), ' inconsistent genealogy detected involving:\n',
         paste('\t - ',
            sapply(cycle, function(c){paste0('"', c, '"', collapse = ", ")}),
            collapse = "\n"
         ))
  }
  logger$log('Check pedigree file DONE')


  logger$log("DONE, return output.")
  # rename colunms
  colnames(ped) <- c('ind', 'parent1', 'parent2')

  return(list(
    data = ped,
    graph = g # TODO check if it is really necessary to return the graph
  ))
}





#' Read GWAS analysis result file (`.json`)
#'
#' @param file path of the json file generated by the function `gwas` containing GWAS result
#'
#' @return `list` of 2 elements `gwas` (data.frame) and `metadata` (list)
readGWAS <- function(file) {
  logger <- logger$new("r-readGWAS()")

  logger$log("Read result file ... ")
  if (!file.exists(file)) {
    stop("GWAS file do not exists")
  }
  gwasRes <- readLines(file)
  logger$log("Read result file DONE ")
  logger$log("Convert Json to data.frame ... ")
  gwasRes <- jsonlite::fromJSON(gwasRes)
  logger$log("Convert Json to data.frame DONE ")
  if (class(gwasRes) != "list" || names(gwasRes) != c("gwas", "metadata") ||
      class(gwasRes$gwas) != "data.frame") {
    stop("Error: provided file do not seems to be generated by r-geno-tools-Engine.")
  }
  logger$log("DONE, return output.")
  gwasRes
}

#' saveGWAS save gwas result in a temporary file
#'
#' @param gwasRes data.frame return by `gwas` function
#' @param metadata list of metadata of the gwas results
#' @param dir  if `filename` is NULL, directory where to save the data,
#' by default it is a temporary directory
#' @param file file path where to save the data. If the file already exists, it
#' will be overwritten. Default NULL
#'
#' @return path of the created filed
saveGWAS <- function(gwasRes, metadata, dir = NULL, file = NULL) {
  logger <- logger$new("r-saveGWAS()")
  if (is.null(file)) {
    if (is.null(dir)) {
       dir <- tempdir()
    }
    logger$log('Check dir ...')
    if (!dir.exists(dir)) {
      logger$log('Error: "dir" directory should exists')
      stop('Error: "dir" directory should exists')
    }
    logger$log('Check dir DONE')
    file <- tempfile(fileext = ".json", tmpdir = dir)

  } else {
    logger$log('Check file ...')
    if (length(file) != 1) {
      logger$log('Error: only one file name should be provided')
      stop('Error: only one file name should be provided')
    }
    if (file.exists(file)) {
      logger$log('Warning: "file" directory already exists. This file will be overwritten.')
    } else {
      file.create(file)
    }
    logger$log('Check file DONE')
  }

  gwasList <- list(gwas = gwasRes,
                   metadata = metadata)
  writeLines(jsonlite::toJSON(gwasList,
                              complex = "list",
                              pretty = T,
                              digits = NA,
                              na = 'string'),
             con = file)
  return(file)
}

#' Filter individuals and remove monomorphic markers
#'
#' @param gDta output of `downloadGenoData` or `readGenoData` functions
#' @param pDta output of `downloadPhenoData` or `readPhenoData` functions
#'
#' @details The function remove the monomorphic markers and
#' @return List of 2 elements: `genoData` (a bed matrix), `phenoData` (a data.frame)
prepareData <- function(gDta, pDta) {
  logger <- logger$new("r-prepareData()")
  # Remove from geno data individuals that are not in phenotypic data-set
  logger$log("Remove from geno data individuals that are not in phenotypic data-set ...")
  gDta <- gaston::select.inds(gDta, id %in% rownames(pDta))
  logger$log("Remove from geno data individuals that are not in phenotypic data-set DONE")


  # reorder phenotypic data with id in bed matrix
  logger$log("reorder matrix ...")
  pDta <- pDta[gDta@ped$id,]
  logger$log("reorder matrix DONE")


  # remove monomorphic markers
  logger$log("remove monomorphic markers ...")
  gDta <- gaston::select.snps(gDta, maf > 0)
  logger$log("remove monomorphic markers DONE")

  logger$log("DONE, return output.")

  list(genoData = gDta, phenoData = pDta)
}
