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





#' Read geno data from a file
#'
#' @param file VCF file path (ext `.vcf` or `.vcf.gz`)
#'
#' @return `gaston::bed.matrix`
readGenoData <- function(file) {
  logger <- logger$new("r-readGenoData()")
  logger$log("Check file extention ... ")

  if (!file.exists(file)) {
    stop("File do not exists") 
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
#' @param row.names [default 1] a single number giving the column of the table which contains the row names
#' @param ... Further arguments to be passed to `read.csv`
#'
#' @return `data.frame`
readPhenoData <- function(file, row.names = 1, ...) {
  logger <- logger$new("r-readPhenoData()")
  # read data file:
  if (!file.exists(file)) {
    stop("File do not exists") 
  }
  logger$log("Read phenotypic file ... ")
  dta <- read.csv(file,
                  row.names = row.names,
                  ...)
  logger$log("Read phenotypic file DONE ")
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


#' Read GWAS analysis result file (`.json`)
#'
#' @param file path of the json file generated by the function `gwas` containing GWAS result
#'
#' @return `list` of 2 elements `gwas` (data.frame) and `metadata` (list)
readGWAS <- function(file) {
  logger <- logger$new("r-readGWAS()")

  logger$log("Read result file ... ")
  if (!file.exists(file)) {
    stop("File do not exists") 
  }
  gwasRes <- readLines(file)
  logger$log("Read result file DONE ")
  logger$log("Convert Json to data.frame ... ")
  gwasRes <- jsonlite::fromJSON(gwasRes)
  logger$log("Convert Json to data.frame DONE ")
  
  if (class(gwasRes) != "list" || names(gwasRes) != c("gwas", "metadata") ||
      class(gwasRes$gwas) != "data.frame") {
    stop("Error: provided file do not seems to be generated by GWAS-Engine.")
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

#' Filter individuals and calculate genetic relatinoal matrix
#'
#' @param gDta output of `downloadGenoData` or `readGenoData` functions
#' @param pDta output of `downloadPhenoData` or `readPhenoData` functions
#'
#' @return List of 3 elements: `genoData`, `phenoData`, `grMatrix`
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


  # remove monomorphic geno
  logger$log("remove monomorphic geno ...")
  gDta <- gaston::select.snps(gDta, maf > 0)
  logger$log("remove monomorphic geno DONE")


  # calculate genetic relational matrix
  logger$log("calculate genetic relatinoal matrix ...")
  grm <- gaston::GRM(gDta)
  logger$log("calculate genetic relatinoal matrix DONE")

  logger$log("DONE, return output.")

  list(genoData = gDta, phenoData = pDta, grMatrix = grm)
}
