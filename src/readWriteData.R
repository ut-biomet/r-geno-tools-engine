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
  tryCatch({
    download.file(url, localFile)
  }, error = function(err) {
    logger$log("Download geno file FAILED.\n\tURL is:", url,
               "\nError message is:\n\t", err$message,
               redis = TRUE,
               status = "FAILED",
               action_type = "GET_MARKER_DTA")
    return(NULL)
  })

  dta <- readGenoData(localFile)

  dta
}

#' Download phenotypic data
#'
#' @param url url of the phenotypic data file (csv file)
#'
#' @return `data.frame`
downloadPhenoData <- function(url){
  logger <- logger$new("r-getPhenoData()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "pheno",
                        tmpdir = tempdir(),
                        fileext = ".csv")

  # download data file:
  logger$log("Download phenotypic file ... ")
  tryCatch({
    download.file(url, localFile)
  }, error = function(err) {
    logger$log("Download geno file FAILED.\n\tURL is:", url,
               "\nError message is:\n\t", err$message,
               redis = TRUE,
               status = "FAILED",
               action_type = "GET_PHENO_DTA")
    return(NULL)
  })

  dta <- readPhenoData(localFile)

  logger$log("DONE, return output.")
  dta
}



#' Download data for GWAS model
#'
#' @param genoUrl
#' @param phenoUrl
#'
#' @return List
downloadData <- function(genoUrl, phenoUrl){
  logger <- logger$new("r-loadData()")
  logger$log("get geno data ...")
  mDta <- downloadGenoData(genoUrl)
  logger$log("get geno data DONE")

  logger$log("get pheno data ...")
  pDta <- downloadPhenoData(phenoUrl)
  logger$log("get pheno data DONE")


  # Remove from geno data individuals that are not in phenotypic data-set
  logger$log("Remove from geno data individuals that are not in phenotypic data-set ...")
  mDta <- select.inds(mDta, id %in% rownames(pDta))
  logger$log("Remove from geno data individuals that are not in phenotypic data-set DONE")


  # reorder phenotypic data with id in bed matrix
  logger$log("reorder matrix ...")
  pDta <- pDta[mDta@ped$id,]
  logger$log("reorder matrix DONE")


  # remove monomorphic geno
  logger$log("remove monomorphic geno ...")
  mDta <- select.snps(mDta, maf > 0)
  logger$log("remove monomorphic geno DONE")


  # calculate genetic relational matrix
  logger$log("calculate genetic relatinoal matrix ...")
  grm <- GRM(mDta)
  logger$log("calculate genetic relatinoal matrix DONE")

  logger$log("DONE, return output.")

  list(genoData = mDta, phenoData = pDta, grMatrix = grm)
}


#' Download a gwas reults
#'
#' @param url url of the model data file (json file)
#'
#' @return `data.frame`
downloadGWAS <- function(url){
  logger <- logger$new("r-loadModel()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "downloadedModel",
                        tmpdir = tempdir(),
                        fileext = ".json")
  logger$log("Download model file ... ")
  download.file(url, localFile)

  logger$log("Read model file ... ")
  gwasRes <- readLines(localFile)
  logger$log("Read model file DONE ")
  logger$log("Convert Json to data.frame ... ")
  gwasRes <- fromJSON(gwasRes)
  gwasRes$gwas <- as.data.frame(gwasRes$gwas)
  logger$log("Convert Json to data.frame DONE ")

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


#' Title
#'
#' @param file file path
#' @param row.names [default 1] a single number giving the column of the table which contains the row names
#' @param ... Further arguments to be passed to `read.csv`
#'
#' @return `data.frame`
readPhenoData <- function(file, row.names = 1, ...) {
  logger <- logger$new("r-readPhenoData()")
  # read data file:
  logger$log("Read phenotypic file ... ")
  dta <- read.csv(file,
                  row.names = row.names,
                  ...)
  logger$log("Read phenotypic file DONE ")
  logger$log("DONE, return output.")
  dta
}


#' Download data for GWAS model
#'
#' @param genoFile
#' @param phenoFile
#'
#' @return List
readData <- function(genoFile, phenoFile){
  logger <- logger$new("r-readData()")
  logger$log("get geno data ...")
  mDta <- getGenoData(genoFile)
  logger$log("get geno data DONE")

  logger$log("get pheno data ...")
  pDta <- readPhenoData(phenoFile)
  logger$log("get pheno data DONE")


  # Remove from geno data individuals that are not in phenotypic data-set
  logger$log("Remove from geno data individuals that are not in phenotypic data-set ...")
  mDta <- select.inds(mDta, id %in% rownames(pDta))
  logger$log("Remove from geno data individuals that are not in phenotypic data-set DONE")


  # reorder phenotypic data with id in bed matrix
  logger$log("reorder matrix ...")
  pDta <- pDta[mDta@ped$id,]
  logger$log("reorder matrix DONE")


  # remove monomorphic geno
  logger$log("remove monomorphic geno ...")
  mDta <- select.snps(mDta, maf > 0)
  logger$log("remove monomorphic geno DONE")


  # calculate genetic relational matrix
  logger$log("calculate genetic relatinoal matrix ...")
  grm <- GRM(mDta)
  logger$log("calculate genetic relatinoal matrix DONE")

  logger$log("DONE, return output.")

  list(genoData = mDta, phenoData = pDta, grMatrix = grm)
}


#' Title
#'
#' @param file json file generated by the function `gwas` containing GWAS result
#'
#' @return `data.frame`
readGWAS <- function(file) {
  logger <- logger$new("r-readGWAS()")

  logger$log("Read model file ... ")
  gwasRes <- readLines(file)
  logger$log("Read model file DONE ")
  logger$log("Convert Json to data.frame ... ")
  gwasRes <- jsonlite::fromJSON(gwasRes)
  logger$log("Convert Json to data.frame DONE ")

  logger$log("DONE, return output.")
  gwasRes
}

#' saveGWAS save gwas result in a temporary file
#'
#' @param gwas data.frame return by `gwas` function
#'
#' @return path of the created file
#' @export
#'
#' @examples
saveGWAS <- function(gwas) {
  logger <- logger$new("r-saveGWAS()")

  file <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(gwas, dataframe = "rows", pretty = T),
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
  gDta <- select.inds(gDta, id %in% rownames(pDta))
  logger$log("Remove from geno data individuals that are not in phenotypic data-set DONE")


  # reorder phenotypic data with id in bed matrix
  logger$log("reorder matrix ...")
  pDta <- pDta[gDta@ped$id,]
  logger$log("reorder matrix DONE")


  # remove monomorphic geno
  logger$log("remove monomorphic geno ...")
  gDta <- select.snps(gDta, maf > 0)
  logger$log("remove monomorphic geno DONE")


  # calculate genetic relational matrix
  logger$log("calculate genetic relatinoal matrix ...")
  grm <- GRM(gDta)
  logger$log("calculate genetic relatinoal matrix DONE")

  logger$log("DONE, return output.")

  list(genoData = gDta, phenoData = pDta, grMatrix = grm)
}
