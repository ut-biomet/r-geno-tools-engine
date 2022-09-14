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
  logger <- Logger$new("r-downloadGenoData()")
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

#' Download phased geno data
#'
#' @param url url of the geno data file (.vcf.gz file)
#'
#' @return list of 2: `haplotypes` a matrix of the individuals haplotypes
#'  and `SNPcoord`, data frame of the SNP coordinates.
downloadPhasedGeno <- function(url) {
  logger <- Logger$new("r-downloadPhasedGeno()")
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
  dta <- readPhasedGeno(localFile)
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
  logger <- Logger$new("r-downloadPhenoData()")
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
  logger <- Logger$new("r-downloadData()")
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
  logger <- Logger$new("r-downloadGWAS()")
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
#' @param url url of the pedigree data file (`csv` file).
#' @param unknown_string [default: ""] a character vector of strings which are
#' to be interpreted as "unknown parent". By default: missing value in the file.
#' @param header [default: TRUE] a logical value indicating whether the file
#' contains the names of the variables as its first line. The default value is
#' TRUE. In any cases, the column 1 will be interpreted as the individual id,
#' column 2 as the first parent, column 3 as the second parent.
#'
#' @return List of 2: `data` pedigree data, `graph` "igraph" object of the pedigree graph.
downloadPedData <- function(url, unknown_string = "", header = TRUE) {
  logger <- Logger$new("r-downloadPedData()")
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



#' Download relationship matrix
#'
#' @param url url of the result data file (csv or json file)
#' @param format format of the input file. Either "csv" or "json" (optional, by
#' default it will use "json").
#' @return matrix
downloadRelMat <- function(url, format = tools::file_ext(url)) {
  logger <- Logger$new("r-downloadRelMat()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "downloadedResult",
                        tmpdir = tempdir(),
                        fileext = paste0(".", format))
  logger$log("Download result file ... ")
  download.file(url, localFile, quiet = TRUE)

  logger$log("Read result file ... ")
  relMat <- readRelMat(localFile, format = format)
  logger$log("Read result file DONE ")

  logger$log("DONE, return output.")
  relMat
}


#' Download crossing table
#'
#' @param url url of the crossing table data file (`csv` file of 2 or 3 columns).
#' @param header [default: TRUE] a logical value indicating whether the file
#' contains the names of the variables as its first line. The default value is
#' TRUE. In any cases, the column 1 and 2 will be interpreted as the parents id.
#' The optional third column will be interpreted as the offspring base name.
#'
#' @return data.frame with the crossing table information.
downloadCrossTable <- function(url, header = TRUE) {
  logger <- Logger$new("r-downloadCrossTable()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "downloadedResult",
                        tmpdir = tempdir(),
                        fileext = ".csv")
  logger$log("Download result file ... ")
  download.file(url, localFile, quiet = TRUE)

  logger$log("Read result file ... ")
  crossTable <- readCrossTable(localFile, header)
  logger$log("Read result file DONE ")

  logger$log("DONE, return output.")
  crossTable
}


#' download SNP coordinates `.csv` file
#'
#' @param url url of the SNPs coordinates file (`csv` file). This `.csv` file can have 4 named columns:
#' - `chr`: Chromosome holding the SNP
#' - `physPos`: SNP physical position on the chromosome
#' - `linkMapPos`: SNP linkage map position on the chromosome in Morgan
#' - `SNPid`: SNP's IDs
#'
#' If `SNPid` columns is missing or have missing values, the SNPid will be
#' automatically imputed using the convention `chr@physPos` therefore columns
#' `chr` and `physPos` should not have any missing values
#'
#' @return data.frame of 4 columns: 'chr', 'physPos', 'linkMapPos', 'SNPid'
downloadSNPcoord <- function(url) {
  logger <- Logger$new("r-downloadSNPcoord()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "downloadedResult",
                        tmpdir = tempdir(),
                        fileext = ".csv")
  logger$log("Download result file ... ")
  download.file(url, localFile, quiet = TRUE)

  logger$log("Read result file ... ")
  snpCoord <- readSNPcoord(localFile)
  logger$log("Read result file DONE ")

  logger$log("DONE, return output.")
  snpCoord
}


#' Download marker effects file
#'
#' @param url url of the marker effects file (`csv` file). This `.csv` file should
#' have 2 named columns:
#' - `SNPid`: Marker id
#' - `effects`: effect of the corresponding marker
#'
#' @return data.frame of 1 columns named `effects` with the marker ids as
#' row names.
downloadMarkerEffects <- function(url) {
  logger <- Logger$new("r-downloadMarkerEffects()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "markerEffects",
                        tmpdir = tempdir(),
                        fileext = ".csv")
  logger$log("Download marker effects file ... ")
  download.file(url, localFile, quiet = TRUE)

  logger$log("Read marker effects file ... ")
  markerEffects <- readMarkerEffects(localFile)
  logger$log("Read marker effects file DONE ")

  logger$log("DONE, return output.")
  markerEffects
}

#' Download progeny BLUP estimation file
#'
#' @param url url of the progeny BLUP estimation file generated by
#' r-geno-tools-engine containing the blup estimations of the progenies of some
#' crosses (`json` file).
#'
#' @return data.frame of 4 columns named "ind1", "ind2", "blup_var", "blup_exp"
downloadProgBlupEstim <- function(url) {
  logger <- Logger$new("r-downloadProgBlupEstim()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "progBlupEstim",
                        tmpdir = tempdir(),
                        fileext = ".json")
  logger$log("Download progenies' blup estimations file ... ")
  download.file(url, localFile, quiet = TRUE)

  logger$log("Read progenies' blup estimations file ... ")
  projBlups <- readProgBlupEstim(localFile)
  logger$log("Read progenies' blup estimations file DONE ")

  logger$log("DONE, return output.")
  projBlups
}




#' Read geno data from a file
#'
#' @param file VCF file path (ext `.vcf` or `.vcf.gz`)
#'
#' @return `gaston::bed.matrix`
readGenoData <- function(file) {
  logger <- Logger$new("r-readGenoData()")
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
    stop("Genotypic file should be in VCF format. `.vcf`, or `.vcf.gz`. ",
         "Provided file extention is ", ext)
  }

  logger$log("Read geno file ... ")

  dta <- gaston::read.vcf(file,
                          verbose = TRUE,
                          convert.chr = FALSE)

  # impute missing SNP ids
  missingIDlines <- is.na(dta@snps$id)
  dta@snps$id[missingIDlines] <- paste0(
    dta@snps$chr[missingIDlines],
    '@',
    dta@snps$pos[missingIDlines])
  dta@snps$id[is.na(dta@snps$chr[missingIDlines]) | is.na(dta@snps$pos[missingIDlines])] <- NA

  logger$log("Read geno file DONE ")
  logger$log("DONE, return output.")

  dta
}


#' Read phased genetic data from a file
#'
#' @param file phased VCF file path (ext `.vcf` or `.vcf.gz`)
#'
#' @return list of 2: `haplotypes` a matrix of the individuals haplotypes
#'  and `SNPcoord`, data frame of the SNP coordinates.
readPhasedGeno <- function(file) {
  logger <- Logger$new("r-readPhasedGeno()")
  logger$log("Check file extention ... ")

  if (!file.exists(file)) {
    stop("Phased genotypic file do not exists")
  }

  ext <- tools::file_ext(file)
  if (identical(ext, "gz")) {
    ext <- paste0(
      tools::file_ext(tools::file_path_sans_ext(file, compression = FALSE)),
      ".", ext)
  }

  if (!(identical(ext, "vcf.gz") || identical(ext, "vcf"))) {
    stop("Phased genotypic file should be in VCF format. `.vcf`, or `.vcf.gz`. Provided file extention is ", ext)
  }

  logger$log("Read geno file ... ")
  vcf <- vcfR::read.vcfR(file, verbose = FALSE)
  logger$log("Read phased geno file DONE")

  colnames(vcf@fix)

  # impute missing SNP ids
  missingIDlines <- is.na(vcf@fix[,'ID'])
  vcf@fix[missingIDlines, 'ID'] <- paste0(
    vcf@fix[missingIDlines, 'CHROM'],
    '@',
    vcf@fix[missingIDlines, 'POS'])
  vcf@fix[is.na(vcf@fix[missingIDlines,'CHROM']) | is.na(vcf@fix[missingIDlines,'POS'])
          ,'ID'] <- NA


  # get SNP information
  logger$log("Extract SNP information...")
  fixVcf <- as.data.frame(vcfR::getFIX(vcf), stringsAsFactors = FALSE)
  SNPcoord <- fixVcf[,c("CHROM", "POS", "ID")]
  colnames(SNPcoord) <- c("chr", "physPos", "SNPid")
  SNPcoord$physPos <- as.integer(SNPcoord$physPos)
  logger$log("Extract SNP information DONE")




  # Check if genotypes are phased
  logger$log("Check pahsing ...")
  if (all(vcf@gt[,"FORMAT"] == "GT")) {
    # quick return if FORMAT == GT
    gt <- vcf@gt
    gt <- gt[, colnames(gt) != "FORMAT"]
    row.names(gt) <- vcf@fix[,"ID"]
  } else {
    gt <- vcfR::extract.gt(vcf, element = 'GT')
  }
  if (!all(grepl("|", gt, fixed = TRUE))) {
    errMsg <- "VCF file should be phased for all variant and all individuals, (`|` separator for GT field)."
    logger$log("ERROR:", errMsg)
    stop(errMsg)
  }
  logger$log("Check pahsing DONE")


  # extract haplotypes
  logger$log("Extract haplotypes...")
  indNames <- colnames(gt)
  markersNames <- rownames(gt)

  # paste all values in one vector (make strsplit faster)
  gt <- paste(gt, collapse = "|")
  gt <- strsplit(gt, split = "|", fixed = TRUE)
  gt <- unlist(gt)
  gt <- as.integer(gt)
  # values of gt are mixed so we need to reorder:
  gt <- c(gt[seq(1, length(gt), 2)], gt[seq(2, length(gt), 2)])
  gt <- matrix(gt,
               nrow = length(markersNames),
               byrow = FALSE)
  rownames(gt) <- markersNames
  colnames(gt) <- paste(rep(indNames, 2),
                        rep(1:2, each = length(indNames)),
                        sep = "_")
  logger$log("Extract haplotypes DONE")

  logger$log("DONE, return output.")

  list(haplotypes = gt,
       SNPcoord = SNPcoord)

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
  logger <- Logger$new("r-readPhenoData()")

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
  dta <- dta[, -ind.names, drop = FALSE]
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
  logger <- Logger$new("r-readData()")
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
#' @param unknown_string [default: ""] a character vector of strings which are
#' to be interpreted as "unknown parent". By default: missing value in the file.
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
#' @return List of 2: `data` pedigree data, `graph` "igraph" object of the pedigree graph.
readPedData <- function(file, unknown_string = "", header = TRUE) {
  logger <- Logger$new('r-readPedData')

  logger$log('Read pedigree file ...')
  if (!file.exists(file)) {
    stop("pedigree file do not exists")
  }
  if (!identical(tools::file_ext(file), 'csv')) {
    stop('Pedigree file should be a `.csv` file.')
  }
  ped <- read.csv(file,
                  na.strings = unknown_string,
                  header = header,
                  stringsAsFactors = FALSE,
                  comment.char = '#')
  ped <- data.frame(lapply(ped, as.character))
  logger$log('Read pedigree file DONE')


  logger$log('Check pedigree file ...')
  # file dimention
  if (ncol(ped) != 3) {
    stop('Pedigree file should have exactly 3 columns', ncol(ped), 'detected.')
  }
  if (nrow(ped) == 0) {
    stop('Pedigree file should have at least one row')
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


#' Read a relationship matrix file
#'
#' @param file path of the file generated by the function `saveRelMat`
#' containing relationship matrix
#' @param format format of the input file. Either "csv" or "json" (optional, by
#' default it will use the `file` extension).
#'
#' @details The metadata of the file are not kept.
#'
#' @return matrix
readRelMat <- function(file, format = tools::file_ext(file)) {
  logger <- Logger$new("r-readRelMat()")


  logger$log('Check file format ...')
  formatErr <- FALSE
  if (length(format) != 1) {
    formatErr <- TRUE
  } else if (!format %in% c('csv', 'json')) {
    formatErr <- TRUE
  }
  if (formatErr) {
    logger$log('Error: File format misspecified. It should be either "csv",
               or "json".')
    stop('Error: File format misspecified. It should be either "csv",
               or "json".')
  }
  logger$log('Check file format DONE')


  if (!file.exists(file)) {
    stop("Relationship matrix file do not exists")
  }

  if (format == 'csv') {
    logger$log("Read relationship matrix `csv` file ... ")
    relMat <- read.csv(file = file,
                       header = TRUE,
                       check.names = FALSE,
                       row.names = 1,
                       comment.char = "#")
    relMat <- as.matrix(relMat)
    logger$log("Read relationship matrix `csv` file DONE")
  }
  if (format == 'json') {
    logger$log("Read relationship matrix `json` file ... ")
    relMatRaw <- jsonlite::fromJSON(readLines(file))
    relMat <- as.matrix(relMatRaw$relMat)
    logger$log("Read relationship matrix `json` file DONE")
  }

  logger$log("Check loaded relationship matrix ...")
  checkRelMat(relMat)
  logger$log("Check loaded relationship matrix DONE")

  logger$log("DONE, return output.")
  relMat
}



#' Read crossing table
#'
#' @param file path of the crossing table data file (`csv` file of 2 or 3 columns).
#' @param header [default: TRUE] a logical value indicating whether the file
#' contains the names of the variables as its first line. The default value is
#' TRUE. In any cases, the column 1 and 2 will be interpreted as the parents id.
#' The optional third column will be interpreted as the offspring base name.
#'
#' @return data.frame with the crossing table information.
readCrossTable <- function(file, header = TRUE) {
  logger <- Logger$new('r-readCrossTable')

  logger$log('Read crossing table file ...')
  if (!file.exists(file)) {
    stop("crossing table file do not exists")
  }
  if (!identical(tools::file_ext(file), 'csv')) {
    stop('crossing table file should be a `.csv` file.')
  }
  crossTable <- read.csv(file,
                  header = header,
                  stringsAsFactors = FALSE)
  logger$log('Read crossing table file DONE')


  logger$log('Check crossing table file ...')
  # file dimension
  if (!ncol(crossTable) %in% c(2, 3)) {
    stop('crossing table file should have 2 or 3 columns',
         ncol(crossTable), 'detected.')
  }
  if (nrow(crossTable) == 0) {
    stop('crossing table file should have at least one row')
  }

  # missing values
  if (any(is.na(crossTable))) {
    stop('crossing table should not have any missing values')
  }

  # duplicated rows
  dupRows <- duplicated(crossTable[, c(1,2)])
  if (any(dupRows)) {
    warning(sum(dupRows),
            ' duplicated patental pair(s) found in the crossing table file.',
            ' They will be removed.')
    crossTable <- crossTable[!dupRows,]
  }

  # rename columns
  colnames(crossTable)[c(1,2)] <- c('ind1', 'ind2')
  if (ncol(crossTable) == 3) {
    colnames(crossTable)[3] <- 'names'
  }

  # Get simulated individuals names if not specified
  if (is.null(crossTable$names)) {
    logger$log("Generate simulated individuals names...")
    crossTable$names <- paste0(crossTable$ind1, '_X_', crossTable$ind2)
    logger$log("Generate simulated individuals names DONE")
  }

  return(crossTable)
}



#' Read SNP coordinates `.csv` file
#'
#' @param file path of the SNPs coordinates file (`csv` file). This `.csv` file can have 4 named columns:
#' - `chr`: Chromosome holding the SNP
#' - `physPos`: SNP physical position on the chromosome
#' - `linkMapPos`: SNP linkage map position on the chromosome in Morgan
#' - `SNPid`: SNP's IDs
#'
#' If `SNPid` columns is missing or have missing values, the SNPid will be
#' automatically imputed using the convention `chr@physPos` therefore columns
#' `chr` and `physPos` should not have any missing values
#'
#' @return data.frame of 4 columns: 'chr', 'linkMapPos', 'SNPid'
readSNPcoord <- function(file) {
  logger <- Logger$new('r-readSNPcoord()')

  logger$log('Read snps coordinates file ...')
  if (!file.exists(file)) {
    stop("snps coordinates file do not exists")
  }
  if (!identical(tools::file_ext(file), 'csv')) {
    stop('snps coordinates file should be a `.csv` file.')
  }
  SNPcoord <- read.csv(file,
                       header = TRUE,
                       stringsAsFactors = FALSE)
  logger$log('Read snps coordinates file DONE')


  logger$log('Check snps coordinates file ...')
  if (nrow(SNPcoord) == 0) {
    stop('snps coordinates file should have at least one row')
  }

  refColNames <- sort(c('chr', 'physPos', 'SNPid', 'linkMapPos'))
  missingVar <- !refColNames %in% colnames(SNPcoord)
  names(missingVar) <- refColNames

  if (missingVar['linkMapPos']) {
    # linkMapPos is mandatory
    stop('snps coordinates file should have the header `linkMapPos`specifying ',
         'the linkage map position in Morgan.')
  }
  if (missingVar['SNPid'] & (missingVar['chr'] | missingVar['physPos'])) {
    # we need either SNPid or (chr and physPos)
    stop('snps coordinates file should have a header specifying either ',
         '`SNPid` (id of the markers) and/or both `chr` and `physPos` ',
         '(chromosome and physical position).')
  }
  # add missing values for missing variables
  if (any(missingVar)) {
    SNPcoord[[names(which(missingVar))]] <- NA
  }
  SNPcoord <- SNPcoord[, c("chr", "physPos", "SNPid", "linkMapPos")]

  # check missing values
  if (any(is.na(SNPcoord$linkMapPos))) {
    # no missing values for linkMapPos
    stop('snps coordinates should not have any missing values for the linkage ',
         'map positions')
  }
  if (any(is.na(SNPcoord$SNPid) & (is.na(SNPcoord$chr) | is.na(SNPcoord$physPos)))) {
    # no missing values SNPid or (chr and physPos)
    stop('snps coordinates should not have any missing values for either ',
         '`SNPid` (id of the markers) or any of `chr` and `physPos` ',
         '(chromosome and physical position).')
  }


  # impute missing SNPid
  missingIDlines <- is.na(SNPcoord$SNPid)
  SNPcoord$SNPid[missingIDlines] <- paste0(
    SNPcoord$chr[missingIDlines],
    '@',
    SNPcoord$physPos[missingIDlines]
    )


  # remove duplicated rows
  dupRows <- duplicated(SNPcoord)
  if (any(dupRows)) {
    warning(sum(dupRows),
            'duplicated rows found in the snps coordinates file. ',
            'They will be removed.')
    SNPcoord <- SNPcoord[!dupRows,]
  }

  # check unicity of SNPids
  duplicatedIds <- which(duplicated(SNPcoord$SNPid))
  if (length(duplicatedIds) != 0) {
    msg <- paste(
    length(duplicatedIds),
    'duplicated SNPs\' id detected in the SNP coordinate file:',
    paste(SNPcoord$SNPid[duplicatedIds], collapse = ', '))
    stop(msg)
  }
  logger$log('Check snps coordinates file DONE')

  return(SNPcoord)
}


#' Read GWAS analysis result file (`.json`)
#'
#' @param file path of the json file generated by the function `gwas` containing GWAS result
#'
#' @return `list` of 2 elements `gwas` (data.frame) and `metadata` (list)
readGWAS <- function(file) {
  logger <- Logger$new("r-readGWAS()")

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



#' Read marker effects file
#'
#' @param file path of the marker effects file (`csv` file). This `.csv` file should
#' have 2 named columns:
#' - `SNPid`: Marker id
#' - `effects`: effect of the corresponding marker
#'
#' @return data.frame of 1 columns named `effects` with the marker ids as
#' row names.
readMarkerEffects <- function(file) {
  logger <- Logger$new("r-readMarkerEffects()")

  logger$log('Marker effects table file ...')
  if (!file.exists(file)) {
    stop("Marker effects file do not exists")
  }
  if (!identical(tools::file_ext(file), 'csv')) {
    stop('Marker effects file should be a `.csv` file.')
  }
  markerEffects <- read.csv(file,
                            header = TRUE,
                            stringsAsFactors = FALSE)
  logger$log('Marker effects table file DONE')


  logger$log('Check marker effects coordinates file ...')
  expectedColumns <- c('SNPid', 'effects')
  if (!all(colnames(markerEffects) %in% expectedColumns)) {
    stop('Marker effects file should have a header specifying those',
         ' columns: `',
         paste(expectedColumns, collapse = '`, `'),
         '`. The detected columns names are: ',
         paste(colnames(markerEffects), collapse = '`, `'), '`.'
    )
  }
  if (nrow(markerEffects) == 0) {
    stop('Marker effects file should have at least one row')
  }

  # missing values
  if (any(is.na(markerEffects))) {
    stop('Marker effects should not have any missing values')
  }

  # duplicated rows
  dupRows <- duplicated(markerEffects)
  if (any(dupRows)) {
    warning(sum(dupRows),
            'duplicated rows found in the marker effects file. ',
            'They will be removed.')
    SNPcoord <- SNPcoord[!dupRows,]
  }

  # check unicity of SNPids
  duplicatedIds <- which(duplicated(markerEffects$SNPid))
  if (length(duplicatedIds) != 0) {
    msg <- paste(
      length(duplicatedIds),
      'duplicated SNPs\' id detected in the marker effects file:',
      paste(SNPcoord$SNPid[duplicatedIds], collapse = ', '))
    stop(msg)
  }
  logger$log('Check marker effects file DONE')

  # reshape the markerEffects data.frame
  row.names(markerEffects) <- markerEffects$SNPid
  markerEffects <- markerEffects[, 'effects', drop = FALSE]
  return(markerEffects)
}



#' Read progeny BLUP estimation file
#'
#' @param file path of the progeny BLUP estimation file generated by
#' r-geno-tools-engine containing the blup estimations of the progenies of some
#' crosses (`json` file).
#'
#' @return data.frame of 4 columns named "ind1", "ind2", "blup_var", "blup_exp"
readProgBlupEstim <- function(file) {
  logger <- Logger$new("r-readProgBlupEstim()")

  logger$log("Read result file ... ")
  if (!file.exists(file)) {
    stop("progeny estimation file do not exists")
  }
  projBlups <- readLines(file)
  logger$log("Read result file DONE ")
  logger$log("Convert Json to data.frame ... ")
  projBlups <- jsonlite::fromJSON(projBlups)
  logger$log("Convert Json to data.frame DONE ")
  if (class(projBlups) != "data.frame"
      || !identical(colnames(projBlups),
                   c("ind1", "ind2", "blup_var", "blup_exp"))) {
    stop("Error: provided file do not seems to be generated by r-geno-tools-Engine.")
  }
  logger$log("DONE, return output.")
  projBlups
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
#' @return path of the created file
saveGWAS <- function(gwasRes, metadata, dir = NULL, file = NULL) {
  logger <- Logger$new("r-saveGWAS()")
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








#' Save relationship matrix in file
#'
#' @param relMat relationship matrix created with `pedRelMat`
#' @param metadata list of metadata of the relationship matrix (optional).
#' @param dir  if `file` is NULL, directory where to save the data,
#' by default it is a temporary directory
#' @param file file path where to save the data. If the file already exists, it
#' will be overwritten. Default NULL (it will create a new "csv" file)
#' @param format format of the output file. Either "csv" or "json". (optional, by
#' default it will use the `file` extension, if file is NULL, "csv").
#'
#' @return path of the created file
saveRelMat <- function(relMat,
                       metadata = NULL,
                       dir = NULL,
                       file = NULL,
                       format = tools::file_ext(file)){

  logger <- Logger$new("r-saveRelMat()")

  logger$log('Check relationship matrix ...')
  checkRelMat(relMat)
  logger$log('Check relationship matrix DONE')

  if (is.null(file)) {
    if (is.null(dir)) {
      dir <- tempdir()
    }
    if (!dir.exists(dir)) {
      logger$log('Check dir ...')
      logger$log('Error: "dir" directory should exists')
      stop('Error: "dir" directory should exists')
      logger$log('Check dir DONE')
    }
    if (length(format) == 0) {
      format <- "csv"
    }
    file <- tempfile(fileext = paste0(".", format), tmpdir = dir)
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

  logger$log('Check file format ...')
  formatErr <- FALSE
  if (length(format) != 1) {
    formatErr <- TRUE
  } else if (!format %in% c('csv', 'json')) {
    formatErr <- TRUE
  }
  if (formatErr) {
    logger$log('Error: File format misspecified. It should be either "csv",
               or "json".')
    stop('Error: File format misspecified. It should be either "csv",
               or "json".')
  }
  logger$log('Check file format DONE')


  if (format == "csv") {
    logger$log('Write relationship matrix in `.csv` file ...')
    if (!is.null(metadata)) {
      metadata <- lapply(metadata, as.character)
      writeLines(paste0("#", names(metadata), "=", unlist(metadata),
                        collapse = "\n"),
                 con = file)
    }
    supThisWarning({
      write.table(x = relMat,
                  file = file,
                  sep = ",",
                  row.names = TRUE,
                  col.names = TRUE,
                  append = !is.null(metadata))
    }, "appending column names to file")
    logger$log('Write relationship matrix in `.csv` file DONE')
  }
  if (format == "json") {
    logger$log('Write relationship matrix in `.json` file ...')
    relMatList <- list(relMat = as.data.frame(relMat),
                       metadata = metadata)
    writeLines(jsonlite::toJSON(relMatList,
                                dataframe = "rows",
                                pretty = T,
                                digits = NA,
                                na = 'string'),
               con = file)
    logger$log('Write relationship matrix in `.json` file DONE')
  }

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
  logger <- Logger$new("r-prepareData()")
  # Remove from geno data individuals that are not in phenotypic data-set
  logger$log("Remove from geno data individuals that are not in phenotypic data-set ...")
  gDta <- gaston::select.inds(gDta, id %in% rownames(pDta))
  logger$log("Remove from geno data individuals that are not in phenotypic data-set DONE")


  # reorder phenotypic data with id in bed matrix
  logger$log("reorder matrix ...")
  pDta <- pDta[gDta@ped$id,,drop=F]
  logger$log("reorder matrix DONE")


  # remove monomorphic markers
  logger$log("remove monomorphic markers ...")
  gDta <- gaston::select.snps(gDta, maf > 0)
  logger$log("remove monomorphic markers DONE")

  logger$log("DONE, return output.")

  list(genoData = gDta, phenoData = pDta)
}


#' Save phased genotypes of simulatied population to vcf.gz file
#'
#' @param file file path where to save the data. If the file already exists, it
#' will be overwritten.
#' @param pop simulated population (`breedSimulatR`'s population)
#' @param SNPcoord snp coordinate of the genotypes. (data.frame with `chr`, `physPos`, and `SNPid` columns)
#'
saveVcf <- function(file, pop, SNPcoord){

  # Fixed region
  data <- SNPcoord[, c("chr", "physPos", "SNPid")]
  data$physPos[is.na(data$physPos)] <- '.'

  colnames(data) <- c("#CHROM", "POS", "ID")
  data <- data[order(data$POS),]
  data <- data[order(data$`#CHROM`),]

  data$REF <- "."
  data$ALT <- "."
  data$QUAL <- "."
  data$FILTER <- "PASS"
  data$INFO <- "."

  # Genotype region
  data$FORMAT <- "GT"
  gt <- vapply(pop$inds, function(ind){
    hap <- do.call(cbind, ind$haplo$values)
    x <- paste(hap[1,], hap[2,], sep = "|")
    names(x) <- colnames(hap)
    x
  }, vector(mode = "character",
            length = length(pop$inds[[1]]$haplo$allelDose)))
  gt <- as.data.frame(gt[data$ID,])
  data <- cbind(data, gt)

  # Meta region
  meta <- paste("##fileformat=VCFv4.3",
                "##source=\"R-Geno-tool-engine\", data in this file are simulated.",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                sep = "\n")

  # write file
  f <- gzfile(file, "w")
  writeLines(text = meta, con = f)
  close(f)
  data.table::fwrite(x = data,
                     file = file,
                     append = TRUE,
                     sep = "\t",
                     quote = FALSE,
                     row.names = FALSE,
                     col.names = TRUE)
}





#' Save a R data.frame as json file
#'
#' @param df data.frame
#' @param file file path where to save the data. If the file already exists, it
#' will be overwritten.
#'
#' @return path of the created file
save_dataFrame_as_json <- function(df, file){
  logger <- Logger$new("r-save_dataFrame_as_json()")

    logger$log('Check file ...')
    if (length(file) != 1) {
      logger$log('Error: only one file name should be provided')
      stop('Error: only one file name should be provided')
    }

    logger$log("Check output file extention ...")
    ext <- tools::file_ext(file)
    if (ext != "json") {
      stop('The output file must end by `.json`')
    }
    logger$log("Check output file extention DONE")

    if (file.exists(file)) {
      logger$log('Warning: "file" directory already exists. This file will be overwritten.')
    } else {
      file.create(file)
    }
    logger$log('Check file DONE')

    writeLines(jsonlite::toJSON(df,
                                dataframe = "rows",
                                pretty = T,
                                digits = NA,
                                na = 'string'),
               con = file)
    return(file)
}
