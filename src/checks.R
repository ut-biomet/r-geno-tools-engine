# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# function to check r-geno-tool-engine objects





checkRelMat <- function(relMat) {
  errMessage <- paste('Bad relationship information,',
                      'please be sure the relationship matrix file',
                      'have been created using `r-geno-tool-engine`.')

  if (!is.numeric(relMat)) {
    stop(errMessage, '`relMat` is not numeric.')
  }
  cols <- sort(colnames(relMat))
  rows <- sort(row.names(relMat))
  if (!identical(cols, rows)) {
    stop(errMessage, "`relMat`'s rows and columns names are different.")
  }
  # if (!isSymmetric(relMat)) {
  #   stop(errMessage, '`relMat` is not symmetric.')
  # }
  if (!all.equal.numeric(relMat, t(relMat))) {
    stop(errMessage, '`relMat` is not approximately symmetric.')
  }

  return(TRUE)
}


#' Check compatibility between 2 snps coordinates data set
#' and keep only genotypes SNPs
#'
#' @param user_SNPcoord SNPs coordinates data coming from the user (`.csv` file)
#' @param vcf_SNPcoord SNPs coordinates data coming from the `.vcf` file
#'
#' @description If the `user_SNPcoord` data miss some SNPs defined in the `.vcf`
#' file, an error is raised. If the `user_SNPcoord` data have additional SNPs,
#' those SNPs will be removed.
#' If the data between those two data-set are inconsistent, an error is raised.
#'
#' @return The filtered `user_SNPcoord` data.frame
checkAndFilterSNPcoord <- function(user_SNPcoord, vcf_SNPcoord) {
  # SNP ids
  user_SNPid <- user_SNPcoord$SNPid
  vcf_SNPid <- vcf_SNPcoord$SNPid
  missingSNP <- which(!vcf_SNPid %in% user_SNPid)
  if (length(missingSNP) != 0) {
    msg <- paste('The SNPs coordinate file miss', length(missingSNP),
                 'genotype\'s SNPs:',
                 paste(vcf_SNPid[missingSNP], collapse = ', '))
    stop(msg)
  }

  additionalSNP <- which(!user_SNPid %in% vcf_SNPid)
  if (length(additionalSNP) != 0) {
    msg <- paste('The SNPs coordinate file have', length(additionalSNP),
                 'SNPs not defined in the `.vcf` file,',
                 'those SNP will not be considered:',
                 paste(user_SNPid[additionalSNP], collapse = ', '))
    warning(msg)
    user_SNPcoord <- user_SNPcoord[-additionalSNP,]
  }

  # Check they have the same data
  user_SNPcoord_tmp <- user_SNPcoord[order(user_SNPcoord$SNPid),]
  vcf_SNPcoord <- vcf_SNPcoord[order(vcf_SNPcoord$SNPid),]

  linkMapColumnId <- which(colnames(user_SNPcoord) == 'linkMapPos')
  user_SNPcoord_tmp <- user_SNPcoord_tmp[, -linkMapColumnId]

  columns <- c('chr', 'SNPid', 'physPos')
  vcf_SNPcoord <- vcf_SNPcoord[, columns]
  user_SNPcoord_tmp <- user_SNPcoord_tmp[, columns]
  diffId <- which(user_SNPcoord_tmp != vcf_SNPcoord)

  if (length(diffId) != 0) {
    diffLineId <- diffId %% nrow(user_SNPcoord_tmp)
    errMsg <- paste('The SNPs coordinate file and the genotype files have',
                    'inconsistent data.',
                    length(diffLineId),
                    'incompatibilities detected about SNPs:',
                    paste(vcf_SNPcoord$SNPid[diffLineId], collapse = ', '))
    stop(errMsg)
  }

  return(user_SNPcoord)
}



#' Check individuals in the crossing table are in the haplotype data
#'
#' @param crossTable the crossing table
#' @param haplo the haplotype data given by the function `readPhasedGeno`
#'
#' @return NULL, raise error if missing individuals are detected.
checkIndNamesConsistency <- function(crossTable, haplo) {
  crossTableInds <- unique(c(crossTable$ind1, crossTable$ind2))
  haploInds <- gsub('_[12]$', '', colnames(haplo)[1:(ncol(haplo)/2)])

  missIndsId <- which(!crossTableInds %in% haploInds)
  if (length(missIndsId) != 0) {
    msg <- paste(
      length(missIndsId),
      'individuals are defined in the crossing table',
      'but not in the genoype file:',
      paste(crossTableInds[missIndsId], collapse = ', ')
    )
    stop(msg)
  }
  invisible(NULL)
}



#' Check the SNPs' coordinates are consistent with the chromosomes information.
#'
#' @param chrInfo output of `readChrInfo`
#' @param SNPcoord output of `readSNPcoord`
#'
#' @return NULL, raise error if not consistent
checkChrInfoConsistency <- function(chrInfo, SNPcoord) {

  # chromosomes names
  if (!identical(sort(chrInfo$name), sort(unique(SNPcoord$chr)))) {
    stop("Chromosomes names in chromosomes information and SNPs coordinates",
         " files do not match.")
  }

  # chromosomes length
  colnames(chrInfo)[1] <- 'chr'
  suppressMessages({SNPcoord <- dplyr::full_join(SNPcoord, chrInfo,)})
  exceedingSNPs <- which(SNPcoord$physPos > SNPcoord$length_phys)
  if (length(exceedingSNPs) != 0) {
    msg <- paste(
      length(exceedingSNPs),
      'SNPs exceed their chromosome\'s physical length:',
      paste(SNPcoord$SNPid[exceedingSNPs], collapse = ', '))
    stop(msg)
  }

  exceedingSNPs <- which(SNPcoord$linkMapPos > SNPcoord$length_morgan)
  if (length(exceedingSNPs) != 0) {
    msg <- paste(
      length(exceedingSNPs),
      'SNPs exceed their chromosome\'s linkage map length:',
      paste(SNPcoord$SNPid[exceedingSNPs], collapse = ', '))
    stop(msg)
  }

  invisible(NULL)
}
