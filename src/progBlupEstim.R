# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# functions related to progenies's genetic values variance and expected values.


#' Calculate the recombination rate matrix for each couple of SNP
#'
#' The recombination rate is alculated using the "haldane inverse" function
#'
#' @param SNPcoord SNP coordinate data.frame return by `readSNPcoord`
#'
#' @return named list of matrices. List names are the chromosomes' names.
#' Matrices' row and columns names are the SNP ids
#' (of the corresponding chromosome). Matrices values are the recombination rate
#' between the corresponding SNPs.
calcRecombRate <- function(SNPcoord) {
  # calculate recombination rate for each chromosome
  r <- lapply(unique(SNPcoord$chr), function(chr){
    subsetSNPcoord <- SNPcoord[SNPcoord$chr == chr,]
    r <- as.matrix(dist(subsetSNPcoord$linkMapPos))
    r <- 0.5 * (1 - exp(-2 * r)) # haldane inverse function
    colnames(r) <- row.names(r) <- subsetSNPcoord$SNPid
    r
  })
  names(r) <- unique(SNPcoord$chr)
  r
}



#' Calculate the genetic variance-covariance of matrix the progenies of 2 given
#' parents
#'
#' @param SNPcoord SNP coordinate data.frame return by `readSNPcoord`
#' @param r recombination rate matrices return by `calcRecombRate`
#' @param haplo haplotypes of individuals ("haplotypes" element of the list
#' return by `readPhasedGeno` function)
#' @param p1.id id of the first parent
#' @param p2.id id of the second parent
#'
#' @return named list of matrices. List names are the chromosomes' names.
#' Matrices' row and columns names are the SNP ids
#' (of the corresponding chromosome). Matrices values are the genetic covariance
#' between the corresponding SNPs for the progeny of the given parents.
calcProgenyGenetCovar <- function(SNPcoord, r, haplo, p1.id, p2.id) {
  geneticCovar <- lapply(unique(SNPcoord$chr), function(chr){
    subsetSNPcoord <- SNPcoord[SNPcoord$chr == chr,]
    SNPs <- row.names(r[[chr]])
    nSNPs <- length(SNPs)
    haplo_parent1 <- haplo[SNPs, p1.id]
    haplo_parent2 <- haplo[SNPs, p2.id]

    matHaplo_1.1 <- matrix(haplo_parent1[,1], nSNPs, nSNPs, byrow = T)
    matHaplo_2.1 <- matrix(haplo_parent1[,2], nSNPs, nSNPs, byrow = T)
    t_matHaplo_1.1 <- matrix(haplo_parent1[,1], nSNPs, nSNPs, byrow = F)
    t_matHaplo_2.1 <- matrix(haplo_parent1[,2], nSNPs, nSNPs, byrow = F)
    z1 <- matHaplo_1.1 * t_matHaplo_1.1 + matHaplo_2.1 * t_matHaplo_2.1 -
      (matHaplo_1.1 * t_matHaplo_2.1 + t_matHaplo_1.1 * matHaplo_2.1)

    matHaplo_1.2 <- matrix(haplo_parent2[,1], nSNPs, nSNPs, byrow = T)
    matHaplo_2.2 <- matrix(haplo_parent2[,2], nSNPs, nSNPs, byrow = T)
    t_matHaplo_1.2 <- matrix(haplo_parent2[,1], nSNPs, nSNPs, byrow = F)
    t_matHaplo_2.2 <- matrix(haplo_parent2[,2], nSNPs, nSNPs, byrow = F)
    z <- z1 + matHaplo_1.2 * t_matHaplo_1.2 + matHaplo_2.2 * t_matHaplo_2.2 -
      (matHaplo_1.2 * t_matHaplo_2.2 + t_matHaplo_1.2 * matHaplo_2.2)

    covar <- 0.25 * (1 - 2 * r[[chr]]) * (z)
    row.names(covar) <- colnames(covar) <- SNPs
    covar
  })

  names(geneticCovar) <- unique(SNPcoord$chr)
  geneticCovar
}



#' Calculate the BULP variance of the progeny
#'
#' @param SNPcoord SNP coordinate data.frame return by `readSNPcoord`
#' @param markerEffects data.frame of the markers effects return by `readMarkerEffects`
#' @param geneticCovar list of the genetic variance covariance matrices return by
#' `calcProgenyGenetCovar`
#'
#' @return numeric
calcProgenyBlupVariance <- function(SNPcoord, markerEffects, geneticCovar) {
  blupVar <- lapply(unique(SNPcoord$chr), function(chr){
    eff <- markerEffects[row.names(geneticCovar[[chr]]), ]
    as.numeric(t(eff) %*% geneticCovar[[chr]] %*% eff)
  })
  blupVar <- sum(unlist(blupVar))
  blupVar
}

#' Calculate the BULP expected values of the progeny
#'
#' @param SNPcoord SNP coordinate data.frame return by `readSNPcoord`
#' @param haplo haplotypes of individuals ("haplotypes" element of the list
#' return by `readPhasedGeno` function)
#' @param p1.id id of the first parent
#' @param p2.id id of the second parent
#' @param markerEffects data.frame of the markers effects return by `readMarkerEffects`
#'
#' @return numeric
calcProgenyBlupExpected <- function(SNPcoord, haplo, p1.id, p2.id, markerEffects) {
  blupExp <- lapply(unique(SNPcoord$chr), function(chr){
    subsetSNPcoord <- SNPcoord[SNPcoord$chr == chr,]
    SNPs <- subsetSNPcoord$SNPid
    haplo_p1 <- haplo[SNPs, p1.id]
    haplo_p2 <- haplo[SNPs, p2.id]

    eff <- markerEffects[SNPs, ]

    0.5 * eff %*% (haplo_p1[,1] + haplo_p1[,2] + haplo_p2[,1] + haplo_p2[,2])
  })
  blupExp <- sum(unlist(blupExp))
  blupExp
}
