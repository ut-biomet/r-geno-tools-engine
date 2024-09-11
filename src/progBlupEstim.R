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



#' Calculate the genetic variance-covariance matrix of the progenies of 2 given
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
#' @param makrerEffects (output of `extract_additive_effects` function) list of
#' 2 elements:
#' `intercept`: named vector of the intercepts
#' `SNPeffects`: data.frame of the additive effects, 1 columns per phenotype with
#' the marker ids as row names.
#' @param geneticCovar list of the genetic variance covariance matrices return by
#' `calcProgenyGenetCovar`
#'
#' @return numeric
calcProgenyBlupVariance <- function(SNPcoord, markerEffects, geneticCovar) {
  apply(markerEffects$SNPeffects, 2, function(markEff){
    blupVar <- lapply(unique(SNPcoord$chr), function(chr){
      eff <- markEff[row.names(geneticCovar[[chr]])]
      as.numeric(t(eff) %*% geneticCovar[[chr]] %*% eff)
    })
    blupVar <- sum(unlist(blupVar))
    blupVar
  })
}

#' Calculate the BULP expected values of the progeny
#'
#' @param SNPcoord SNP coordinate data.frame return by `readSNPcoord`
#' @param haplo haplotypes of individuals ("haplotypes" element of the list
#' return by `readPhasedGeno` function)
#' @param p1.id id of the first parent
#' @param p2.id id of the second parent
#' @param makrerEffects (output of `extract_additive_effects` function) list of
#' 2 elements:
#' `intercept`: named vector of the intercepts
#' `SNPeffects`: data.frame of the additive effects, 1 columns per phenotype with
#' the marker ids as row names.
#'
#' @return list for each trait with a list with
#' - `sum` the global expected value for the trait (taking in account the intercept)
#' - `by_chr` list of the expected value for each chromosome (NOT taking in account the intercept)
calcProgenyBlupExpected <- function(SNPcoord, haplo, p1.id, p2.id, markerEffects) {
  blupExp <- lapply(names(markerEffects$intercept), function(trait){
    markEff <- markerEffects$SNPeffects[, trait, drop = FALSE]
    blupExp <- lapply(unique(SNPcoord$chr), function(chr){
      subsetSNPcoord <- SNPcoord[SNPcoord$chr == chr,]
      SNPs <- subsetSNPcoord$SNPid
      haplo_p1 <- haplo[SNPs, p1.id]
      haplo_p2 <- haplo[SNPs, p2.id]
      eff <- markEff[SNPs,]
      0.5 * eff %*% (haplo_p1[,1] + haplo_p1[,2] + haplo_p2[,1] + haplo_p2[,2])
    })
    names(blupExp) <- unique(SNPcoord$chr)
    list(
      sum = sum(unlist(blupExp)) + markerEffects$intercept[trait],
      by_chr = blupExp
    )
  })
  names(blupExp) <- names(markerEffects$intercept)
  blupExp
}


#' Calculate the variance covariance matrix of the progeny's blups for several
#' traits
#'
#' @param SNPcoord SNP coordinate data.frame return by `readSNPcoord`
#' @param r recombination rate matrices return by `calcRecombRate`
#' @param haplo haplotypes of individuals ("haplotypes" element of the list
#' return by `readPhasedGeno` function)
#' @param p1.id id of the first parent
#' @param p2.id id of the second parent
#' @param makrerEffects (output of `extract_additive_effects` function) list of
#' 2 elements:
#' `intercept`: named vector of the intercepts
#' `SNPeffects`: data.frame of the additive effects, 1 columns per phenotype with
#' the marker ids as row names.
#' @param blupExpectedValues output of `calcProgenyBlupExpected`
#'
#' @return matrix
calcProgenyBlupCovariance <- function(SNPcoord,
                                      r,
                                      haplo,
                                      p1.id,
                                      p2.id,
                                      markerEffects,
                                      blupExpectedValues) {

  EG <- blupExpectedValues
  CovG1G2 <- matrix(nrow = length(EG), ncol = length(EG))
  colnames(CovG1G2) <- names(markerEffects$intercept)
  rownames(CovG1G2) <- names(markerEffects$intercept)
  for (t1 in names(markerEffects$intercept)) {
    for (t2 in names(markerEffects$intercept)) {

      cov <- lapply(unique(SNPcoord$chr), function(chr){

        subsetSNPcoord <- SNPcoord[SNPcoord$chr == chr,]
        SNPs <- subsetSNPcoord$SNPid

        nSNPs <- length(SNPs)
        eff_1 <- markerEffects$SNPeffects[SNPs, t1]
        eff_2 <- markerEffects$SNPeffects[SNPs, t2]

        haplo_parent1 <- haplo[SNPs, p1.id]
        haplo_parent2 <- haplo[SNPs, p2.id]

        x1.m <- matrix(haplo_parent1[,1], nSNPs, nSNPs, byrow = T)
        x1.p <- matrix(haplo_parent1[,2], nSNPs, nSNPs, byrow = T)
        y1.m <- matrix(haplo_parent1[,1], nSNPs, nSNPs, byrow = F)
        y1.p <- matrix(haplo_parent1[,2], nSNPs, nSNPs, byrow = F)

        x2.m <- matrix(haplo_parent2[,1], nSNPs, nSNPs, byrow = T)
        x2.p <- matrix(haplo_parent2[,2], nSNPs, nSNPs, byrow = T)
        y2.m <- matrix(haplo_parent2[,1], nSNPs, nSNPs, byrow = F)
        y2.p <- matrix(haplo_parent2[,2], nSNPs, nSNPs, byrow = F)

        EX1Y1 <- (0.25 * (2*x1.m*y1.m + 2*x1.p*y1.p) -
                    (r[[chr]] / 2) * (x1.m*y1.m - x1.m*y1.p - x1.p*y1.m + x1.p*y1.p))
        EX2Y2 <- (0.25 * (2*x2.m*y2.m + 2*x2.p*y2.p) -
                    (r[[chr]] / 2) * (x2.m*y2.m - x2.m*y2.p - x2.p*y2.m + x2.p*y2.p))
        EX1Y2 <- 0.25 * (x1.m*y2.m + x1.m*y2.p + x1.p*y2.m + x1.p*y2.p)
        EX2Y1 <- 0.25 * (x2.m*y1.m + x2.m*y1.p + x2.p*y1.m + x2.p*y1.p)

        EX1X2Y1Y2 <- EX1Y1 + EX2Y2 + EX1Y2 + EX2Y1

        EG1G2 <- as.numeric(t(eff_1) %*% EX1X2Y1Y2 %*% eff_2)

        cov = EG1G2 - EG[[t1]]$by_chr[[chr]] * EG[[t2]]$by_chr[[chr]]
        cov
      })
      cov <-  sum(unlist(cov))
      if (isTRUE(all.equal(cov, 0))) {
        # sometime cov is not exactly equal to 0 due to computer precision
        # and the result could be negative. This can be a problem for the
        # diagonal of the covariance matrix
        cov <- 0
      }
      CovG1G2[t1, t2] <- sum(unlist(cov))
    }
  }
  CovG1G2
}
