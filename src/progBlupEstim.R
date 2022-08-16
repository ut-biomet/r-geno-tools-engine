# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# functions related to progenies's genetic values variance and expected values.


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



calcProgenyBlupVariance <- function(SNPcoord, markerEffects, geneticCovar) {
  blupVar <- lapply(unique(SNPcoord$chr), function(chr){
    eff <- markerEffects[row.names(geneticCovar[[chr]]), ]
    as.numeric(t(eff) %*% geneticCovar[[chr]] %*% eff)
  })
  blupVar <- sum(unlist(blupVar))
  blupVar
}

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




