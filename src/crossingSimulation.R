# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Functions related to crossing simulation

#' Initialise the simulation
#'
#' @param haplotypes haplotypes of the parents, data.frame with genotype values
#' in row, and individuals'haplotype in columns. The columns name should be
#' `individualName_1` and `individualName_2` for the first/second haplotype of
#' the individual named `individualName`. (the list item named `haplo` return
#' by the function `readPhasedGeno`)
#' @param SNPcoord snp coordinates, data.frame of 4 columns:
#'  - `chr`: Name of the chromosome holding the SNP
#'  - `physPos`: SNP physical position on the chromosome
#'  - `linkMapPos`: SNP linkage map position on the chromosome **in morgan**
#'  - `SNPid`: SNP's IDs
#'
#' @return `breedSimulatR`'s population object
initializeSimulation <- function(haplotypes,
                                 SNPcoord) {
  logger <- Logger$new("r-initializeSimulation()")

  # initialisation
  # transform SNP linkage map position to start at 0 (simple translation)
  for (chr in unique(SNPcoord$chr)) {
    selLines <- SNPcoord$chr == chr # (selected)lines of the current chromosome
    minLMPos <- min(SNPcoord[selLines, "linkMapPos"])
    SNPcoord[selLines, "linkMapPos"] <- SNPcoord$linkMapPos[selLines] - minLMPos
  }

  # get chr size with max SNP coordinates for each chromosome
  SNPcoord$physPos <- seq(nrow(SNPcoord)) # dummy data
  logger$log("Extract chromosomes information ...")
  chrInfo <- aggregate(SNPcoord, by = list(name = SNPcoord$chr), FUN = max)
  chrInfo <- chrInfo[, c("name", "physPos", "linkMapPos")]
  colnames(chrInfo) <- c("name", "max_physPos", "max_linkMapPos")

  # create specie object
  logger$log("Create specie ...")
  specie <- breedSimulatR::specie$new(
    nChr = nrow(chrInfo),
    chrNames = chrInfo$name,
    lchr = nrow(SNPcoord), # dummy data
    lchrCm = chrInfo$max_linkMapPos * 10^2,
    verbose = FALSE
  )
  logger$log("Create specie DONE")

  # create snp information object
  logger$log("Create snp information ...")
  SNPcoord$linkMapPos <- SNPcoord$linkMapPos * 10^2 # must be in centi-morgan
  snpInfo <- breedSimulatR::SNPinfo$new(
    SNPcoord = SNPcoord,
    specie = specie
  )
  logger$log("Create snp information DONE")

  # create population object
  logger$log("Create parents population ...")
  listInds <- vector(mode = "list", length = ncol(haplotypes) / 2)
  indNames <- gsub("_[12]$", "", colnames(haplotypes)[1:(ncol(haplotypes) / 2)])
  names(listInds) <- indNames

  for (indName in indNames) {
    haplo <- t(haplotypes[, paste(indName, c(1, 2), sep = "_")])
    haplo <- breedSimulatR::haplotype$new(
      SNPinfo = snpInfo,
      haplo = haplo
    )
    listInds[[indName]] <- breedSimulatR::individual$new(
      name = indName,
      specie = specie,
      haplo = haplo,
      verbose = FALSE
    )
  }
  pop <- breedSimulatR::population$new(
    name = "", inds = listInds,
    verbose = FALSE
  )
  logger$log("Create parents population DONE")
  return(pop)
}
