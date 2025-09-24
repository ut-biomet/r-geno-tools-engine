# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# unit test file for crossing simulation functions



capture_output({
  # initializeSimulation ----
  phasedGenoFile <- "../../data/geno/breedGame_phasedGeno.vcf.gz"
  SNPcoordFile <- "../../data/SNPcoordinates/breedingGame_SNPcoord.csv"

  test_that("initializeSimulation", {
    g <- readPhasedGeno(phasedGenoFile)
    SNPcoord <- checkAndFilterSNPcoord(readSNPcoord(SNPcoordFile), g$SNPcoord)

    expect_no_error({
      parent_pop <- initializeSimulation(g$haplotypes, SNPcoord)
    })
    expect_is(parent_pop, "population")
    expect_equal(parent_pop$nInd, ncol(g$haplotypes) / 2)

    # chromosome information
    chrInfo_max <- aggregate(SNPcoord, by = list(name = SNPcoord$chr), FUN = max)
    chrInfo_max <- chrInfo_max[, c("name", "physPos", "linkMapPos")]
    colnames(chrInfo_max) <- c("name", "max_physPos", "max_linkMapPos")
    chrInfo_min <- aggregate(SNPcoord, by = list(name = SNPcoord$chr), FUN = min)
    chrInfo_min <- chrInfo_min[, c("name", "physPos", "linkMapPos")]
    colnames(chrInfo_min) <- c("name", "min_physPos", "min_linkMapPos")
    chrInfo <- dplyr::left_join(chrInfo_min, chrInfo_max)
    chrInfo$chrLength <- chrInfo$max_linkMapPos - chrInfo$min_linkMapPos

    expect_equal(sort(parent_pop$specie$chrNames), sort(chrInfo$name))

    chrLength <- chrInfo$chrLength * 10^2
    names(chrLength) <- chrInfo$name
    chrLength <- chrLength[parent_pop$specie$chrNames]
    expect_equal(parent_pop$specie$lchrCm, chrLength)

    # SNP coordinates:
    simulInit_SNPcoord <- parent_pop$inds[[1]]$haplo$SNPinfo$SNPcoord
    simulInit_SNPcoord <- simulInit_SNPcoord[order(simulInit_SNPcoord$SNPid), ]

    # transform SNP linkage map position to start at 0 (simple translation)
    for (chr in unique(SNPcoord$chr)) {
      selLines <- SNPcoord$chr == chr # (selected)lines of the current chromosome
      minLMPos <- min(SNPcoord[selLines, "linkMapPos"])
      SNPcoord[selLines, "linkMapPos"] <- SNPcoord$linkMapPos[selLines] - minLMPos
    }

    expect_equal(sort(simulInit_SNPcoord$SNPid), sort(SNPcoord$SNPid))
    expect_equal(
      simulInit_SNPcoord$linkMapPos,
      SNPcoord[order(SNPcoord$SNPid), "linkMapPos"] * 10^2
    )
    expect_equal(
      simulInit_SNPcoord$chr,
      SNPcoord[order(SNPcoord$SNPid), "chr"]
    )
  })
})
