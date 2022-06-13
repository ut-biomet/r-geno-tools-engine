# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# unit test file for crossing simulation functions



capture_output({

  # initializeSimulation ----
  phasedGenoFile <- '../../data/geno/breedGame_phasedGeno.vcf.gz'
  SNPcoordFile <- '../../data/SNPcoordinates/breedingGame_SNPcoord.csv'

  test_that('initializeSimulation', {
    g <- readPhasedGeno(phasedGenoFile)
    SNPcoord <- checkAndFilterSNPcoord(readSNPcoord(SNPcoordFile), g$SNPcoord)

    expect_error({
      parent_pop <- initializeSimulation(g$haplotypes, SNPcoord)
    }, NA)
    expect_is(parent_pop, "population")
    expect_equal(parent_pop$nInd, ncol(g$haplotypes)/2)

    # chromosome information
    chrInfo <- aggregate(SNPcoord, by = list(name = SNPcoord$chr), FUN = max)
    chrInfo <- chrInfo[, c('name', 'physPos', 'linkMapPos')]
    colnames(chrInfo) <- c('name', 'max_physPos', 'max_linkMapPos')

    expect_equal(sort(parent_pop$specie$chrNames), sort(chrInfo$name))

    max_linkMapPos <- chrInfo$max_linkMapPos * 10^2
    names(max_linkMapPos) <- chrInfo$name
    max_linkMapPos <- max_linkMapPos[parent_pop$specie$chrNames]
    expect_equal(parent_pop$specie$lchrCm, max_linkMapPos)

    max_physPos <- chrInfo$max_physPos
    names(max_physPos) <- chrInfo$name
    max_physPos <- max_physPos[parent_pop$specie$chrNames]
    expect_equal(parent_pop$specie$lchr, max_physPos)


    # SNP coordinates:
    simulInit_SNPcoord <- parent_pop$inds[[1]]$haplo$SNPinfo$SNPcoord
    simulInit_SNPcoord <- simulInit_SNPcoord[order(simulInit_SNPcoord$SNPid),]
    expect_equal(sort(simulInit_SNPcoord$SNPid), sort(SNPcoord$SNPid))
    expect_equal(simulInit_SNPcoord$linkMapPos,
                 SNPcoord[order(SNPcoord$SNPid), 'linkMapPos'] * 10^2)
    expect_equal(simulInit_SNPcoord$physPos,
                 SNPcoord[order(SNPcoord$SNPid), 'physPos'])
    expect_equal(simulInit_SNPcoord$chr,
                 SNPcoord[order(SNPcoord$SNPid), 'chr'])
    })



})
