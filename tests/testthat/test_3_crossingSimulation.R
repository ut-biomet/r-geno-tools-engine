# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# unit test file for crossing simulation functions



capture_output({

  # initializeSimulation ----
  phasedGenoFile <- '../../data/geno/breedGame_phasedGeno.vcf.gz'
  chrInfoFile <- '../../data/chromosomesInformation/breedingGame_chrInfo.csv'
  SNPcoordFile <- '../../data/SNPcoordinates/breedingGame_SNPcoord.csv'

  test_that('initializeSimulation', {
    g <- readPhasedGeno(phasedGenoFile)
    chrInfo <- readChrInfo(chrInfoFile)
    SNPcoord <- checkAndFilterSNPcoord(readSNPcoord(SNPcoordFile), g$SNPcoord)

    expect_error({
      parent_pop <- initializeSimulation(chrInfo, g$haplotypes, SNPcoord)
    }, NA)
    expect_is(parent_pop, "population")
    expect_equal(parent_pop$nInd, ncol(g$haplotypes)/2)

    # chromosome information
    expect_equal(sort(parent_pop$specie$chrNames), sort(chrInfo$name))

    chrLength_cm <- chrInfo$length_morgan * 10^2
    names(chrLength_cm) <- chrInfo$name
    chrLength_cm <- chrLength_cm[parent_pop$specie$chrNames]
    expect_equal(parent_pop$specie$lchrCm, chrLength_cm)

    chrLength_physPos <- chrInfo$length_phys
    names(chrLength_physPos) <- chrInfo$name
    chrLength_physPos <- chrLength_physPos[parent_pop$specie$chrNames]
    expect_equal(parent_pop$specie$lchr, chrLength_physPos)


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
