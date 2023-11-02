# Author: Julien Diot juliendiot@ut-biomet.org
# 2023 The University of Tokyo
#
# Description:
# unit test of progeny blup estimation functions




capture_output({
  file_list <- list(
    oneTrait = list(
      genoFile = '../../data/geno/breedGame_phasedGeno.vcf.gz',
      SNPcoordFile = '../../data/SNPcoordinates/breedingGame_SNPcoord.csv',
      crossTableFile = '../../data/crossingTable/breedGame_small_crossTable.csv',
      markerEffectsFile = '../../data/markerEffects/breedGame_markerEffects.json'
    ),
    twoTraits = list(
      genoFile = '../../data/geno/breedGame_phasedGeno.vcf.gz',
      SNPcoordFile = '../../data/SNPcoordinates/breedingGame_SNPcoord.csv',
      crossTableFile = '../../data/crossingTable/breedGame_small_crossTable.csv',
      markerEffectsFile = '../../data/markerEffects/breedGame_markerEffects_2traits.json'
    )
  )

  for (testName in names(file_list)) {
    files <- file_list[[testName]]
    g <- readPhasedGeno(files$genoFile)
    SNPcoord <- readSNPcoord(files$SNPcoordFile)
    crossTable <- readCrossTable(files$crossTableFile)
    markerEffects <- readMarkerEffects(files$markerEffectsFile)

    r <- calcRecombRate(SNPcoord)
    nCrosses <- nrow(crossTable)

    test_that(paste("calcProgenyBlupCovariance", testName), {
      for (i in seq(nCrosses)) {
        p1.id <- which(grepl(crossTable$ind1[i], colnames(g$haplotypes)))
        p2.id <- which(grepl(crossTable$ind2[i], colnames(g$haplotypes)))
        expect_error({
          blupExp <- calcProgenyBlupExpected(SNPcoord,
                                             g$haplotypes,
                                             p1.id,
                                             p2.id,
                                             markerEffects)
        }, NA)
        expect_error({
          blupCovar <- calcProgenyBlupCovariance(SNPcoord = SNPcoord,
                                                 r = r,
                                                 haplo = g$haplotypes,
                                                 p1.id = p1.id,
                                                 p2.id = p2.id,
                                                 markerEffects = markerEffects,
                                                 blupExpectedValues = blupExp)
        }, NA)
        expect_error({
          geneticCovar <- calcProgenyGenetCovar(SNPcoord, r, g$haplotypes, p1.id, p2.id)
          blupVar <- calcProgenyBlupVariance(SNPcoord, markerEffects, geneticCovar)
        }, NA)
        expect_equal(diag(blupCovar), blupVar)
      }
    })
  }
})
