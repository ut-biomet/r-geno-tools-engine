# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# unit test file for checks functions



capture_output({

  # checkAndFilterSNPcoord ----
  phasedGenoFile <- '../../data/geno/breedGame_phasedGeno.vcf.gz'
  SNPcoordFile <- '../../data/SNPcoordinates/breedingGame_SNPcoord.csv'
  test_that(paste('checkAndFilterSNPcoord'), {
    g <- readPhasedGeno(phasedGenoFile)
    SNPcoord <- readSNPcoord(SNPcoordFile)

    # no error / no filter
    expect_error({
      filteredSNPcoord <- checkAndFilterSNPcoord(SNPcoord, g$SNPcoord)
    }, NA)
    expect_identical(filteredSNPcoord, SNPcoord)

    # additional SNP in user's SNP (should be removed)
    SNPcoord_additionalSNP <- SNPcoord
    SNPcoord_additionalSNP[nrow(SNPcoord) + 1, ] <- list('chr1',
                                                         42,
                                                         'additionalSNP',
                                                         0.00042)
    expect_warning({
      filteredSNPcoord_addSNP <- checkAndFilterSNPcoord(SNPcoord_additionalSNP,
                                                        g$SNPcoord)
    }, paste('The SNPs coordinate file have 1 SNPs not defined in the `.vcf`',
             'file, those SNP will not be considered: additionalSNP'))
    expect_equal(filteredSNPcoord_addSNP, SNPcoord)


    # missing SNP in user's SNP (error)
    SNPcoord_missingSNP <- SNPcoord[-nrow(SNPcoord),]
    expect_error({
      checkAndFilterSNPcoord(SNPcoord_missingSNP, g$SNPcoord)
    }, paste("The SNPs coordinate file miss 1 genotype's SNPs:",
             SNPcoord$SNPid[nrow(SNPcoord)]))


    # inconsistent physical position between user's snp and vcf snp (error)
    SNPcoord_inconsistPhysPos <- SNPcoord
    SNPcoord_inconsistPhysPos[1, 'physPos'] <- SNPcoord[1, 'physPos'] + 1
    expect_error({
      checkAndFilterSNPcoord(SNPcoord_inconsistPhysPos, g$SNPcoord)
    }, paste('The SNPs coordinate file and the genotype files have',
             'inconsistent data.'))

  })



  # checkIndNamesConsistency ----
  crossTableFile <- '../../data/crossingTable/breedGame_crossTable.csv'
  phasedGenoFile <- '../../data/geno/breedGame_phasedGeno.vcf.gz'
  test_that('checkIndNamesConsistency', {
    g <- readPhasedGeno(phasedGenoFile)
    crossTable <- readCrossTable(crossTableFile)

    # no error
    expect_error({
      checkIndNamesConsistency(crossTable, g$haplotypes)
    }, NA)
    # parent not used in the crossing table (no error)
     expect_error({
      checkIndNamesConsistency(crossTable[1:2,], g$haplotypes)
    }, NA)

    # additional individual in the crossing table
    crossTable_addInds <- crossTable
    crossTable_addInds[nrow(crossTable) + 1, ] <- list('doNotExist',
                                                       'F2_0001.0001',
                                                       'toto')
    expect_error({
      checkIndNamesConsistency(crossTable_addInds, g$haplotypes)
    }, paste('1 individuals are defined in the crossing',
             'table but not in the genoype file: doNotExist'))

    crossTable_addInds2 <- crossTable
    crossTable_addInds2[nrow(crossTable) + 1, ] <- list('F2_0001.0001',
                                                        'doNotExist2',
                                                        'toto')
     expect_error({
      checkIndNamesConsistency(crossTable_addInds2, g$haplotypes)
    }, paste('1 individuals are defined in the crossing',
             'table but not in the genoype file: doNotExist2'))
  })


})
