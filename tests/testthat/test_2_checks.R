# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# unit test file for checks functions



capture_output({
  # checkAndFilterSNPcoord ----
  phasedGenoFile <- "../../data/geno/breedGame_phasedGeno.vcf.gz"
  g <- readPhasedGeno(phasedGenoFile)

  SNPcoordFiles <- c("../../data/SNPcoordinates/breedingGame_SNPcoord.csv")
  for (SNPcoordFile in SNPcoordFiles) {
    test_that(paste("checkAndFilterSNPcoord -", SNPcoordFile), {
      SNPcoord <- readSNPcoord(SNPcoordFile)

      # no error / no filter
      expect_no_error({
        filteredSNPcoord <- checkAndFilterSNPcoord(SNPcoord, g$SNPcoord)
      })
      expect_identical(filteredSNPcoord, SNPcoord)

      # additional SNP in user's SNP (should be removed)
      SNPcoord_additionalSNP <- SNPcoord
      SNPcoord_additionalSNP[nrow(SNPcoord) + 1, ] <- list(
        "chr1",
        42,
        "additionalSNP",
        0.00042
      )
      expect_warning(
        {
          filteredSNPcoord_addSNP <- checkAndFilterSNPcoord(
            SNPcoord_additionalSNP,
            g$SNPcoord
          )
        },
        paste(
          "The SNPs coordinate file have 1 SNPs not defined in the `.vcf`",
          "file, those SNP will not be considered: additionalSNP"
        )
      )
      expect_equal(filteredSNPcoord_addSNP, SNPcoord)


      # missing SNP in user's SNP (error)
      SNPcoord_missingSNP <- SNPcoord[-nrow(SNPcoord), ]
      err <- expect_engineError({
        checkAndFilterSNPcoord(SNPcoord_missingSNP, g$SNPcoord)
      })
    })
  }

  SNPcoordFiles <- c("../data/breedingGame_SNPcoord_missingPhysPos.csv",
                     "../data/breedingGame_SNPcoord_missingPhysPos_columns.csv")
  for (SNPcoordFile in SNPcoordFiles) {
    test_that(paste("checkAndFilterSNPcoord special cases -", SNPcoordFile), {
      SNPcoord <- readSNPcoord(SNPcoordFile)

      expect_warning({
        filteredSNPcoord <- checkAndFilterSNPcoord(SNPcoord, g$SNPcoord)
      })
      expect_true(!any(is.na(filteredSNPcoord$physPos)))
    })
  }

  SNPcoordFiles_error <- c("../data/breedingGame_SNPcoord_inconsistent_physPos_against_vcf.csv")
  for (SNPcoordFile in SNPcoordFiles_error) {
    test_that(paste("checkAndFilterSNPcoord error -", SNPcoordFile), {
      SNPcoord <- readSNPcoord(SNPcoordFile)
      err <- expect_engineError({
        filteredSNPcoord <- checkAndFilterSNPcoord(SNPcoord, g$SNPcoord)
      })
    })
  }





  # checkIndNamesConsistency ----
  crossTableFile <- "../../data/crossingTable/breedGame_crossTable.csv"
  phasedGenoFile <- "../../data/geno/breedGame_phasedGeno.vcf.gz"
  test_that("checkIndNamesConsistency", {
    g <- readPhasedGeno(phasedGenoFile)
    crossTable <- readCrossTable(crossTableFile)

    # no error
    expect_no_error({
      checkIndNamesConsistency(crossTable, g$haplotypes)
    })
    # parent not used in the crossing table (no error)
    expect_no_error({
      checkIndNamesConsistency(crossTable[1:2, ], g$haplotypes)
    })

    # additional individual in the crossing table
    crossTable_addInds <- crossTable
    crossTable_addInds[nrow(crossTable) + 1, ] <- list(
      "F1_additional_ind",
      "F2_0001.0001",
      "toto"
    )
    err <- expect_engineError({
      checkIndNamesConsistency(crossTable_addInds, g$haplotypes)
    })

    crossTable_addInds2 <- crossTable
    crossTable_addInds2[nrow(crossTable) + 1, ] <- list(
      "F2_0001.0001",
      "F1_additional_ind2",
      "toto"
    )
    err <- expect_engineError({
      checkIndNamesConsistency(crossTable_addInds2, g$haplotypes)
    })
  })
})
