capture.output({
  # crossingSimulation ----
  crossigSimDataList <- list(
    completeDta = list(
      phasedGenoFile = "../../data/geno/breedGame_phasedGeno.vcf.gz",
      SNPcoordFile = "../../data/SNPcoordinates/breedingGame_SNPcoord.csv",
      crossTableFile = "../../data/crossingTable/breedGame_crossTable.csv",
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    ),
    missingChr = list(
      phasedGenoFile = "../../data/geno/breedGame_phasedGeno.vcf.gz",
      SNPcoordFile = "../data/breedingGame_SNPcoord_missingChr.csv",
      crossTableFile = "../../data/crossingTable/breedGame_crossTable.csv",
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    ),
    missingSNPid = list(
      phasedGenoFile = "../data/breedGame_phasedGeno_missingIds.vcf.gz",
      SNPcoordFile = "../data/breedingGame_SNPcoord_missingIds.csv",
      crossTableFile = "../../data/crossingTable/breedGame_crossTable.csv",
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    ),
    missingPhysPos = list(
      phasedGenoFile = "../data/breedGame_phasedGeno_missingPos.vcf.gz",
      SNPcoordFile = "../data/breedingGame_SNPcoord_missingPhysPos.csv",
      crossTableFile = "../../data/crossingTable/breedGame_crossTable.csv",
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    ),
    missingPhysPosColumns = list(
      phasedGenoFile = "../data/breedGame_phasedGeno_missingPos.vcf.gz",
      SNPcoordFile = "../data/breedingGame_SNPcoord_missingPhysPos_columns.csv",
      crossTableFile = "../../data/crossingTable/breedGame_crossTable.csv",
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    )
  )

  for (i in seq_along(crossigSimDataList)) {
    data <- crossigSimDataList[[i]]
    phasedGenoFile <- data$phasedGenoFile
    SNPcoordFile <- data$SNPcoordFile
    crossTableFile <- data$crossTableFile
    nCross <- data$nCross
    outFile <- data$outFile

    test_that(paste("crossingSimulation -", names(crossigSimDataList)[i]), {
      if (grepl(pattern = "missingPhysPos|missingChr", names(crossigSimDataList)[i])) {
        expect_warning({
          createdFile <- crossingSimulation(
            genoFile = phasedGenoFile,
            crossTableFile = crossTableFile,
            SNPcoordFile = SNPcoordFile,
            nCross = nCross,
            outFile = outFile
          )
        })
      } else {
        expect_no_error({
          createdFile <- crossingSimulation(
            genoFile = phasedGenoFile,
            crossTableFile = crossTableFile,
            SNPcoordFile = SNPcoordFile,
            nCross = nCross,
            outFile = outFile
          )
        })
      }
      expect_equal(createdFile, outFile)
      expect_no_error({
        createdGeno <- readGenoData(createdFile)
      })
      expect_no_error({
        createdPhasedGeno <- readPhasedGeno(createdFile)
      })

      # individual names
      crossTable <- readCrossTable(crossTableFile)
      expect_equal(nrow(createdGeno@ped), nrow(crossTable) * nCross)
      createdIndBaseName <- unique(gsub("-[0-9]*$", "", createdGeno@ped$id))
      expect_equal(sort(createdIndBaseName), sort(crossTable$names))

      # SNPs: check original and new vcf files' snp position are the same
      g <- readPhasedGeno(phasedGenoFile)
      snpPhysPos_newVcf <- createdPhasedGeno$SNPcoord[order(createdPhasedGeno$SNPcoord$SNPid), ]
      snpPhysPos_ref <- g$SNPcoord[order(g$SNPcoord$SNPid), ]
      expect_equal(snpPhysPos_newVcf$SNPid, snpPhysPos_ref$SNPid)
      expect_equal(snpPhysPos_newVcf$physPos, snpPhysPos_ref$physPos)

      if (grepl("missingPhysPos", names(crossigSimDataList)[i])) {
        expect_true(all(is.na(createdPhasedGeno$SNPcoord$physPos)))
      }
    })
  }

  test_that("crossingSimulation with downloaded files", {
    as_url <- function(file) {
      file <- normalizePath(file)
      file <- paste0("file://", file)
      file
    }
    expect_warning(
      {
        createdFileWithoutChrInfo <- crossingSimulation(
          genoUrl = as_url(phasedGenoFile),
          crossTableUrl = as_url(crossTableFile),
          SNPcoordUrl = as_url(SNPcoordFile),
          nCross = nCross,
          outFile = outFile
        )
      },
      'output "file" already exists. This file will be overwritten.'
    )
  })
  unlink(outFile)

  inconsistentSNPFile <- "../data/inconsistent_SNPcoord_2.csv"
  test_that("crossingSimulation inconsistent SNPs", {
    err <- expect_warning({
      createdFile <- crossingSimulation(
        genoFile = phasedGenoFile,
        crossTableFile = crossTableFile,
        SNPcoordFile = inconsistentSNPFile,
        nCross = nCross,
        outFile = outFile
      )
    })
  })
})
