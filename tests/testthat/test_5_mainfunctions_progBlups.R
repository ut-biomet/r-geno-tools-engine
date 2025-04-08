capture.output({
  # calc_progenyBlupEstimation ----
  progBlupDataList <- list(
    completeDta = list(
      genoFile = "../../data/geno/breedGame_phasedGeno.vcf.gz",
      crossTableFile = "../../data/crossingTable/breedGame_small_crossTable.csv",
      SNPcoordFile = "../../data/SNPcoordinates/breedingGame_SNPcoord.csv",
      markerEffectsFile = "../../data/markerEffects/breedGame_markerEffects.csv",
      outFile = tempfile(fileext = ".json")
    ),
    missingChr = list(
      genoFile = "../../data/geno/breedGame_phasedGeno.vcf.gz",
      crossTableFile = "../../data/crossingTable/breedGame_small_crossTable.csv",
      SNPcoordFile = "../data/breedingGame_SNPcoord_missingChr.csv",
      markerEffectsFile = "../../data/markerEffects/breedGame_markerEffects.csv",
      outFile = tempfile(fileext = ".json")
    ),
    missingSNPid = list(
      genoFile = "../data/breedGame_phasedGeno_missingIds.vcf.gz",
      SNPcoordFile = "../data/breedingGame_SNPcoord_missingIds.csv",
      crossTableFile = "../../data/crossingTable/breedGame_small_crossTable.csv",
      markerEffectsFile = "../data/breedGame_markerEffects_imputedSNPid.csv",
      outFile = tempfile(fileext = ".json")
    ),
    missingPhysPos = list(
      genoFile = "../data/breedGame_phasedGeno_missingPos.vcf.gz",
      SNPcoordFile = "../data/breedingGame_SNPcoord_missingPhysPos.csv",
      crossTableFile = "../../data/crossingTable/breedGame_small_crossTable.csv",
      markerEffectsFile = "../../data/markerEffects/breedGame_markerEffects.csv",
      outFile = tempfile(fileext = ".json")
    ),
    missingPhysPosColumns = list(
      genoFile = "../data/breedGame_phasedGeno_missingPos.vcf.gz",
      SNPcoordFile = "../data/breedingGame_SNPcoord_missingPhysPos_columns.csv",
      crossTableFile = "../../data/crossingTable/breedGame_small_crossTable.csv",
      markerEffectsFile = "../../data/markerEffects/breedGame_markerEffects.csv",
      outFile = tempfile(fileext = ".json")
    ),
    completeDta_severalTraits = list(
      genoFile = "../../data/geno/breedGame_phasedGeno.vcf.gz",
      crossTableFile = "../../data/crossingTable/breedGame_small_crossTable.csv",
      SNPcoordFile = "../../data/SNPcoordinates/breedingGame_SNPcoord.csv",
      markerEffectsFile = "../../data/markerEffects/breedGame_markerEffects_2traits.json",
      outFile = tempfile(fileext = ".json")
    ),
  )

  for (i in seq_along(progBlupDataList)) {
    genoFile <- progBlupDataList[[i]]$genoFile
    crossTableFile <- progBlupDataList[[i]]$crossTableFile
    SNPcoordFile <- progBlupDataList[[i]]$SNPcoordFile
    markerEffectsFile <- progBlupDataList[[i]]$markerEffectsFile
    outFile <- progBlupDataList[[i]]$outFile

    test_that(paste("calc_progenyBlupEstimation -", names(progBlupDataList)[i]), {
      if (grepl(pattern = "missingPhysPos|missingChr", names(progBlupDataList)[i])) {
        expect_warning({
          projBlups_list <- calc_progenyBlupEstimation(
            genoFile = genoFile,
            crossTableFile = crossTableFile,
            SNPcoordFile = SNPcoordFile,
            markerEffectsFile = markerEffectsFile,
            outFile = outFile
          )
        })
      } else {
        expect_no_error({
          projBlups_list <- calc_progenyBlupEstimation(
            genoFile = genoFile,
            crossTableFile = crossTableFile,
            SNPcoordFile = SNPcoordFile,
            markerEffectsFile = markerEffectsFile,
            outFile = outFile
          )
        })
      }
      expect_is(projBlups_list, "list")
      lapply(projBlups_list, function(projBlups) {
        expect_is(projBlups, "list")
        expect_identical(
          names(projBlups),
          c("ind1", "ind2", "blup_exp", "blup_var", "cov")
        )
      })
    })
  }


  # draw_progBlupsPlot ----
  progEstimFiles <- c(
    "../../data/results/progenyBlupEstimation.json",
    "../data/progenyBlupEstimation_oldVersion.json"
  )
  for (file in progEstimFiles) {
    test_that(paste("draw_progBlupsPlot", basename(file)), {
      tmpF <- tempfile(fileext = ".html")
      expect_no_error({
        p1 <- draw_progBlupsPlot(
          progEstimFile = file,
          outFile = tmpF
        )
        unlink(tmpF)
      })

      for (s in c("asc", "dec")) {
        expect_no_error({
          p1 <- draw_progBlupsPlot(
            progEstimFile = file,
            sorting = s,
            outFile = tmpF
          )
          unlink(tmpF)
        })
      }
      expect_warning(
        {
          p1 <- draw_progBlupsPlot(
            progEstimFile = file,
            sorting = "??????",
            outFile = tmpF
          )
          unlink(tmpF)
        },
        "Sorting method not recognised, x axis will be sorted alphabetically. Recognised methods are: alpha, asc, dec"
      )
    })
  }

  progEstimFiles_severalTraits <- "../../data/results/progenyBlupEstimation_2traits.json"
  for (file in progEstimFiles_severalTraits) {
    test_that(paste("draw_progBlupsPlot several traits", basename(file)), {
      tmpF <- tempfile(fileext = ".html")
      expect_no_error({
        p1 <- draw_progBlupsPlot(
          progEstimFile = file,
          trait = "trait1",
          outFile = tmpF
        )
        unlink(tmpF)
      })

      for (s in c("asc", "dec")) {
        expect_no_error({
          p1 <- draw_progBlupsPlot(
            progEstimFile = file,
            sorting = s,
            trait = "trait1",
            outFile = tmpF
          )
          unlink(tmpF)
        })
      }
      expect_warning(
        {
          p1 <- draw_progBlupsPlot(
            progEstimFile = file,
            sorting = "??????",
            trait = "trait1",
            outFile = tmpF
          )
          unlink(tmpF)
        },
        "Sorting method not recognised, x axis will be sorted alphabetically. Recognised methods are: alpha, asc, dec"
      )
    })
  }


  # draw_progBlupsPlot_2traits ----
  for (file in progEstimFiles_severalTraits) {
    test_that(paste("draw_progBlupsPlot_2traits several traits", basename(file)), {
      tmpF <- tempfile(fileext = ".html")
      expect_no_error({
        p1 <- draw_progBlupsPlot_2traits(
          progEstimFile = file,
          x_trait = "trait1",
          y_trait = "trait2",
          confidenceLevel = 0.95,
          x_suffix = "",
          y_suffix = "",
          ellipses_npoints = 100,
          outFile = tmpF
        )
        unlink(tmpF)
      })
    })
  }
})
