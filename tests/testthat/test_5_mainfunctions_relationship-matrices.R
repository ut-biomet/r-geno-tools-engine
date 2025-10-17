capture.output({
  # calc_pedRelMat ----
  pedFiles <- c(
    "../../data/pedigree/testPedData_char.csv",
    "../../data/pedigree/testPedData_num.csv",
    "../../data/pedigree/testPedData_missFounder.csv"
  )
  formats <- c("csv", "json")
  for (file in pedFiles) {
    for (format in formats) {
      test_that(paste("calc_pedRelMat", format, basename(file)), {
        expect_no_error({
          suppressWarnings({
            relMat_results <- calc_pedRelMat(
              pedFile = file,
              pedUrl = NULL,
              header = TRUE,
              unknown_string = "",
              outFormat = format
            )
          })
        })

        expect_true(class(relMat_results) == "list")
        expect_identical(
          names(relMat_results),
          c("relMat", "metadata", "file")
        )
        expect_true(is.matrix(relMat_results$relMat))
        expect_identical(
          names(relMat_results$metadata),
          c("info", "date", "nInds", "pedFP")
        )
        expect_true(file.exists(relMat_results$file))
      })
    }
  }

  # calc_genoRelMat ----
  genoFiles <- c("../../data/geno/breedGame_geno.vcf.gz")
  formats <- c("csv", "json")
  for (format in formats) {
    for (file in genoFiles) {
      test_that(paste("calc_genoRelMat", format, basename(file)), {
        expect_no_error({
          relMat_results <- calc_genoRelMat(
            genoFile = file,
            genoUrl = NULL,
            n_markers = 2000,
            outFormat = format
          )
        })

        expect_true(class(relMat_results) == "list")
        expect_identical(
          names(relMat_results),
          c("relMat", "metadata", "file")
        )
        expect_true(is.matrix(relMat_results$relMat))
        expect_identical(
          names(relMat_results$metadata),
          c("info", "date", "nInds", "genoFP")
        )
        expect_true(file.exists(relMat_results$file))
      })
    }
  }



  # calc_combinedRelMat ----
  rm_files <- list(
    breedGameFiles = list(
      ped_rm = c("../../data/results/breedGame_pedRelMat.csv"),
      geno_rm = c("../../data/results/breedGame_genoRelMat.csv")
    )
  )
  formats <- c("csv", "json")
  methods <- c("Legarra", "Martini")
  tau <- list(Legarra = NULL, Martini = 0.4)
  omega <- list(Legarra = NULL, Martini = 0.7)
  for (format in formats) {
    for (method in methods) {
      for (filesName in names(rm_files)) {
        test_that(paste("calc_genoRelMat", format, method, filesName), {
          files <- rm_files[[filesName]]
          expect_no_error({
            relMat_results <- calc_combinedRelMat(
              pedRelMatFile = files$ped_rm,
              genoRelMatFile = files$geno_rm,
              method = method,
              tau = tau[[method]],
              omega = omega[[method]],
              outFormat = format
            )
          })

          expect_true(class(relMat_results) == "list")
          expect_identical(
            names(relMat_results),
            c("relMat", "metadata", "file")
          )
          expect_true(is.matrix(relMat_results$relMat))
          expect_identical(
            names(relMat_results$metadata),
            c(
              "info", "date", "nInds",
              "geno_relMatFP", "ped_relMatFP"
            )
          )
          expect_true(file.exists(relMat_results$file))
        })
      }
    }
  }



  # draw_relMat ----
  relMatFiles <- c(
    "../../data/results/pedigreeRelationship.csv",
    "../../data/results/pedigreeRelationship.json"
  )
  for (file in relMatFiles) {
    for (inter in c(TRUE, FALSE)) {
      test_that(paste("draw_relHeatmap inter:", inter, basename(file)), {
        ext <- ".html"
        if (!inter) {
          ext <- ".png"
        }
        tmpF <- tempfile(fileext = ext)
        expect_no_error({
          pedNet <- draw_relHeatmap(
            relMatFile = file,
            relMatUrl = NULL,
            interactive = inter,
            outFile = tmpF
          )
        })
        expect_true(file.exists(tmpF))
        expect_true(file.info(tmpF)$size > 0)
      })
    }
  }
})
