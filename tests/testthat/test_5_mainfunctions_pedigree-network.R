capture.output({
  # draw_pedNetrowk ----
  pedFiles <- c(
    "../../data/pedigree/testPedData_char.csv",
    "../../data/pedigree/testPedData_num.csv",
    "../../data/pedigree/testPedData_missFounder.csv"
  )
  for (file in pedFiles) {
    test_that(paste("draw_pedNetwork", basename(file)), {
      tmpF <- tempfile(fileext = ".html")
      expect_no_error({
        suppressWarnings({
          pedNet <- draw_pedNetwork(
            pedFile = file,
            pedUrl = NULL,
            header = TRUE,
            unknown_string = "",
            outFile = tmpF
          )
        })
      })
      expect_identical(
        class(pedNet),
        c("forceNetwork", "htmlwidget")
      )
      expect_true(file.exists(tmpF))
      expect_true(file.info(tmpF)$size > 0)
    })
  }
})
