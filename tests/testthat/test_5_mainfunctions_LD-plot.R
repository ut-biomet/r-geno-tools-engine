capture.output({
  # Draw LD plot ----
  test_that("Draw LD plot", {
    expect_no_error({
      imgFile <- draw_ldPlot(
        genoFile = "../../data/geno/testMarkerData01.vcf.gz",
        genoUrl = NULL,
        from = 42,
        to = 62,
        outFile = tempfile(fileext = ".png")
      )
    })
  })


  # Test draw_ldPlot() with wrong parameters:
  goodParams <- list(
    genoFile = "../../data/geno/testMarkerData01.vcf.gz",
    genoUrl = NULL,
    from = 42,
    to = 62,
    outFile = tempfile(fileext = ".png")
  )

  wrongParamsL <- list(
    genoFile = c("doNotExist", NA),
    genoUrl = c("doNotExist", NA),
    from = c("42", 42.1),
    to = c("50", 52.1),
    outFile = list(c("f1", "f2"))
  )

  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i + 1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("draw_ldPlot, WrongParams", p, i, sep = "-")
      test_that(testName, {
        err <- expect_engineError({
          draw_ldPlot(
            genoFile = params[["genoFile"]],
            genoUrl = params[["genoUrl"]],
            from = params[["from"]],
            to = params[["to"]],
            outFile = params[["outFile"]]
          )
        })
      })
    }
  }
})
