# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# test engine dependencies



capture_output({
  test_that("dependencies", {
    expect_true(compareVersion("1.5.7", as.character(packageVersion("gaston"))) != 1)
    expect_true(compareVersion("1.7.2", as.character(packageVersion("jsonlite"))) != 1)
    expect_true(compareVersion("0.3.0", as.character(packageVersion("manhattanly"))) != 1)
    expect_true(compareVersion("0.3.0", as.character(packageVersion("manhattanly"))) != 1)
  })
})
