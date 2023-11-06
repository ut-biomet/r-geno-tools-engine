# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# test engine dependencies



capture_output({
  test_that("dependencies", {
    expect_true(compareVersion("1.5.7", as.character(packageVersion("gaston"))) != 1)
    expect_true(compareVersion("0.6.27", as.character(packageVersion("digest"))) != 1)
    expect_true(compareVersion("1.7.2", as.character(packageVersion("jsonlite"))) != 1)
    expect_true(compareVersion("0.3.0", as.character(packageVersion("manhattanly"))) != 1)
    expect_true(compareVersion("4.0.4", as.character(packageVersion("tools"))) != 1)
    expect_true(compareVersion("2.5.0", as.character(packageVersion("R6"))) != 1)
    expect_true(compareVersion("4.9.3", as.character(packageVersion("plotly"))) != 1)
    expect_true(compareVersion("0.4", as.character(packageVersion("networkD3"))) != 1)
    expect_true(compareVersion("0.5.0", as.character(packageVersion("ellipse"))) != 1)
  })
})
