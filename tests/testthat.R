library(testthat)

# load functions
invisible(
  sapply(FUN = source,
         X = list.files("src", pattern = ".R$", full.names = T))
)


# run tests
# test_file("tests/testthat/test_0_dependencies.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)
# test_file("tests/testthat/test_1_readWriteData.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)
# test_file("tests/testthat/test_2_gwas.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)
# test_file("tests/testthat/test_2_relationshipMatrix.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)
# test_file("tests/testthat/test_2_checks.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)
# test_file("tests/testthat/test_3_crossingSimulation.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)
# test_file("tests/testthat/test_3_filtering.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)
# test_file("tests/testthat/test_4_plots.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)
# test_file("tests/testthat/test_5_mainfunctions.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)
# test_file("tests/testthat/test_6_engineCommands.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

test_dir("tests/testthat",
         stop_on_failure = TRUE,
         stop_on_warning = FALSE)

# clean testOutput dir
f <- list.files("tests/testthat/testOutput",
                full.names = TRUE,
                include.dirs = TRUE)
invisible(lapply(f, unlink, recursive = TRUE))
