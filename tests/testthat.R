library(testthat)

# load functions
invisible(
  sapply(
    FUN = source,
    X = list.files("src", pattern = ".R$", full.names = T)
  )
)
source("tests/engine_error_expectations.R")

UNTESTED_ERROR_CODES <- NULL

# run tests
# test_file("tests/testthat/test_0_dependencies.R",
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_1_readWriteData.R",
#   reporter = StopReporter$new(),
#   stop_on_failure = TRUE,
#   stop_on_warning = FALSE
# )

# test_file("tests/testthat/test_2_gwas.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_2_progBlupEstim.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_2_relationshipMatrix.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_2_checks.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_2_gs_models.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_3_crossingSimulation.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_3_filtering.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_4_plots.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_5_mainfunctions_LD-plot.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_5_mainfunctions_crossing-simulation.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_5_mainfunctions_gs.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_5_mainfunctions_gwas.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_5_mainfunctions_pedigree-network.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_5_mainfunctions_progBlups.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_5_mainfunctions_relationship-matrices.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

# test_file("tests/testthat/test_6_engineCommands.R",
#           reporter = StopReporter$new(),
#           stop_on_failure = TRUE,
#           stop_on_warning = FALSE)

test_dir("tests/testthat",
  stop_on_failure = TRUE,
  stop_on_warning = FALSE
)

if (length(UNTESTED_ERROR_CODES) != 0) {
  warning(paste0(
    "Error codes `",
    paste(UNTESTED_ERROR_CODES, collapse = "`, `"),
    "` have not been tested"
  ))
}
