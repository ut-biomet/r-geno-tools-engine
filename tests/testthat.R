library(testthat)

# load functions
invisible(
  sapply(FUN = source,
         X = list.files("src", pattern = ".R$",full.names = T))
)

# run tests
test_dir("tests/testthat",
         stop_on_failure = TRUE,
         stop_on_warning = TRUE)

# clean testOutput dir
f <- list.files("tests/testthat/testOutput/", full.names = TRUE)
lapply(f, file.remove)
