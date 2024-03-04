.PHONY: tests

deps:
	Rscript --vanilla -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')"
	Rscript --vanilla -e "if (!requireNamespace('renv', quietly = TRUE)) remotes::install_github('rstudio/renv', '$(renv_v)')"
	Rscript -e "renv::restore()"


tests:
	Rscript --vanilla ./tests/testthat.R
