# Author: Julien Diot juliendiot@ut-biomet.org
# 2024 The University of Tokyo
#
# Description:
# relationship matrix tests.

library(testthat)

capture_output({
  test_that("get_model_metrics", {
    for (n in c(1, 2, 10, 100)) {
      mock_actual_values <- rnorm(n)
      mock_predicted_values <- mock_actual_values + rnorm(n)

      expect_no_error({
        metrics <- get_model_metrics(mock_actual_values, mock_predicted_values)
      })

      expect_true(is.list(metrics))
      expect_false(is.null(names(metrics)))
      expect_false(any(names(metrics) == ""))
      expect_false(any(duplicated(names(metrics))))

      for (metric in metrics) {
        expect_true(is.numeric(metric))
      }
    }
  })

  n_inds <- 100
  n_markers <- 1000

  ind_names <- paste0("ind_", seq(1:n_inds))
  markers_names <- paste0("snp_", seq(1:n_markers))


  geno <- matrix(
    sample(c(0, 1, 2), n_inds * n_markers, replace = TRUE),
    ncol = n_markers
  )
  row.names(geno) <- ind_names
  colnames(geno) <- markers_names

  mock_marker_eff_add <- rnorm(n_markers)
  gv_additive <- geno %*% mock_marker_eff_add
  pheno <- data.frame(
    pheno_trait = gv_additive + rnorm(n_inds, 0, sqrt(var(gv_additive))),
    row.names = ind_names
  )

  geno <- gaston::as.bed.matrix(geno)



  data_list <- list(
    list(
      geno = geno,
      pheno = pheno
    )
  )

  for (data in data_list) {
    n_inds <- nrow(data$pheno)
    n_markers <- ncol(data$geno)

    pheno_vector <- data$pheno[, 1]
    rel_mat_add <- calc_additive_rel_mat(data$geno, standardized = TRUE)$rel_mat
    rel_mat_dom <- calc_dominance_rel_mat(data$geno, standardized = TRUE)$rel_mat
    K_add_only <- list(
      rel_mat_add
    )
    K_add_dom <- list(
      rel_mat_add,
      rel_mat_dom
    )

    fit_functions <- list(
      fit_with_gaston = fit_with_gaston,
      fit_with_rainbowr = fit_with_rainbowr
    )

    for (fit_function_name in names(fit_functions)) {
      fit_function <- fit_functions[[fit_function_name]]
      test_that(fit_function_name, {
        expect_no_error({
          model_add_only <- fit_function(pheno_vector, K_add_only)
        })
        expect_no_error({
          model_add_dom <- fit_function(pheno_vector, K_add_dom)
        })

        for (model in list(model_add_only, model_add_dom)) {
          expect_true(is.list(model))
          expect_false(is.null(names(model)))
          expect_identical(names(model), c("logL", "blups", "intercept"))

          expect_true(is.numeric(model$logL))
          expect_equal(length(model$logL), 1)

          expect_true(is.numeric(model$intercept))
          expect_equal(length(model$intercept), 1)

          expect_true(is.matrix(model$blups))
          expect_true(is.numeric(model$blups))
          expect_equal(nrow(model$blups), n_inds)
        }
        expect_equal(ncol(model_add_only$blups), 1)
        expect_equal(ncol(model_add_dom$blups), 2)
      })
    }



    test_that("calc_mark_eff", {
      expect_no_error({
        estim_mark_eff <- calc_mark_eff(
          gaston::as.matrix(data$geno),
          rel_mat_add,
          rnorm(n_inds)
        )
      })
      expect_true(is.vector(estim_mark_eff))
      expect_true(is.numeric(estim_mark_eff))
      expect_equal(length(estim_mark_eff), n_markers)
      expect_identical(
        names(estim_mark_eff),
        colnames(gaston::as.matrix(data$geno))
      )
    })



    test_that("train_gs_model", {
      # Check marker effect are for markers encoded in 0, 1, 2 or (dom 0,1,0)
      expect_no_error({
        model_add <- train_gs_model(data$pheno, data$geno, with_dominance = FALSE)
      })
      expect_no_error({
        model_dom <- train_gs_model(data$pheno, data$geno, with_dominance = TRUE)
      })

      for (model in list(model_add, model_dom)) {
        expect_true(is.list(model))
        expect_identical(names(model), c("intercept", "eff"))

        expect_true(is.numeric(model$intercept))
        expect_equal(length(model$intercept), 1)

        expect_true(is.data.frame(model$eff))
        expect_identical(colnames(model$eff), c("additive", "dominance"))

        expect_identical(
          row.names(model$eff),
          colnames(gaston::as.matrix(data$geno))
        )
      }

      expect_true(all(is.na(model_add$eff$dominance)))
    })




    test_that("predict_gs_model", {
      model_add <- train_gs_model(data$pheno, data$geno, with_dominance = FALSE)
      model_dom <- train_gs_model(data$pheno, data$geno, with_dominance = TRUE)

      for (model in list(model_add, model_dom)) {
        expect_no_error({
          predictions <- predict_gs_model(data$geno, model)
        })

        expect_true(is.data.frame(predictions))
        expect_equal(ncol(predictions), 1)
        expect_equal(nrow(predictions), nrow(data$geno))
        expect_false(any(is.na(predictions[, 1])))
      }

      expect_false(isTRUE(all.equal(
        predict_gs_model(data$geno, model_add),
        predict_gs_model(data$geno, model_dom)
      )))
    })





    n_folds <- c(2, 10)
    n_repetitions <- c(1, 2, 3)
    with_dominance <- c(TRUE, FALSE)
    test_cases <- expand.grid(n_folds, n_repetitions, with_dominance)
    colnames(test_cases) <- c("n_folds", "n_repetitions", "with_dominance")

    apply(test_cases, 1, function(test_case) {
      test_name <- paste0(
        "cross_validation_evaluation, n_folds = ",
        test_case["n_folds"],
        ", n_repetitions = ",
        test_case["n_repetitions"],
        ", with_dominance = ",
        test_case["with_dominance"]
      )
      test_that(test_name, {
        expect_no_error({
          evaluation <- cross_validation_evaluation(
            data$pheno,
            data$geno,
            with_dominance = test_case["with_dominance"],
            n_folds = test_case["n_folds"],
            n_repetitions = test_case["n_repetitions"]
          )
        })
        expect_true(is.list(evaluation))
        expect_identical(
          names(evaluation),
          c("predictions", "metrics")
        )
        expect_true(is.data.frame(evaluation$predictions))
        expect_identical(
          colnames(evaluation$predictions),
          c("ind", "actual", "predicted", "fold", "repetition"),
        )
        expect_equal(
          nrow(evaluation$predictions),
          as.numeric(test_case["n_repetitions"]) * nrow(data$pheno),
        )
        expect_false(any(is.na(evaluation$predictions)))



        expect_true(is.data.frame(evaluation$metrics))
        expect_equal(
          nrow(evaluation$metrics),
          as.numeric(test_case["n_repetitions"])
        )
        expect_true(all(c("repetition") %in% colnames(evaluation$metrics)))
        expect_false(any(is.na(evaluation$metrics)))
      })
    })
  }




  test_that("train_gs_model returns result encoded in allele dose", {
    # Test the markers effects are returned for the genotypes
    # encoded in allele doses.
    # I can't think of a deterministic method to test that, instead I
    # compare the predictions made on standardized genotype against those
    # on allele dose genotype but this may not be accurate.
    pheno_train <- readPhenoData("../../data/genomic_selection/pheno_train.csv")
    geno_train <- readGenoData("../../data/genomic_selection/geno_G1.vcf.gz")
    model <- train_gs_model(pheno_train, geno_train, with_dominance = FALSE)

    geno_test <- readGenoData("../../data/genomic_selection/geno_G2.vcf.gz")
    geno_test_raw <- calc_additive_geno(geno_test,
      standardized = FALSE
    )
    geno_test_std <- calc_additive_geno(geno_test,
      standardized = TRUE
    )
    pheno_test <- readPhenoData("../../data/genomic_selection/pheno_test.csv")
    pheno_test$gv <- pheno_test$gt + 100

    pred_with_std_geno <- (
      model$intercept
        + geno_test_std %*% model$eff$additive
    )
    pred_with_raw_geno <- (
      model$intercept
        + geno_test_raw %*% model$eff$additive
    )

    metric_std <- get_model_metrics(pheno_test$gv, pred_with_std_geno)
    metric_raw <- get_model_metrics(pheno_test$gv, pred_with_raw_geno)

    expect_true(metric_raw$rmse < metric_std$rmse)
    expect_true(metric_raw$corel_pearson > metric_std$corel_pearson)
    expect_true(metric_raw$corel_spearman > metric_std$corel_spearman)
    expect_true(metric_raw$r2 > metric_std$r2)
  })
})
