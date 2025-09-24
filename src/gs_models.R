# Author: Julien Diot juliendiot@ut-biomet.org
# 2024 The University of Tokyo
#
# Description:
# Functions for builing/predicting GS models



#' Fit a GS model using `gaston` package
#'
#' @param pheno **1 dimensional** vector of phenotypic values
#' @param K list of 1 or 2 variance covariance matrices to consider for the
#' random effect.
#'
#' @details This is a wrapper around `gaston::lmm.aireml()` see
#' `?gaston::lmm.aireml()` for more information
#'
#' @return list of 3 elements
#' - `logL`: Value of log-likelihood of the model (cf. `?gaston::lmm.aireml()`)
#' - `blups`: `length(pheno)` by `length(K)` matrix of the blups of the
#' individuals (one random effects per column).
#' - `intercept`: value of the intercept of the model
fit_with_gaston <- function(pheno, K) {
  mod_gaston <- gaston::lmm.aireml(
    Y = pheno,
    K = K,
    verbose = FALSE
  )

  if (length(K) == 1) {
    blups <- matrix(mod_gaston$BLUP_omega, ncol = 1)
  } else {
    blups <- matrix(NA, nrow = length(pheno), ncol = length(K))
    for (i in seq(length(K))) {
      blups[, i] <- mod_gaston$tau[i] * K[[i]] %*% mod_gaston$Py
    }
  }

  return(
    list(
      logL = mod_gaston$logL,
      blups = blups,
      intercept = mod_gaston$BLUP_beta
    )
  )
}

#' Fit a GS model using `RAINBOWR` package
#'
#' @param pheno **1 dimensional** vector of phenotypic values
#' @param K list of 1 or 2 variance covariance matrices used for the random
#' effect.
#'
#' @details This is a wrapper around `RAINBOWR::EM3.general()` see
#' `?RAINBOWR::EM3.general()` for more information
#'
#' @return list of 3 elements
#' - `logL`: Value of log-likelihood of the model (cf. `?RAINBOWR::EM3.general()`)
#' - `blups`: `length(pheno)` by `length(K)` matrix of the blups of the
#' individuals (one random effects per column).
#' - `intercept`: value of the intercept of the model
fit_with_rainbowr <- function(pheno, K) {
  ZETA <- lapply(K, function(k) {
    list(
      Z = diag(nrow(k)),
      K = k
    )
  })

  mod_rainbowr <- RAINBOWR::EM3.general(
    package = "RAINBOWR",
    y = pheno,
    ZETA = ZETA
  )

  blups <- matrix(mod_rainbowr$u.each, byrow = FALSE, ncol = length(K))

  return(
    list(
      logL = mod_rainbowr$LL,
      blups = blups,
      intercept = mod_rainbowr$beta
    )
  )
}



#' Calculate markers effects from the BLUPS, the genetic matrix and the
#' variance covariance matrix of the blups (ie. genetic relationship matrices)
#'
#' @param geno_mat `n_inds` by `n_marker` genomic matrix
#' @param rel_mat `n_inds` by `n_inds` relationship matrix associated to
#' `geno_mat`
#' @param blups `n_inds` vector of the individuals blups
#' effect.
#'
#' @details Calculation is made by
#' \deqn{(\frac{X}{n_{mark}})(\frac{X}{n_{mark}})'Z^{-1}B}
#'
#' with:
#' - $X$ the genetic matrix
#' - $n_{mark}$ the number of markers (ie. number of columns of $X$)
#' - $Z$ the relationship matrix of $X$
#' - $B$ the blups of the individuals
#'
#' @return vector of the estimated markers effects
calc_mark_eff <- function(geno_mat, rel_mat, blups) {
  mark_eff <- crossprod(geno_mat / ncol(geno_mat), solve(rel_mat)) %*% blups
  # mark_eff[is.na(mark_eff)] <- 0
  return(t(mark_eff)[1, ])
}




#' Train GS model optionally with dominance
#'
#' @param pheno 1 column data.frame of the phenotype
#' @param geno bed.matrix returned by `readGenoData()`
#' @param with_dominance (default = FALSE) control if the model should include
#' dominance effects
#'
#' @details
#' IMPORTANT ! It will return the markers effects for the genotypes encoded:
#' - **in allele dose (0, 1, 2) for the additive effects**
#' - **as (0, 1, 0) for the dominance effects**
#'
#' @return list of 2 elements:
#' - `intercept`: the intercept of the model
#' - `eff`: data.frame of 2 columns `additive` and `dominance` with
#' their corresponding markers effects
train_gs_model <- function(pheno, geno, with_dominance = FALSE) {
  stopifnot(identical(row.names(pheno), geno@ped$id))
  stopifnot(all(!is.na(pheno[, 1])))
  stopifnot(all(geno@snps$callrate == 1))

  # remove monomorphic markers
  geno <- gaston::select.snps(geno, maf != 0)

  relationship_matrices <- list()

  additive_std <- calc_additive_rel_mat(geno, standardized = TRUE)
  relationship_matrices[[1]] <- additive_std$rel_mat

  if (with_dominance) {
    dominance_std <- calc_dominance_rel_mat(geno, standardized = TRUE)
    dominance_not_std <- calc_dominance_rel_mat(geno, standardized = FALSE)
    check_dominance_model_is_applicable(dominance_not_std)

    relationship_matrices[[2]] <- dominance_std$rel_mat
  }

  model <- fit_with_gaston(pheno[[1]], relationship_matrices)

  if (is.nan(model$logL) || is.infinite(model$logL)) {
    warning(
      "Failed optimisation with gaston for GS model fitting.",
      "Fallback to optimisation with RAINBOWR."
    )
    model <- fit_with_rainbowr(pheno[[1]], relationship_matrices)
  }

  estim_mark_eff_add_std <- calc_mark_eff(
    additive_std$geno_mat,
    additive_std$rel_mat,
    model$blups[, 1]
  )
  # valid only if estim_mark_eff_add_std have been calculated with standardized
  # genetic matrix !!!!
  estim_mark_eff_add <- estim_mark_eff_add_std / sqrt(2 * geno@p * (1 - geno@p))
  intercept_adjustment_add <- sum(estim_mark_eff_add * 2 * geno@p)

  estim_mark_eff_dom <- rep(NA, length(estim_mark_eff_add))
  names(estim_mark_eff_dom) <- names(estim_mark_eff_add)
  if (with_dominance) {
    estim_mark_eff_dom <- calc_mark_eff(
      dominance_not_std$geno_mat,
      dominance_not_std$rel_mat,
      model$blups[, 2]
    )
  }

  eff <- data.frame(
    row.names = geno@snps$id
  )
  eff$additive <- estim_mark_eff_add[row.names(eff)]
  eff$dominance <- estim_mark_eff_dom[row.names(eff)]

  return(
    list(
      intercept = model$intercept - intercept_adjustment_add,
      eff = eff
    )
  )
}



#' Make GS predictions
#'
#' @param geno bed.matrix returned by `readGenoData()` on wich we want to make
#' prediction
#' @param estim_mark_eff (list, output of `train_gs_model()`) estimated markers
#' effects and intercept.
#'
#' @return data.frame of one column containing the prediction
predict_gs_model <- function(geno, estim_mark_eff) {
  estimated_markers <- row.names(estim_mark_eff$eff)
  geno <- gaston::select.snps(geno, id %in% estimated_markers)
  estim_mark_eff$eff <- estim_mark_eff$eff[geno@snps$id, ]

  additive_gv <- calc_additive_geno(geno, standardized = FALSE) %*% estim_mark_eff$eff$additive

  dominance_gv <- 0
  dominance_model <- (any(!is.na(estim_mark_eff$eff$dominance)) &
    any(estim_mark_eff$eff$dominance != 0, na.rm = TRUE))
  if (dominance_model) {
    dominance_gv <- calc_dominance_geno(geno, standardized = FALSE) %*% estim_mark_eff$eff$dominance
  }

  predicted_gv <- as.data.frame(estim_mark_eff$intercept + additive_gv + dominance_gv)
  return(predicted_gv)
}




#' Repeated K-Folds cross-validation of a GS model
#'
#' @param pheno 1 column data.frame of the phenotype
#' @param geno bed.matrix returned by `readGenoData()`
#' @param with_dominance (default = FALSE) control if the model should include
#' dominance effects
#' @param n_folds (default `10`) number of folds for each repetition
#' @param n_repetitions (default `5`) number of repetition
#'
#' @return list of 2 elements:
#' - `predictions`: data.frame of the predicted values during the
#' cross-validation. Available columns:
#'   - `ind` individual id
#'   - `actual` phenotype value in the training data
#'   - `predicted` predicted values
#'   - `fold` cross-validation fold id
#'   - `repetition` cross-validation repetition id
#' - `metrics`: data.frame of the calculated model metrics during the
#' cross-validation. Available columns are:
#'   - `fold` cross-validation fold id
#'   - `repetition` cross-validation repetition id
#'   - `...` values of the metrics returned by `get_model_metrics()` (one column
#'   per returned list's elements with the same name)
cross_validation_evaluation <- function(pheno,
                                        geno,
                                        with_dominance = TRUE,
                                        n_folds = 10,
                                        n_repetitions = 5) {
  stopifnot(n_folds != 1)

  pheno <- pheno[geno@ped$id, 1, drop = FALSE]

  folds <- caret::createMultiFolds(pheno[, 1], k = n_folds, times = n_repetitions)

  eval_results <- lapply(names(folds), function(fold_id) {
    train_ids <- folds[[fold_id]]

    pheno_train <- pheno[train_ids, 1, drop = FALSE]
    geno_train <- geno[train_ids, ]

    pheno_test <- pheno[-train_ids, 1, drop = FALSE]
    geno_test <- geno[-train_ids, ]

    estim_mark_eff <- train_gs_model(pheno_train, geno_train, with_dominance)
    pheno_pred <- predict_gs_model(geno_test, estim_mark_eff)

    metrics <- get_model_metrics(pheno_test[[1]], pheno_pred[[1]])
    list(
      data = data.frame(
        ind = row.names(pheno_test),
        actual = pheno_test[[1]],
        predicted = pheno_pred[[1]],
        fold = fold_id
      ),
      metrics = metrics
    )
  })

  predictions <- do.call(rbind, lapply(eval_results, function(fold) {
    fold$data
  }))
  predictions$repetition <- sub(".*\\.", "", predictions$fold)

  metrics <- do.call(rbind, lapply(eval_results, function(fold) {
    dta <- as.data.frame(fold$metrics)
    dta$fold <- unique(fold$data$fold)
    dta$repetition <- sub(".*\\.", "", dta$fold)
    dta
  }))

  list(
    predictions = predictions,
    metrics = metrics
  )
}


#' Calculate model metrics
#'
#' @param actual vector of actual values
#' @param predictions vector of model's predicted values
#'
#' @return list of metrics values:
#' - `rmse` root mean square error
#' - `corel_pearson` Pearson's corelation coefficient
#' - `corel_spearman` Spearman's corelation coefficient
#' - `r2` R squared
get_model_metrics <- function(actual, predictions) {
  error <- actual - predictions
  list(
    rmse = sqrt(mean(error^2)),
    corel_pearson = cor(actual, predictions, method = "pearson"),
    corel_spearman = cor(actual, predictions, method = "spearman"),
    r2 = 1 - sum(error^2) / sum((actual - mean(actual))^2)
  )
}



#' Check if a dominance model is applicable with the provided genetic data
#'
#' @param dominance list returned by
#' @param homozygous_threshold A threshold used to identify individuals or SNPs
#' with a homozygosity proportion exceeding this value. These will be counted
#' and reported in the error message if applicable.
#' This parameter does not influence the function's behavior, only the error
#' message it can raise.
#'
#' @details
#' The dominance model need the dominance relationship matrix to be invertible
#' in order to be able to calculate the dominance effects. This function will
#' return an error if it is not the case.
#' The error will contain information about the homozygousity of the markers and
#' individuals as if the data have too many homozygous individuals/markers it
#' will probably not be suited for a dominance model. The ``
#'
#'
#' @return `TRUE` or raise an `engineError`
check_dominance_model_is_applicable <- function(dominance, homozygous_threshold = 0.95) {
  if (det(dominance$rel_mat) == 0) {
    # dominance$rel_mat is not invertible
    # Sometime for computing precision reason, `det(dominance$rel_mat)` can be
    # 0 but the matrix can yet be invertible (`solve` function is more
    # robust than `det` function)
    inversible <- tryCatch(
      {
        solve(dominance$rel_mat)
        TRUE
      },
      error = function(err) {
        return(FALSE)
      }
    )

    if (inversible) {
      return(TRUE)
    }

    snp_homozygous_proportion <- 1 - colSums(dominance$geno_mat) / nrow(dominance$geno_mat)
    homozygous_snp <- snp_homozygous_proportion > homozygous_threshold
    n_homozygous_snp <- sum(homozygous_snp)
    homozygous_snp_proportion <- n_homozygous_snp / ncol(dominance$geno_mat)

    ind_homozygous_proportion <- 1 - rowSums(dominance$geno_mat) / ncol(dominance$geno_mat)
    homozygous_ind <- ind_homozygous_proportion > homozygous_threshold
    n_homozygous_ind <- sum(homozygous_ind)

    engineError(
      message = paste(
        "Dominace model is not applicable with the provided data as",
        "the dominance relationship matrix is not invertible.",
        "FYI,", n_homozygous_ind, "individuals have a",
        "homozygousity higher than", homozygous_threshold * 100,
        "%, and", signif(homozygous_snp_proportion * 100, 2), "% of",
        "the markers have a homozygousity proportion higher than",
        homozygous_threshold * 100, "%."
      ),
      extra = list(
        code = errorCode("DOMINANCE_MODEL_NOT_APPLICABLE"),
        n_homozygous_ind = n_homozygous_ind,
        n_homozygous_snp = n_homozygous_snp,
        n_ind = nrow(dominance$geno_mat),
        n_snp = ncol(dominance$geno_mat),
        homozygous_threshold = homozygous_threshold
      )
    )
  }
  return(TRUE)
}
