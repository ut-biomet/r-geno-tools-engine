# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Relationship matrices calculation functions

#' Pedigree Relationship Matrix calculation
#'
#' @param ped List return by `readPedData` function
#'
#' @return matrix
#' @author Hiroyoshi Iwata, Julien Diot
pedRelMat <- function(ped) {
  logger <- Logger$new("r-pedRelMat()")

  ### Check input ----
  logger$log('Check inputs ...')
  if (!is.list(ped)) {
    bad_argument(ped, must_be = "a list", not = ped, "type")
  }
  expected_names <- c("data", "graph")
  if (!identical(names(ped), expected_names)) {
    bad_argument("names(ped)", must_be = expected_names, not = names(ped))
  }
  if (!is.data.frame(ped$data)) {
    bad_argument("ped$data", must_be = "a data.frame", not = ped$data, "type")
  }
  expected_columns <- c('ind', 'parent1', 'parent2')
  if (!identical(colnames(ped$data), expected_columns)) {
    bad_argument("colnames(ped$data)", must_be = expected_columns, not = colnames(ped$data))
  }
  logger$log('Check inputs DONE')


  ### Create the look-up table ----
  ### this table is the same than the pedigree table with id given to the
  ### individuals according to their generation.
  logger$log('Create look-up table ...')

  lut <- as.data.frame(matrix(rep(NA, 3 * nrow(ped$data)), ncol = 3))
  colnames(lut) <- c('id', 'parent1', 'parent2')

  # initialize id of unknown parents
  lut$parent1[is.na(ped$data$parent1)] <- 0
  lut$parent2[is.na(ped$data$parent2)] <- 0
  gen <- rep(0, nrow(ped$data)) # vector of generation

  i <- 0
  while (sum(is.na(lut)) != 0 && i <= (nrow(ped$data) + 1)) {
    i <- i + 1 # to not be stuck in the loop

    # select ind with both known parents
    selector <- is.na(lut$id) & !is.na(lut$parent1) & !is.na(lut$parent2)
    if (!any(selector)) {
      # should never come here (pedigree checked in `readPedData`)
      stop("Can't create the look-up table, no Individuals with known parents")
    }

    # give id to those individuals
    lut$id[selector] <- max(c(lut$id, 0), na.rm = T) + seq(sum(selector))
    # update generation count
    gen[selector] <- max(gen) + 1

    # fill the look-up table with these new ids
    selectedId <- lut$id[selector]
    selectedNames <- ped$data$ind[selector]
    for (j in seq(sum(selector))) {
      lut$parent1[ped$data$parent1 == selectedNames[j]] <- selectedId[j]
      lut$parent2[ped$data$parent2 == selectedNames[j]] <- selectedId[j]
    }
  }
  lut$name <- ped$data$ind
  lut <- lut[order(lut$id),] # re-order the table
  logger$log('Create look-up table DONE')



  ### Calculate relationship matrix
  logger$log('Calculate relationship matrix ...')
  relMat <- matrix(NA, nrow = nrow(ped$data), ncol = nrow(ped$data))
  for (i in seq(nrow(lut))) {
    # l <- which(lut$id == i) # l is always = to i because lut is ordered.
    # p1 <- lut[l, 'parent1']
    # p2 <- lut[l, 'parent2']

    p1 <- lut$parent1[i]
    p2 <- lut$parent2[i]

    for (j in seq(i)) {
      if (j != i) {
        if (p1 != 0 & p2 != 0) {
          relMat[i, j] <- relMat[j, i] <- 0.5 * (relMat[j, p1] + relMat[j, p2])
        } else if (p1 == 0 & p2 == 0) {
          relMat[i, j] <- relMat[j, i] <- 0
        } else if (p1 == 0) {
          relMat[i, j] <- relMat[j, i] <- 0.5 * relMat[j, p2]
        } else if (p2 == 0) {
          relMat[i, j] <- relMat[j, i] <- 0.5 * relMat[j, p1]
        }
      } else {
        if (p1 != 0 & p2 != 0) {
          relMat[i, i] <- 1 + 0.5 * relMat[p1, p2]
        } else {
          relMat[i, i] <- 1
        }
      }
    }
  }
  row.names(relMat) <- colnames(relMat) <- lut$name
  # because lut is ordered, the above line is similar to:
  #   names <- sapply(seq(nrow(relMat)), function(i){lut[lut$id == i, 'name']})
  #   row.names(relMat) <- colnames(relMat) <- name


  # reset initial order
  relMat <- relMat[ped$data$ind, ped$data$ind]

  logger$log('Calculate relationship matrix DONE')

  ### output ----
  logger$log("DONE, return output.")
  return(relMat)
}






#' Additive Genomic Matrix
#'
#' @param geno `gaston::bed.matrix` return by `readGenoData` function
#' @param standardized `boolean` (default TRUE) control if the returned genetic
#' matrix should be standardized.
#'
#' @details The standardization is made with: (x - 2*p) / sqrt(2*p*(1 - p)) with x the
#' genetic values in alleles dose and p the allelic frequency.
#'
#' @return matrix
calc_additive_geno <- function(geno, standardized = TRUE){
  if (standardized) {
    gaston::standardize(geno) <- "p"
    # this is equivalent of `(x - 2*p) / sqrt(2*p*(1 - p))`
    # with x in allele dose
  } else {
    gaston::standardize(geno) <- "none"
  }
  return(gaston::as.matrix(geno))
}

#' Additive Genomic Relationship Matrix calculation
#'
#' @param geno `gaston::bed.matrix` return by `readGenoData` function
#' @param standardized `boolean` (default TRUE) control if the calculation
#' should be done on the standardized additive genetic matrix.
#'
#' @return list of 2 elements:
#' - `geno_mat`: the genetic matrix used for calculating the relationship matrix
#' - `rel_mat`: the relationship matrix
calc_additive_rel_mat <- function(geno, standardized = TRUE) {
  if (standardized) {
    return(
      list(
        geno_mat = calc_additive_geno(geno, standardized = TRUE),
        rel_mat = gaston::GRM(geno, autosome.only = FALSE)
      )
    )
  }

  X <- calc_additive_geno(geno, standardized = FALSE)
  list(
    geno_mat = X,
    rel_mat = tcrossprod(X) / ncol(X)
  )
}

#' Dominance Genomic Matrix
#'
#' @param geno `gaston::bed.matrix` return by `readGenoData` function
#' @param standardized `boolean` (default TRUE) control if the returned genetic
#' matrix should be standardized.
#'
#' @details The standardization is made with the matrix with entries
#' `p/(1 - p)`, `-1`, `(1-p)/p` according to the values `0`, `1`, `2` in the
#' genetic matrix, with p the allele frequency (cf. `?gaston::DM`)
#'
#' @return matrix
calc_dominance_geno <- function(geno, standardized = TRUE){
  gaston::standardize(geno) <- "none"

  if (standardized) {
    D <- apply(gaston::as.matrix(geno), MARGIN = 2, function(x) {
      af <- sum(x) / (2*length(x))
      i <- af/(1 - af)
      j <- -1
      k <- (1 - af)/af

      (
        i * as.numeric(x == 0)
        + j * as.numeric(x == 1)
        + k * as.numeric(x == 2)
      )
    })
    colnames(D) <- geno@snps$id
    row.names(D) <- geno@ped$id
    return(D)
  }

  return(-abs(gaston::as.matrix(geno) - 1) + 1)
}

#' Dominance Genomic Relationship Matrix calculation
#'
#' @param geno `gaston::bed.matrix` return by `readGenoData` function
#' @param standardized `boolean` (default TRUE) control if the calculation
#' should be done on the standardized dominance genetic matrix.
#'
#' @return list of 2 elements:
#' - `geno_mat`: the genetic matrix used for calculating the relationship matrix
#' - `rel_mat`: the relationship matrix
calc_dominance_rel_mat <- function(geno, standardized = TRUE) {
  if (standardized) {
    return(
      list(
        geno_mat = calc_dominance_geno(geno, standardized = TRUE),
        rel_mat = gaston::DM(geno, autosome.only = FALSE)
      )
    )
  }

  D <- calc_dominance_geno(geno, standardized = FALSE)
  return(
    list(
      geno_mat = D,
      rel_mat = tcrossprod(D) / ncol(D)
    )
  )
}


#' Combined (pedigree + genomic) Relationship Matrix
#'
#' Correct a pedigree relationship matrix by using genomic relationship matrix.
#'
#' @param ped_rm pedigree relationship matrix (matrix from function pedRelMat)
#' @param geno_rm genomic relationship matrix (matrix from function genoRelMat)
#' @param method method to use, either "Legarra" or "Martini"
#' @param tau tau parameter of the Martini's method
#' @param omega omega parameter of the Martini's method
#'
#' @return matrix
#' @author Hiroyoshi Iwata, Julien Diot
combinedRelMat <- function(ped_rm,
                           geno_rm,
                           method = 'Legarra',
                           tau = NULL,
                           omega = NULL) {

  logger <- Logger$new("r-combineRelMat()")

  ### Check input ----
  logger$log('Check inputs ...')
  checkRelMat(ped_rm)
  checkRelMat(geno_rm)

  if (method == 'Legarra') {
    tau <- 1
    omega <- 1
  }

  additionalInds <- rownames(geno_rm)[!rownames(geno_rm) %in% rownames(ped_rm)]
  if (length(additionalInds) != 0) {
    warning(paste(
      length(additionalInds),
      'individuals of the genomic relationship matrix',
      'are not in the pedigree relationship matrix:',
      paste(additionalInds, collapse = ', '),
      '\nThe combined relationship matrix will only include individuals of the',
      'pedigree relationship matrix.'
    ))
    geno_rm <- geno_rm[row.names(geno_rm) %in% colnames(ped_rm),
                       colnames(geno_rm) %in% colnames(ped_rm)]
  }
  logger$log('Check inputs DONE')

  # calculate genetic relationship matrix
  logger$log('Calculate combined relationship matrix ...')

  # the following is an alternative to `AGHmatrix` by Iwata-sensei
  A.name <- rownames(ped_rm)
  G.name <- rownames(geno_rm)

  is.inG <- A.name %in% G.name
  A1.name <- A.name[!is.inG]
  A2.name <- A.name[is.inG]

  A00 <- ped_rm[c(A1.name, A2.name), c(A1.name, A2.name)]
  G22 <- geno_rm[A2.name, A2.name]

  A12 <- A00[A1.name, A2.name]
  A22 <- A00[A2.name, A2.name]

  G22.inv <- solve(G22)
  A22.inv <- solve(A22)
  H22 <- solve(tau * G22.inv + (1 - omega) * A22.inv)

  T11 <- A12 %*% A22.inv %*% (H22 - A22) %*% A22.inv %*% t(A12)
  T12 <- A12 %*% A22.inv %*% (H22 - A22)
  hrm <- A00 + rbind(cbind(T11, T12), cbind(t(T12), H22 - A22))


  rownames(hrm) <- colnames(hrm) <- rownames(A00)
  hrm <- hrm[row.names(ped_rm), colnames(ped_rm)]

  # `hrm` is not exactly symetric due to computer precision
  hrm <- (hrm + t(hrm)) / 2 # make hrm symetric by taking the mean value

  logger$log('Calculate combined relationship matrix DONE')

  ### output ----
  logger$log("DONE, return output.")
  return(hrm)
}
