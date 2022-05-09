# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# relationship matrix tests.

capture_output({

  # pedRelMat -----
  pedFiles <- c('../../data/pedigree/testPedData_char.csv',
                '../../data/pedigree/testPedData_num.csv',
                '../../data/pedigree/testPedData_missFounder.csv')
  for (file in pedFiles) {
    test_that(paste("pedRelMat", basename(file)), {
      suppressWarnings({
        ped <- readPedData(file)
      })
      expect_error({
        relMat <- pedRelMat(ped)
      }, NA)
      expect_true(is.matrix(relMat))
      expect_true(is.numeric(relMat))
      expect_true(isSymmetric(relMat))
      expect_true(!any(is.na(relMat)))
      expect_true(all(relMat >= 0))
      expect_true(all(relMat <= 2))
      expect_true(!is.null(colnames(relMat)))
      expect_true(!is.null(row.names(relMat)))
      expect_identical(colnames(relMat), row.names(relMat))
      expect_identical(ped$data$ind, colnames(relMat))
    })
    test_that(paste("pedRelMat vs AGHmatrix", basename(file)),{
      skip_if_not_installed("AGHmatrix")
      suppressWarnings({
        ped <- readPedData(file)
      })
      exp_relMat <- AGHmatrix::Amatrix(replace(ped$data, is.na(ped$data), 0))
      relMat <- pedRelMat(ped)
      expect_identical(colnames(relMat), colnames(exp_relMat))
      expect_identical(row.names(relMat), row.names(exp_relMat))
      expect_identical(relMat, exp_relMat)
    })
  }


  # genoRelMat ----
  genoFiles <- c('../../data/geno/breedGame_geno.vcf.gz')
  for (file in genoFiles) {
    test_that(paste("genoRelMat", basename(file)), {
      geno <- readGenoData(file)
      # reduce `geno` size to be faster
      geno <- gaston::select.inds(geno, id %in% sample(geno@ped$id, 100))

      expect_error({
        relMat <- genoRelMat(geno)
      }, NA)
      expect_true(is.matrix(relMat))
      expect_true(is.numeric(relMat))
      expect_true(isSymmetric(relMat))
      expect_true(!any(is.na(relMat)))
      expect_true(!is.null(colnames(relMat)))
      expect_true(!is.null(row.names(relMat)))
      expect_identical(colnames(relMat), row.names(relMat))
      expect_identical(geno@ped$id, colnames(relMat))
    })
  }


  # combinedRelMat ----
  A <- pedRelMat(readPedData('../../data/pedigree/breedGame_pedigree.csv'))
  sampledInds <- sample(row.names(A), 100)
  A <- A[sampledInds, sampledInds]
  geno <- readGenoData('../../data/geno/breedGame_geno.vcf.gz')
  geno <- gaston::select.inds(geno, id %in% sampledInds)
  G <- genoRelMat(geno)

  paramList <- list(
    default = list(ped_rm = A, geno_rm = G),
    legarra = list(ped_rm = A, geno_rm = G, method = 'Legarra'),
    martini = list(ped_rm = A, geno_rm = G, method = 'Martini',
                   tau = 0.5, omega = 0.5),
    martini_randParams = list(ped_rm = A, geno_rm = G, method = 'Martini',
                   tau = runif(1, 0, 20), omega = runif(1, -20, 1)),
    warn_missIndInPed = list(ped_rm = A[-(1:10), -(1:10)], geno_rm = G),
    warn_tauOmeg_leggara = list(ped_rm = A, geno_rm = G, method = 'Legarra',
                                tau = 0.5, omega = 0.5)
  )

  for (paramName in names(paramList)) {
    test_that(paste('combinedRelMat:', paramName), {

      params <- paramList[[paramName]]
      if (grepl('warn_', paramName)) {
        expect_warning({
          relMat <- do.call(combinedRelMat, params)
        })
      } else {
        expect_error({
          relMat <- do.call(combinedRelMat, params)
        }, NA)
      }
      expect_true(is.matrix(relMat))
      expect_true(is.numeric(relMat))
      expect_identical(dim(relMat), dim(params$ped_rm))
      expect_true(!any(is.na(relMat)))
      expect_true(isSymmetric(relMat)) # exact symmetry
      # check approximate symmetry (diff< ~1.5e-8.) instead
      # expect_true(all.equal.numeric(relMat, t(relMat)))

      expect_true(!is.null(colnames(relMat)))
      expect_true(!is.null(row.names(relMat)))
      expect_identical(colnames(relMat), row.names(relMat))
      expect_identical(colnames(params$ped_rm), colnames(relMat))
    })
    if (!grepl('warn_', paramName)) {
      test_that(paste("combineRelMat vs AGHmatrix:", paramName),{
        skip_if_not_installed("AGHmatrix")
        params <- paramList[[paramName]]
        suppressWarnings({
          relMat <- do.call(combinedRelMat, params)
        })
        agh_params <- list()
        agh_params$A <- params$ped_rm
        agh_params$G <- params$geno_rm
        agh_params$method <- "Martini"
        agh_params$omega <- params$omega
        agh_params$tau <- params$tau
        exp_relMat <- do.call(AGHmatrix::Hmatrix, agh_params)
        exp_relMat <- exp_relMat[row.names(agh_params$A),
                                 colnames(agh_params$A)]

        expect_identical(colnames(relMat), colnames(exp_relMat))
        expect_identical(row.names(relMat), row.names(exp_relMat))
        expect_true(all.equal.numeric(exp_relMat, relMat))
      })
    }
  }

  test_that(paste('combinedRelMat: wrong method'), {
    expect_error({
      relMat <- combinedRelMat(ped_rm = A, geno_rm = G, method = 'doNotExist')
    }, 'Method should be either "Legarra" or "Marini".')
    expect_error({
      relMat <- combinedRelMat(ped_rm = A, geno_rm = G,
                               method = c('Legarra', 'Martini'))
    }, 'Only one method should be provided to `combineRelMat` function.')
  })

  test_that(paste('combinedRelMat: no common inds'), {
    A2 <- A
    colnames(A2) <- row.names(A2) <- paste0('toto', seq(nrow(A)))

    expect_error({
      relMat <- combinedRelMat(ped_rm = A2, geno_rm = G)
    }, paste('No common individuals between genomic',
             'and pedigree relationship matrices.'))
  })


  test_that(paste('combinedRelMat: wron tau omega'), {
    expect_error({
      relMat <- combinedRelMat(ped_rm = A,
                               geno_rm = G,
                               method = 'Martini',
                               tau = 0,
                               omega = 1)
    }, 'The combination `tau`= 0, and `omega` = 1 is not possible.')
    expect_error({
      relMat <- combinedRelMat(ped_rm = A,
                               geno_rm = G,
                               method = 'Martini',
                               tau = -1,
                               omega = 1)
    }, '`tau` must be a positive number.')
    expect_error({
      relMat <- combinedRelMat(ped_rm = A,
                               geno_rm = G,
                               method = 'Martini',
                               tau = 0,
                               omega = 1.5)
    }, '`omega` must be lower than one.')
  })

})
