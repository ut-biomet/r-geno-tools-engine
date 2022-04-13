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
      # expect_true(isSymmetric(relMat)) # exact symmetry fail
      # check approximate symmetry (diff< ~1.5e-8.) instead
      expect_true(all.equal.numeric(relMat, t(relMat)))

      expect_true(!is.null(colnames(relMat)))
      expect_true(!is.null(row.names(relMat)))
      expect_identical(colnames(relMat), row.names(relMat))
      expect_identical(colnames(params$ped_rm), colnames(relMat))
    })
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



})
