capture.output({
  # GS model ----
  tests_cases <- list(
    test_01 = list(
      geno = "../../data/genomic_selection/geno_G1.vcf.gz",
      pheno = "../../data/genomic_selection/pheno_train.csv",
      trait = "pheno"
    ),
    missing_genotype = list(
      geno = "../../data/genomic_selection/geno_G1.vcf.gz",
      pheno = "../data/GS_pheno_train_missing_genotype.csv",
      trait = "pheno"
    ),
    missing_phenotype = list(
      geno = "../../data/genomic_selection/geno_G1.vcf.gz",
      pheno = "../data/GS_pheno_train_missing_phenotype.csv",
      trait = "pheno"
    ),
    duplicated_snp_ids = list(
      geno = "../data/geno_G1_duplicated_ids.vcf.gz",
      pheno = "../../data/genomic_selection/pheno_train.csv",
      trait = "pheno"
    ),
    unused_trait_badly_formatted = list(
      geno = "../data/geno_G1_duplicated_ids.vcf.gz",
      pheno = "../../data/genomic_selection/pheno_train_with_pheno2_trait_badly_formated.csv",
      trait = "pheno"
    )
  )

  for (test in names(tests_cases)) {
    for (with_dominance in c(T, F)) {
      test_that(paste0("train_gs_model_main ", test, " dominance: ", with_dominance), {
        expect_no_error({
          model <- train_gs_model_main(
            genoFile = tests_cases[[test]]$geno,
            phenoFile = tests_cases[[test]]$pheno,
            genoUrl = NULL,
            phenoUrl = NULL,
            trait = tests_cases[[test]]$trait,
            with_dominance = with_dominance,
            thresh_maf = 0,
            outFile = tempfile(fileext = ".json")
          )
        })
        expect_no_error({
          model <- train_gs_model_main(
            genoFile = tests_cases[[test]]$geno,
            phenoFile = tests_cases[[test]]$pheno,
            genoUrl = NULL,
            phenoUrl = NULL,
            trait = tests_cases[[test]]$trait,
            with_dominance = with_dominance,
            thresh_maf = 0,
            n_markers = 2000,
            outFile = tempfile(fileext = ".json")
          )
        })
        expect_equal(1, 1) # testthat v3.2.2 have a bug with `expect_no_error`
        # if it succeed the test is considered as "skipped".
      })
    }
  }


  tests_bad_cases <- list(
    no_heterozygot = list(
      geno = "../../data/geno/testMarkerData01.vcf.gz",
      pheno = "../../data/pheno/testPhenoData01.csv",
      trait = "Flowering.time.at.Arkansas"
    )
  )
  for (test in names(tests_bad_cases)) {
    test_that(paste0("train_gs_model_main Bad case ", test), {
      expect_engineError({
        model <- train_gs_model_main(
          genoFile = tests_bad_cases[[test]]$geno,
          phenoFile = tests_bad_cases[[test]]$pheno,
          genoUrl = NULL,
          phenoUrl = NULL,
          trait = tests_bad_cases[[test]]$trait,
          with_dominance = TRUE,
          thresh_maf = 0,
          outFile = tempfile(fileext = ".json")
        )
      })
    })
  }


  # GS cross-validation ----
  tests_cases <- list(
    test_01 = list(
      geno = "../../data/genomic_selection/geno_G1.vcf.gz",
      pheno = "../../data/genomic_selection/pheno_train.csv",
      trait = "pheno"
    ),
    missing_genotype = list(
      geno = "../../data/genomic_selection/geno_G1.vcf.gz",
      pheno = "../data/GS_pheno_train_missing_genotype.csv",
      trait = "pheno"
    ),
    missing_phenotype = list(
      geno = "../../data/genomic_selection/geno_G1.vcf.gz",
      pheno = "../data/GS_pheno_train_missing_phenotype.csv",
      trait = "pheno"
    )
  )

  for (test in names(tests_cases)) {
    test_that(paste0("cross_validation_evaluation_main ", test), {
      for (with_dominance in c(T, F)) {
        expect_no_error({
          evaluation <- cross_validation_evaluation_main(
            genoFile = tests_cases[[test]]$geno,
            phenoFile = tests_cases[[test]]$pheno,
            genoUrl = NULL,
            phenoUrl = NULL,
            trait = tests_cases[[test]]$trait,
            with_dominance = with_dominance,
            n_folds = 10,
            n_repetitions = 5,
            thresh_maf = 0,
            outFile = tempfile(fileext = ".json")
          )
        })
        expect_no_error({
          evaluation <- cross_validation_evaluation_main(
            genoFile = tests_cases[[test]]$geno,
            phenoFile = tests_cases[[test]]$pheno,
            genoUrl = NULL,
            phenoUrl = NULL,
            trait = tests_cases[[test]]$trait,
            with_dominance = with_dominance,
            n_folds = 10,
            n_repetitions = 5,
            thresh_maf = 0,
            n_markers = 2000,
            outFile = tempfile(fileext = ".json")
          )
        })
      }
    })
  }


  # GS evaluation plot ----
  tests_cases <- list(
    test_01 = list(
      geno = "../../data/genomic_selection/geno_G1.vcf.gz",
      pheno = "../../data/genomic_selection/pheno_train.csv",
      trait = "pheno"
    ),
    missing_genotype = list(
      geno = "../../data/genomic_selection/geno_G1.vcf.gz",
      pheno = "../data/GS_pheno_train_missing_genotype.csv",
      trait = "pheno"
    ),
    missing_phenotype = list(
      geno = "../../data/genomic_selection/geno_G1.vcf.gz",
      pheno = "../data/GS_pheno_train_missing_phenotype.csv",
      trait = "pheno"
    )
  )

  for (test in names(tests_cases)) {
    test_that(paste0("cross_validation_evaluation_main ", test), {
      for (with_dominance in c(T, F)) {
        evaluation <- cross_validation_evaluation_main(
          genoFile = tests_cases[[test]]$geno,
          phenoFile = tests_cases[[test]]$pheno,
          trait = tests_cases[[test]]$trait,
          with_dominance = with_dominance,
          n_folds = 10,
          n_repetitions = 5,
          thresh_maf = 0,
          outFile = tempfile(fileext = ".json")
        )
        expect_no_error({
          draw_evaluation_plot(
            evaluation$file,
            tempfile(fileext = ".html")
          )
        })
      }
    })
  }

  # GS predictions ----

  tests_cases <- list(
    additive = list(
      geno = "../../data/genomic_selection/geno_G2.vcf.gz",
      markerEffects = "../../data/results/GS_model_additive.json"
    ),
    dominance = list(
      geno = "../../data/genomic_selection/geno_G2.vcf.gz",
      markerEffects = "../../data/results/GS_model_dominance.json"
    ),
    several_traits = list(
      geno = "../../data/geno/breedGame_phasedGeno.vcf.gz",
      markerEffects = "../../data/markerEffects/breedGame_markerEffects_2traits.json"
    )
  )

  for (test in names(tests_cases)) {
    test_that(paste0("train_gs_model_main ", test), {
      expect_no_error({
        predictions <- predict_gs_model_main(
          genoFile = tests_cases[[test]]$geno,
          genoUrl = NULL,
          markerEffectsFile = tests_cases[[test]]$markerEffects,
          markerEffectsUrl = NULL,
          outFile = tempfile(fileext = ".csv")
        )
      })
    })
  }
})
