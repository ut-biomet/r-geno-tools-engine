# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Test engine commands






capture_output({
  Sys.setenv(RGENOROOT = normalizePath('../..'))

  rGenoCommand <- paste0("Rscript ",
                         Sys.getenv('RGENOROOT'),
                         "/r-geno-tools-engine.R")

  # help ----
  test_that('engine help', {
    expect_error({
      x <- system(rGenoCommand,
                  intern = FALSE,
                  ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  # all commands help ----
  x <- system(rGenoCommand,
              intern = TRUE)
  x <- x[which(grepl('Available commands', x)) + 1]
  x <- gsub('[^\\w,\\-]', '', x, perl = TRUE)
  cmds <- c('', strsplit(x, split = ',')[[1]])
  cmds <- paste(cmds, '--help')
  for (cmd in cmds) {
    test_that(paste('help:', cmd), {
      expect_error({
        x <- system(paste(rGenoCommand, cmd),
                    intern = FALSE,
                    ignore.stdout = TRUE)
      }, NA)
      expect_equal(x, 0)
    })
  }

  # gwas ----
  test_that('gwas', {
    cmd <- paste(rGenoCommand, 'gwas',
                 '--genoFile "$RGENOROOT/data/geno/testMarkerData01.vcf.gz"',
                 '--phenoFile "$RGENOROOT/data/pheno/testPhenoData01.csv"',
                 '--trait "Flowering.time.at.Arkansas"',
                 '--test "score"',
                 '--fixed 0',
                 '--response "quantitative"',
                 '--thresh-maf 0.05',
                 '--thresh-callrate 0.95',
                 '--outFile "$RGENOROOT/tests/testthat/testOutput/gwasRes.json"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })


  # gwas-manplot ----
  test_that('gwas-manplot', {

    cmd <- paste(rGenoCommand, 'gwas-manplot',
                 '--gwasFile "$RGENOROOT/tests/testthat/testOutput/gwasRes.json"',
                 '--adj-method "bonferroni"',
                 '--thresh-p 0.05',
                 '--filter-nPoints 3000',
                 '--outFile "$RGENOROOT/tests/testthat/testOutput/manPlot.html"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })




  # gwas-adjResults ----
  test_that('gwas-adjResults', {

    cmd <- paste(rGenoCommand, 'gwas-adjresults',
                 '--gwasFile "$RGENOROOT/tests/testthat/testOutput/gwasRes.json"',
                 '--adj-method "bonferroni"',
                 '--filter-nPoints 3000',
                 '--outFile "$RGENOROOT/tests/testthat/testOutput/gwasRes_adj.json"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })




  # ldplot ----
  test_that('ldplot', {

    cmd <- paste(rGenoCommand, 'ldplot',
                 '--genoFile "$RGENOROOT/data/geno/testMarkerData01.vcf.gz"',
                 '--from 42',
                 '--to 62',
                 '--outFile "$RGENOROOT/tests/testthat/testOutput/ldplot.png"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })



  # relmat-ped ----
  test_that('relmat-ped', {

    cmd <- paste(rGenoCommand, 'relmat-ped',
                 '--pedFile "$RGENOROOT/data/pedigree/testPedData_char.csv"',
                 '--outFile "$RGENOROOT/tests/testthat/testOutput/pedRelMat.json"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  # relmat-geno ----
  test_that('relmat-geno', {

    cmd <- paste(rGenoCommand, 'relmat-geno',
                 '--genoFile "$RGENOROOT/data/geno/breedGame_geno.vcf.gz"',
                 '--outFile "$RGENOROOT/tests/testthat/testOutput/genoRelMat.json"')
    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  # relmat-combined ----
  test_that('relmat-combined', {

    cmd <- paste(rGenoCommand, 'relmat-combined',
                 '--ped-relmatFile $RGENOROOT/data/results/breedGame_pedRelMat.csv',
                 '--geno-relmatFile $RGENOROOT/data/results/breedGame_genoRelMat.csv',
                 # '--combine-method Legarra',
                 '--combine-method Martini',
                 '--tau 1.3',
                 '--omega 0.7',
                 '--outFile "$RGENOROOT/tests/testthat/testOutput/genoRelMat.json"')
    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })




  # relmat-heatmap ----
  test_that('relmat-heatmap', {

    cmd <- paste(rGenoCommand, 'relmat-heatmap',
                 '--relmatFile "$RGENOROOT/tests/testthat/testOutput/pedRelMat.json"',
                 '--outFile "$RGENOROOT/tests/testthat/testOutput/relMat.png"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })



  # pedNetwork ----
  test_that('pedNetwork', {

    cmd <- paste(rGenoCommand, 'pedNetwork',
                 '--pedFile "$RGENOROOT/data/pedigree/testPedData_char.csv"',
                 '--outFile "$RGENOROOT/tests/testthat/testOutput/pedNet.html"')


    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })



  # crossing-simulation ----
  test_that('crossing-simulation', {

    cmd <- paste(
      rGenoCommand, 'crossing-simulation',
      '--genoFile "$RGENOROOT/data/geno/breedGame_phasedGeno.vcf.gz"',
      '--crossTableFile "$RGENOROOT/data/crossingTable/breedGame_crossTable.csv"',
      '--SNPcoordFile "$RGENOROOT/data/SNPcoordinates/breedingGame_SNPcoord.csv"',
      '--outFile "$RGENOROOT/tests/testthat/testOutput/crossSim.vcf.gz"'
    )

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  # progeny-blup-calculation ----
  test_that('progeny-blup-calculation', {
    cmd <- paste(
      rGenoCommand, 'progeny-blup-calculation',
      '--genoFile "$RGENOROOT/data/geno/breedGame_phasedGeno.vcf.gz"',
      '--crossTableFile "$RGENOROOT/data/crossingTable/breedGame_small_crossTable.csv"',
      '--SNPcoordFile "$RGENOROOT/data/SNPcoordinates/breedingGame_SNPcoord.csv"',
      '--markerEffectsFile "$RGENOROOT/data/markerEffects/breedGame_markerEffects.csv"',
      '--outFile "$RGENOROOT/tests/testthat/testOutput/progBlups.json"'
    )

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  test_that('progeny-blup-calculation (multi-trait marker effect)', {
    cmd <- paste(
      rGenoCommand, 'progeny-blup-calculation',
      '--genoFile "$RGENOROOT/data/geno/breedGame_phasedGeno.vcf.gz"',
      '--crossTableFile "$RGENOROOT/data/crossingTable/breedGame_small_crossTable.csv"',
      '--SNPcoordFile "$RGENOROOT/data/SNPcoordinates/breedingGame_SNPcoord.csv"',
      '--markerEffectsFile "$RGENOROOT/data/markerEffects/breedGame_markerEffects_2traits.json"',
      '--outFile "$RGENOROOT/tests/testthat/testOutput/progBlups_2traits.json"'
    )

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })


  # progeny-blup-plot ----
  test_that('progeny-blup-plot', {
    cmd <- paste(
      rGenoCommand, 'progeny-blup-plot',
      '--progeniesBlupFile "$RGENOROOT/tests/testthat/testOutput/progBlups.json"',
      '--y-axis-name "Phenotypic trait"',
      '--error-bar-interval 0.95',
      '--outFile "$RGENOROOT/tests/testthat/testOutput/blupPlot.html"'
    )

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  test_that('progeny-blup-plot (multi-trait file)', {
    cmd <- paste(
      rGenoCommand, 'progeny-blup-plot',
      '--progeniesBlupFile "$RGENOROOT/tests/testthat/testOutput/progBlups_2traits.json"',
      '--y-axis-name "Phenotypic trait"',
      '--error-bar-interval 0.95',
      '--trait "trait1"',
      '--outFile "$RGENOROOT/tests/testthat/testOutput/blupPlot.html"'
    )

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  test_that('progeny-blup-plot-2-traits', {
    cmd <- paste(
      rGenoCommand, 'progeny-blup-plot-2-traits',
      '--progeniesBlupFile "$RGENOROOT/tests/testthat/testOutput/progBlups_2traits.json"',
      '--x-trait "trait1"',
      '--y-trait "trait2"',
      '--confidence-level 0.95',
      '--outFile "$RGENOROOT/tests/testthat/testOutput/blupPlot_2traits.html"'
    )

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })





})
