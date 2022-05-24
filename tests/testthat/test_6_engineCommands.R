# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Test engine commands






capture_output({
  setwd('../..')

  # help ----
  test_that('engine help', {
    expect_error({
      x <- system("Rscript ./r-geno-tools-engine.R",
                  intern = FALSE,
                  ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  # all commands help ----
  x <- system("Rscript ./r-geno-tools-engine.R",
              intern = TRUE)
  x <- x[which(grepl('Available commands', x)) + 1]
  x <- gsub('[^\\w,\\-]', '', x, perl = TRUE)
  cmds <- c('', strsplit(x, split = ',')[[1]])
  cmds <- paste(cmds, '--help')
  for (cmd in cmds) {
    test_that(paste('help:', cmd), {
      expect_error({
        x <- system(paste("Rscript ./r-geno-tools-engine.R", cmd),
                    intern = FALSE,
                    ignore.stdout = TRUE)
      }, NA)
      expect_equal(x, 0)
    })
  }

  # gwas ----
  test_that('gwas', {
    cmd <- paste('Rscript ./r-geno-tools-engine.R gwas',
                 '--genoFile "data/geno/testMarkerData01.vcf"',
                 '--phenoFile "data/pheno/testPhenoData01.csv"',
                 '--trait "Flowering.time.at.Arkansas"',
                 '--test "score"',
                 '--fixed 0',
                 '--response "quantitative"',
                 '--thresh-maf 0.05',
                 '--thresh-callrate 0.95',
                 '--outFile "tests/testthat/testOutput/gwasRes.json"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })


  # gwas-manplot ----
  test_that('gwas-manplot', {

    cmd <- paste('Rscript ./r-geno-tools-engine.R gwas-manplot',
                 '--gwasFile "tests/testthat/testOutput/gwasRes.json"',
                 '--adj-method "bonferroni"',
                 '--thresh-p 0.05',
                 '--filter-nPoints 3000',
                 '--outFile "tests/testthat/testOutput/manPlot.html"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })




  # gwas-adjResults ----
  test_that('gwas-adjResults', {

    cmd <- paste('Rscript ./r-geno-tools-engine.R gwas-adjresults',
                 '--gwasFile "tests/testthat/testOutput/gwasRes.json"',
                 '--adj-method "bonferroni"',
                 '--filter-nPoints 3000',
                 '--outFile "tests/testthat/testOutput/gwasRes_adj.json"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })




  # ldplot ----
  test_that('ldplot', {

    cmd <- paste('Rscript ./r-geno-tools-engine.R ldplot',
                 '--genoFile "data/geno/testMarkerData01.vcf"',
                 '--from 42',
                 '--to 62',
                 '--outFile "tests/testthat/testOutput/ldplot.png"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })



  # relmat-ped ----
  test_that('relmat-ped', {

    cmd <- paste('Rscript ./r-geno-tools-engine.R relmat-ped',
                 '--pedFile "data/pedigree/testPedData_char.csv"',
                 '--outFile "tests/testthat/testOutput/pedRelMat.json"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  # relmat-geno ----
  test_that('relmat-geno', {

    cmd <- paste('Rscript ./r-geno-tools-engine.R relmat-geno',
                 '--genoFile "data/geno/breedGame_geno.vcf.gz"',
                 '--outFile "tests/testthat/testOutput/genoRelMat.json"')
    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

  # relmat-combined ----
  test_that('relmat-combined', {

    cmd <- paste('Rscript ./r-geno-tools-engine.R relmat-combined',
                 '--ped-relmatFile data/results/breedGame_pedRelMat.csv',
                 '--geno-relmatFile data/results/breedGame_genoRelMat.csv',
                 # '--combine-method Legarra',
                 '--combine-method Martini',
                 '--tau 1.3',
                 '--omega 0.7',
                 '--outFile "tests/testthat/testOutput/genoRelMat.json"')
    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })




  # relmat-heatmap ----
  test_that('relmat-heatmap', {

    cmd <- paste('Rscript ./r-geno-tools-engine.R relmat-heatmap',
                 '--relmatFile "tests/testthat/testOutput/pedRelMat.json"',
                 '--outFile "tests/testthat/testOutput/relMat.png"')

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })



  # pedNetwork ----
  test_that('pedNetwork', {

    cmd <- paste('Rscript ./r-geno-tools-engine.R pedNetwork',
                 '--pedFile "data/pedigree/testPedData_char.csv"',
                 '--outFile "tests/testthat/testOutput/pedNet.html"')


    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })



  # crossing-simulation ----
  test_that('crossing-simulation', {

    cmd <- paste(
      'Rscript ./r-geno-tools-engine.R crossing-simulation',
      '--genoFile "data/geno/breedGame_phasedGeno.vcf.gz"',
      '--crossTableFile "data/crossingTable/breedGame_crossTable.csv"',
      '--SNPcoordFile "data/SNPcoordinates/breedingGame_SNPcoord.csv"',
      '--outFile "tests/testthat/testOutput/crossSim.vcf.gz"'
    )

    expect_error({
      x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
    }, NA)
    expect_equal(x, 0)
  })

})
