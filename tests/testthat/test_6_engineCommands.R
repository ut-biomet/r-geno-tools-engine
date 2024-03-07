# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Test engine commands





if (!nzchar(Sys.getenv("RGENOROOT"))) {
  Sys.setenv(RGENOROOT = normalizePath('../..'))
}
if (!nzchar(Sys.getenv("ROOT_DATA_FILES"))) {
  Sys.setenv(ROOT_DATA_FILES = Sys.getenv("RGENOROOT"))
}

# clean testOutput dir
f <- list.files(paste0(Sys.getenv("ROOT_DATA_FILES"),
                       "/tests/testthat/testOutput"),
                full.names = TRUE,
                include.dirs = TRUE)
invisible(lapply(f, unlink, recursive = TRUE))


rGenoCommand <- paste0(Sys.getenv("RGENOROOT"),
                       "/r-geno-tools-engine.R")

# help ----
test_that('engine help', {
  expect_no_warning({
    x <- system(rGenoCommand,
                intern = FALSE,
                ignore.stdout = TRUE)
  })
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
    expect_no_warning({
      x <- system(paste(rGenoCommand, cmd),
                  intern = FALSE,
                  ignore.stdout = TRUE)
    })
    expect_equal(x, 0)
  })
}

# gwas ----
test_that('gwas', {
  cmd <- paste(rGenoCommand, 'gwas',
               '--genoFile "$ROOT_DATA_FILES/data/geno/testMarkerData01.vcf.gz"',
               '--phenoFile "$ROOT_DATA_FILES/data/pheno/testPhenoData01.csv"',
               '--trait "Flowering.time.at.Arkansas"',
               '--test "score"',
               '--response "quantitative"',
               '--thresh-maf 0.05',
               '--thresh-callrate 0.95',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/gwasRes_score.json"')

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})

test_that('gwas', {
  cmd <- paste(rGenoCommand, 'gwas',
               '--genoFile "$ROOT_DATA_FILES/data/geno/testMarkerData01.vcf.gz"',
               '--phenoFile "$ROOT_DATA_FILES/data/pheno/testPhenoData01.csv"',
               '--trait "Flowering.time.at.Arkansas"',
               '--test "wald"',
               '--fixed 0',
               '--response "quantitative"',
               '--thresh-maf 0.05',
               '--thresh-callrate 0.95',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/gwasRes_wald.json"')

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})


# gwas-manplot ----
test_that('gwas-manplot', {

  # html
  cmd <- paste(rGenoCommand, 'gwas-manplot',
               '--gwasFile "$ROOT_DATA_FILES/tests/testthat/testOutput/gwasRes_score.json"',
               '--adj-method "bonferroni"',
               '--thresh-p 0.05',
               '--filter-nPoints 3000',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/manPlot.html"')

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)

  # png
  cmd <- paste(rGenoCommand, 'gwas-manplot',
               '--gwasFile "$ROOT_DATA_FILES/tests/testthat/testOutput/gwasRes_score.json"',
               '--adj-method "bonferroni"',
               '--thresh-p 0.05',
               '--filter-nPoints 3000',
               '--no-interactive',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/manPlot.png"')

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})




# gwas-adjResults ----
test_that('gwas-adjResults', {

  cmd <- paste(rGenoCommand, 'gwas-adjresults',
               '--gwasFile "$ROOT_DATA_FILES/tests/testthat/testOutput/gwasRes_score.json"',
               '--adj-method "bonferroni"',
               '--filter-nPoints 3000',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/gwasRes_adj.json"')

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})




# ldplot ----
test_that('ldplot', {

  cmd <- paste(rGenoCommand, 'ldplot',
               '--genoFile "$ROOT_DATA_FILES/data/geno/testMarkerData01.vcf.gz"',
               '--from 42',
               '--to 62',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/ldplot.png"')

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})



# relmat-ped ----
test_that('relmat-ped', {

  cmd <- paste(rGenoCommand, 'relmat-ped',
               '--pedFile "$ROOT_DATA_FILES/data/pedigree/testPedData_char.csv"',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/pedRelMat.json"')


  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})

# relmat-geno ----
test_that('relmat-geno', {

  cmd <- paste(rGenoCommand, 'relmat-geno',
               '--genoFile "$ROOT_DATA_FILES/data/geno/breedGame_geno.vcf.gz"',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/genoRelMat.json"')
  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})

# relmat-combined ----
test_that('relmat-combined', {

  cmd <- paste(rGenoCommand, 'relmat-combined',
               '--ped-relmatFile $ROOT_DATA_FILES/data/results/breedGame_pedRelMat.csv',
               '--geno-relmatFile $ROOT_DATA_FILES/data/results/breedGame_genoRelMat.csv',
               # '--combine-method Legarra',
               '--combine-method Martini',
               '--tau 1.3',
               '--omega 0.7',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/combinedRelMat.json"')
  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})




# relmat-heatmap ----
test_that('relmat-heatmap', {

  cmd <- paste(rGenoCommand, 'relmat-heatmap',
               '--relmatFile "$ROOT_DATA_FILES/tests/testthat/testOutput/pedRelMat.json"',
               '--no-interactive',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/relMat.png"')

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})



# pedNetwork ----
test_that('pedNetwork', {

  cmd <- paste(rGenoCommand, 'pedNetwork',
               '--pedFile "$ROOT_DATA_FILES/data/pedigree/testPedData_char.csv"',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/pedNet.html"')


  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})



# crossing-simulation ----
test_that('crossing-simulation', {

  cmd <- paste(
    rGenoCommand, 'crossing-simulation',
    '--genoFile "$ROOT_DATA_FILES/data/geno/breedGame_phasedGeno.vcf.gz"',
    '--crossTableFile "$ROOT_DATA_FILES/data/crossingTable/breedGame_crossTable.csv"',
    '--SNPcoordFile "$ROOT_DATA_FILES/data/SNPcoordinates/breedingGame_SNPcoord.csv"',
    '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/crossSim.vcf.gz"'
  )

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})

# progeny-blup-calculation ----
test_that('progeny-blup-calculation', {
  cmd <- paste(
    rGenoCommand, 'progeny-blup-calculation',
    '--genoFile "$ROOT_DATA_FILES/data/geno/breedGame_phasedGeno.vcf.gz"',
    '--crossTableFile "$ROOT_DATA_FILES/data/crossingTable/breedGame_small_crossTable.csv"',
    '--SNPcoordFile "$ROOT_DATA_FILES/data/SNPcoordinates/breedingGame_SNPcoord.csv"',
    '--markerEffectsFile "$ROOT_DATA_FILES/data/markerEffects/breedGame_markerEffects.csv"',
    '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/progBlups.json"'
  )

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})

test_that('progeny-blup-calculation (multi-trait marker effect)', {
  cmd <- paste(
    rGenoCommand, 'progeny-blup-calculation',
    '--genoFile "$ROOT_DATA_FILES/data/geno/breedGame_phasedGeno.vcf.gz"',
    '--crossTableFile "$ROOT_DATA_FILES/data/crossingTable/breedGame_small_crossTable.csv"',
    '--SNPcoordFile "$ROOT_DATA_FILES/data/SNPcoordinates/breedingGame_SNPcoord.csv"',
    '--markerEffectsFile "$ROOT_DATA_FILES/data/markerEffects/breedGame_markerEffects_2traits.json"',
    '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/progBlups_2traits.json"'
  )

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})


# progeny-blup-plot ----
test_that('progeny-blup-plot', {
  cmd <- paste(
    rGenoCommand, 'progeny-blup-plot',
    '--progeniesBlupFile "$ROOT_DATA_FILES/tests/testthat/testOutput/progBlups.json"',
    '--y-axis-name "Phenotypic trait"',
    '--error-bar-interval 0.95',
    '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/blupPlot_1trait.html"'
  )

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})

test_that('progeny-blup-plot (multi-trait file)', {
  cmd <- paste(
    rGenoCommand, 'progeny-blup-plot',
    '--progeniesBlupFile "$ROOT_DATA_FILES/tests/testthat/testOutput/progBlups_2traits.json"',
    '--y-axis-name "Phenotypic trait"',
    '--error-bar-interval 0.95',
    '--trait "trait1"',
    '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/blupPlot_1trait_with_multi-trait-input.html"'
  )

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})

test_that('progeny-blup-plot-2-traits', {
  cmd <- paste(
    rGenoCommand, 'progeny-blup-plot-2-traits',
    '--progeniesBlupFile "$ROOT_DATA_FILES/tests/testthat/testOutput/progBlups_2traits.json"',
    '--x-trait "trait1"',
    '--y-trait "trait2"',
    '--confidence-level 0.95',
    '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/blupPlot_2traits.html"'
  )

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE)
  })
  expect_equal(x, 0)
})


test_that('Expected error RAW', {

  cmd <- paste(rGenoCommand, 'pedNetwork',
               '--pedFile "/doNotExist"',
               '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/pedNet.html"')

  expect_no_warning({
    x <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  })
  expect_equal(x, 1)
})


test_that('Expected error JSON', {
  args <- paste('pedNetwork',
                '--pedFile "/doNotExist"',
                '--outFile "$ROOT_DATA_FILES/tests/testthat/testOutput/pedNet.html"',
                '--json-errors')

  tmp_stderr <- tempfile(fileext = ".json")
  tmp_stdout <- tempfile(fileext = ".log")

  expect_no_warning({
    x <- system2(rGenoCommand, args, stderr = tmp_stderr, stdout = tmp_stdout)
  })

  expect_equal(x, 42)

  expect_no_error({
    jsonlite::fromJSON(tmp_stderr)
  })

  unlink(tmp_stderr)
  unlink(tmp_stdout)
})
