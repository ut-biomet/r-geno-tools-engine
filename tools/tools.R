# Set of utility functions to use for this projects
# those functions are not meant to be integrated in the engine
# but for development purpose only



#' Write engine's documentation
#'
#' @param srcDir path of R sources folder (default "./src")
#' @param docDir path of documentation folder (default "./doc")
#'
#' @return NULL
#' @details Write engine's functions documentation in a README.md file located in `docDir`.
writeDoc <- function(srcDir = "./src",
                     docDir = "./doc"){
  stopifnot(dir.exists(srcDir))
  stopifnot(dir.exists(docDir))
  tmpDir <- file.path(docDir, ".tmp")
  dir.create(tmpDir, showWarnings = FALSE)

  outFile <- file.path(docDir, "README.md")
  file.remove(outFile)
  file.create(outFile)

  sourcefiles <- list.files(srcDir, pattern = ".R$", full.names = TRUE)
  for (sourcefile in sourcefiles) {
    source_env = roxygen2::env_file(sourcefile)
    source(sourcefile, local = source_env, keep.source = TRUE)
    rd_blocks = roxygen2::parse_file(sourcefile, env = source_env)
    help_topics = roxygen2::roclet_process(roxygen2::rd_roclet(),
                                           rd_blocks,
                                           env = source_env,
                                           # env = NULL,
                                           dirname(sourcefile))
    rd_code = lapply(help_topics, format)

    for (funName in names(rd_code)) {
      writeLines(rd_code[[funName]],
                 con = file.path(tmpDir, funName))
      Rd2md::Rd2markdown(file.path(tmpDir, funName), outFile, append = TRUE)
      file.remove(file.path(tmpDir, funName))
    }

    # log functions with missing doc:
    definedFunctions <- names(source_env)
    definedFunctions <- definedFunctions[definedFunctions != '.packageName']
    documentedFunctions <- tools::file_path_sans_ext(names(rd_code))
    missIds <- which(!definedFunctions %in% documentedFunctions)
    if (any(missIds)) {
    message("Missing roxygen doc in file ", sourcefile, ":\n\t`",
            paste(definedFunctions[missIds], collapse = '`, `'), "`")
    }
  }
  unlink(tmpDir, recursive = TRUE)

  NULL
}


#' Create results example files
#'
#'
#' @return NULL
#' @details Use predefined input/output files
createResultExample <- function() {

  invisible(
    sapply(FUN = source,
      X = list.files("./src", pattern = ".R$", full.names = TRUE))
  )
  set.seed(42)

  # create gwas results ----
  cat('create gwas results ----\n')
  genoFile <- 'data/geno/testMarkerData01.vcf'
  phenoFile <- 'data/pheno/testPhenoData01.csv'

  gwas <- run_gwas(genoFile,
                   phenoFile,
                   trait = "Flowering.time.at.Arkansas",
                   test = 'score',
                   thresh_maf = 0.05,
                   thresh_callrate = 0.5,
                   outFile = 'data/results/gwasResult.json')

  # create adjusted results ----
  cat('create adjusted results ----\n')
  gwasFile <- gwas$file

  gwas <- run_resAdjustment(gwasFile,
                            adj_method = "bonferroni",
                            outFile = 'data/results/ajdResults.json')

  # create Manhattan plot  ----
  cat('create Manhattan plot  ----\n')
  gwasFile <- gwas$file

  draw_manhattanPlot(gwasFile,
                     adj_method = 'bonferroni',
                     thresh_p = 0.05,
                     interactive = TRUE,
                     outFile = 'data/results/manplot.html')
  draw_manhattanPlot(gwasFile,
                     adj_method = 'bonferroni',
                     thresh_p = 0.05,
                     interactive = FALSE,
                     outFile = 'data/results/manplot.png')


  # create LdPlot ----
  cat('create LdPlot ----\n')
  draw_ldPlot(genoFile,
              from = 42,
              to = 62,
              outFile = 'data/results/ldplot.png')

  # create pedigree network ----
  cat('create pedigree network ----\n')
  pedFile <- 'data/pedigree/testPedData_char.csv'

  draw_pedNetwork(pedFile,
                  outFile = 'data/results/pedigreeNetwork.html')

  # create pedigree relationship ----
  cat('create pedigree relationship ----\n')
  calc_pedRelMat(pedFile,
                 outFile = 'data/results/pedigreeRelationship.csv')
  calc_pedRelMat(pedFile,
                 outFile = 'data/results/pedigreeRelationship.json')

  calc_pedRelMat(pedFile = 'data/pedigree/breedGame_pedigree.csv',
                 outFile = 'data/results/breedGame_pedRelMat.csv')
  calc_pedRelMat(pedFile = 'data/pedigree/breedGame_pedigree.csv',
                 outFile = 'data/results/breedGame_pedRelMat.json')

  # create genomic relationship ----
  cat('create genomic relationship ----\n')
  calc_genoRelMat(genoFile = 'data/geno/breedGame_geno.vcf.gz',
                  outFile = 'data/results/breedGame_genoRelMat.csv')
  calc_genoRelMat(genoFile = 'data/geno/breedGame_geno.vcf.gz',
                  outFile = 'data/results/breedGame_genoRelMat.json')

  # create combined relationship ----
  cat('create combined relationship ----\n')
  calc_combinedRelMat(pedRelMatFile = 'data/results/breedGame_pedRelMat.csv',
                      genoRelMatFile = 'data/results/breedGame_genoRelMat.csv',
                      method = 'Legarra',
                      outFile = 'data/results/breedGame_combinedRelMat.csv')

  # create relationship heatmap ----
  cat('create relationship heatmap ----\n')
  draw_relHeatmap(relMatFile = 'data/results/pedigreeRelationship.csv',
                  interactive = TRUE,
                  outFile = 'data/results/relationshipHeatmap.html')
  draw_relHeatmap(relMatFile = 'data/results/pedigreeRelationship.csv',
                  interactive = FALSE,
                  outFile = 'data/results/relationshipHeatmap.png')


  # create progeny blup estimation results ----
  cat('create progeny blup estimation results ----\n')
  genoFile <- 'data/geno/breedGame_phasedGeno.vcf.gz'
  crossTableFile <- 'data/crossingTable/breedGame_small_crossTable.csv'
  SNPcoordFile <- 'data/SNPcoordinates/breedingGame_SNPcoord.csv'
  markerEffectsFile <- 'data/markerEffects/breedGame_markerEffects.csv'
  outFile <- 'data/results/progenyBlupEstimation.json'

  calc_progenyBlupEstimation(
    genoFile = genoFile,
    crossTableFile = crossTableFile,
    SNPcoordFile = SNPcoordFile,
    markerEffectsFile = markerEffectsFile,
    outFile = outFile
  )

  outFile_plot <- 'data/results/progenyBlupEstimation_plot.html'
  draw_progBlupsPlot(progEstimFile = outFile,
                     errorBarInterval= 0.95,
                     y_axisName = "Genetic values",
                     sorting = 'alpha',
                     trait = "trait1",
                     outFile = outFile_plot)

  markerEffectsFile <- 'data/markerEffects/breedGame_markerEffects_2traits.json'
  outFile <- 'data/results/progenyBlupEstimation_2traits.json'
  calc_progenyBlupEstimation(
    genoFile = genoFile,
    crossTableFile = crossTableFile,
    SNPcoordFile = SNPcoordFile,
    markerEffectsFile = markerEffectsFile,
    outFile = outFile
  )


  outFile_plot <- 'data/results/progenyBlupEstimation_plot_2traits.html'
  draw_progBlupsPlot_2traits(progEstimFile = outFile,
                             x_trait = 'trait1',
                             y_trait = 'trait2',
                             confidenceLevel = 0.95,
                             x_suffix = "",
                             y_suffix = "",
                             ellipses_npoints = 100,
                             outFile = outFile_plot)



  # clear folder
  unlink(list.dirs('data/results', recursive = FALSE),
         recursive = TRUE)

  NULL
}
