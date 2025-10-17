#!/usr/bin/env Rscript

# Author: Julien Diot juliendiot@ut-biomet.org
# Tue 24 Aug 2021 The University of Tokyo
#
# Description:
# Main script of R-geno-tools-enigne


options(keep.source = TRUE) # to get the error locations
ERROR_AS_JSON <- FALSE

if (nzchar(Sys.getenv("RGENOROOT"))) {
  projectRoot <- Sys.getenv("RGENOROOT")
} else {
  projectRoot <- getwd()
}

if ("renv" %in% rownames(installed.packages()) && !isTRUE(as.logical(Sys.getenv("NO_RENV")))) {
  renv::load(projectRoot)
}

withCallingHandlers(
  engineError = function(err) {
    if (ERROR_AS_JSON) {
      jsonErr <- jsonlite::toJSON(
        list(
          message = err$message,
          extra = err$extra
        ),
        auto_unbox = TRUE
      )
      write(jsonErr, stderr())
      q("no", status = 42)
    } else {
      stop(err)
    }
  },
  error = function(err) {
    if (ERROR_AS_JSON) {
      # extract stack-trace:
      error_locations <- limitedLabels(sys.calls())
      error_locations <- error_locations[-1] # drop 1st
      error_locations[1] <- paste0("r-geno-tools-engine.R", error_locations[1])

      start_with_R_file <- grepl("^[\\w,\\s-]+\\.[Rr]", error_locations, perl = TRUE)
      error_locations <- error_locations[start_with_R_file]
      last_loc <- error_locations[length(error_locations)]
      last_fun <- regmatches(
        last_loc,
        gregexec("^.*: ([\\w.]+)\\(",
          last_loc,
          perl = TRUE
        )
      )[[1]][2, 1] # extract catching group
      if (last_fun == ".handleSimpleError") {
        last_loc <- gsub("\\.handleSimpleError.*", replacement = deparse(err$call), last_loc)
      }
      error_locations[length(error_locations)] <- last_loc
      error_locations <- paste0(error_locations, collapse = "\n")

      # write error info in json
      jsonErr <- jsonlite::toJSON(
        list(
          message = err$message,
          extra = list(
            "code" = "UNEXPECTED_ERROR",
            "stacktrace" = error_locations
          )
        ),
        auto_unbox = TRUE
      )
      write(jsonErr, stderr())
      q("no", status = 42)
    } else {
      stop(err)
    }
  },
  warning = function(warn) {
    # catch warnings and send them to stdout instead of stderr
    cat("-- WARNING -- ", warn$message, "\n")
    rlang::cnd_muffle(warn)
  },
  {
    # Create parser ----
    library(argparse)
    formatter_class <- "lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=42, width=100)"
    main_parser <- ArgumentParser(
      prog = "r-geno-tools-engine",
      description = "Toolbox for some genomic analysis written in R.",
      formatter_class = formatter_class,
      # argument_default = "True",
      usage = "%(prog)s"
    )
    main_subparsers <- main_parser$add_subparsers(
      prog = "r-geno-tools-engine",
      title = "Available commands",
      dest = "command"
    )

    # load argument descriptions
    source(file.path(projectRoot, "src/commandArgs.R"))
    parserList <- c() # to keep names of subcommands

    ## GWAS gwas ----
    gwas_parser <- main_subparsers$add_parser("gwas",
      help = "Do a gwas analysis",
      description = "Run a GWAS analysis",
      formatter_class = formatter_class,
      argument_default = "True"
    )
    gwas_parser$add_argument(arg$genoFile$flag,
      help = arg$genoFile$help,
      type = arg$genoFile$type,
      required = TRUE
    )
    gwas_parser$add_argument(arg$phenoFile$flag,
      help = arg$phenoFile$help,
      type = arg$phenoFile$type,
      required = TRUE
    )
    gwas_parser$add_argument(arg$trait$flag,
      help = arg$trait$help,
      type = arg$trait$type,
      required = TRUE
    )
    gwas_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    gwas_parser$add_argument(arg$test$flag,
      help = arg$test$help,
      type = arg$test$type,
      default = arg$test$default
    )
    gwas_parser$add_argument(arg$fixed$flag,
      help = arg$fixed$help,
      type = arg$fixed$type,
      default = arg$fixed$default
    )
    gwas_parser$add_argument(arg$response$flag,
      help = arg$response$help,
      type = arg$response$type,
      default = arg$response$default
    )
    gwas_parser$add_argument(arg$thresh_maf$flag,
      help = arg$thresh_maf$help,
      type = arg$thresh_maf$type,
      default = arg$thresh_maf$default
    )
    gwas_parser$add_argument(arg$thresh_callrate$flag,
      help = arg$thresh_callrate$help,
      type = arg$thresh_callrate$type,
      default = arg$thresh_callrate$default
    )
    gwas_parser$add_argument(arg$n_markers$flag,
      help = arg$n_markers$help,
      type = arg$n_markers$type,
      default = arg$n_markers$default
    )
    gwas_parser$add_argument(arg$n_markers_tolerance$flag,
      help = arg$n_markers_tolerance$help,
      type = arg$n_markers_tolerance$type,
      default = arg$n_markers_tolerance$default
    )
    gwas_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## GWAS manplot ----
    manplot_parser <- main_subparsers$add_parser("gwas-manplot",
      help = "Draw a Manhattan Plot",
      description = "Draw a Manhattan Plot",
      formatter_class = formatter_class
    )
    manplot_parser$add_argument(arg$gwasFile$flag,
      help = arg$gwasFile$help,
      type = arg$gwasFile$type,
      required = TRUE
    )
    manplot_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    manplot_parser$add_argument(arg$adj_method$flag,
      help = arg$adj_method$help,
      type = arg$adj_method$type,
      default = arg$adj_method$default
    )
    manplot_parser$add_argument(arg$thresh_p$flag,
      help = arg$thresh_p$help,
      type = arg$thresh_p$type,
      default = arg$thresh_p$default
    )
    manplot_parser$add_argument(arg$chr$flag,
      help = arg$chr$help,
      type = arg$chr$type,
      default = arg$chr$default
    )
    manplot_parser$add_argument(arg$no_interactive$flag,
      help = arg$no_interactive$help,
      action = "store_true"
    )
    manplot_parser$add_argument(arg$filter_pAdj$flag,
      help = arg$filter_pAdj$help,
      type = arg$filter_pAdj$type,
      default = arg$filter_pAdj$default
    )
    manplot_parser$add_argument(arg$filter_nPoints$flag,
      help = arg$filter_nPoints$help,
      type = arg$filter_nPoints$type,
      default = arg$filter_nPoints$default
    )
    manplot_parser$add_argument(arg$filter_quant$flag,
      help = arg$filter_quant$help,
      type = arg$filter_quant$type,
      default = arg$filter_quant$default
    )
    manplot_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## GWAS adjResults ----
    adjResults_parser <- main_subparsers$add_parser("gwas-adjresults",
      help = "Adjust GWAS p-values",
      description = "Adjust GWAS p-values and filter the results",
      formatter_class = formatter_class
    )
    adjResults_parser$add_argument(arg$gwasFile$flag,
      help = arg$gwasFile$help,
      type = arg$gwasFile$type,
      required = TRUE
    )
    adjResults_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    adjResults_parser$add_argument(arg$adj_method$flag,
      help = arg$adj_method$help,
      type = arg$adj_method$type,
      default = arg$adj_method$default
    )
    adjResults_parser$add_argument(arg$filter_pAdj$flag,
      help = arg$filter_pAdj$help,
      type = arg$filter_pAdj$type,
      default = arg$filter_pAdj$default
    )
    adjResults_parser$add_argument(arg$filter_nPoints$flag,
      help = arg$filter_nPoints$help,
      type = arg$filter_nPoints$type,
      default = arg$filter_nPoints$default
    )
    adjResults_parser$add_argument(arg$filter_quant$flag,
      help = arg$filter_quant$help,
      type = arg$filter_quant$type,
      default = arg$filter_quant$default
    )
    adjResults_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )



    ## LDplot ----
    ldplot_parser <- main_subparsers$add_parser("ldplot",
      help = "Get plot of the linkage disequilibrium between consecutive markers",
      description = "Get plot of the linkage disequilibrium between consecutive markers",
      formatter_class = formatter_class
    )
    ldplot_parser$add_argument(arg$genoFile$flag,
      help = arg$genoFile$help,
      type = arg$genoFile$type,
      required = TRUE
    )
    ldplot_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    ldplot_parser$add_argument(arg$from$flag,
      help = arg$from$help,
      type = arg$from$type,
      required = TRUE
    )
    ldplot_parser$add_argument(arg$to$flag,
      help = arg$to$help,
      type = arg$to$type,
      required = TRUE
    )
    ldplot_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )



    ## relmat pedigree ----
    pedRelMat_parser <- main_subparsers$add_parser("relmat-ped",
      help = "Pedigree relationship matrix",
      description = "Caclulate a pedigree relationship matrix.\\n\\nThe output file shoud end by either `.csv` or `.json`.",
      formatter_class = formatter_class
    )
    pedRelMat_parser$add_argument(arg$pedFile$flag,
      help = arg$pedFile$help,
      type = arg$pedFile$type,
      required = TRUE
    )
    pedRelMat_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    pedRelMat_parser$add_argument(arg$u_string$flag,
      help = arg$u_string$help,
      type = arg$u_string$type,
      default = arg$u_string$default
    )
    pedRelMat_parser$add_argument(arg$no_header$flag,
      help = arg$no_header$help,
      action = "store_true"
    )
    pedRelMat_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## relmat genomic ----
    genoRelMat_parser <- main_subparsers$add_parser("relmat-geno",
      help = "Genomic relationship matrix",
      description = "Caclulate a genomic relationship matrix.\\n\\nThe output file shoud end by either `.csv` or `.json`.",
      formatter_class = formatter_class
    )
    genoRelMat_parser$add_argument(arg$genoFile$flag,
      help = arg$genoFile$help,
      type = arg$genoFile$type,
      required = TRUE
    )
    genoRelMat_parser$add_argument(arg$n_markers$flag,
      help = arg$n_markers$help,
      type = arg$n_markers$type,
      default = arg$n_markers$default
    )
    genoRelMat_parser$add_argument(arg$n_markers_tolerance$flag,
      help = arg$n_markers_tolerance$help,
      type = arg$n_markers_tolerance$type,
      default = arg$n_markers_tolerance$default
    )
    genoRelMat_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    genoRelMat_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )


    ## relmat combined ----
    combinedRelMat_parser <- main_subparsers$add_parser(
      "relmat-combined",
      help = "Combined relationship matrix",
      description = paste(
        "Correct a pedigree relationship matrix using a genomic relationship matrix.",
        "\\n\\nThe output file shoud end by either `.csv` or `.json`.",
        "\\n\\nIndividuals in the genomic relationship matrix which do not belong ",
        "to the pedigree relationship matrix will be ignored.",
        "\\n\\nRefference of the method can be found at:",
        "\\n\t- \"Martini, JW, et al. 2018 The effect of the H-1 scaling factors tau and omega on the structure of H in the single-step procedure. Genetics Selection Evolution 50(1), 16\"",
        "\\n\t-\"Legarra, A, et al. 2009 A relationship matrix including full pedigree and genomic information. Journal of Dairy Science 92, 4656–4663\"",
        "\\n \\n \\n "
      ),
      formatter_class = formatter_class
    )
    combinedRelMat_parser$add_argument(arg$ped_relmatFile$flag,
      help = arg$ped_relmatFile$help,
      type = arg$ped_relmatFile$type,
      required = TRUE
    )
    combinedRelMat_parser$add_argument(arg$geno_relmatFile$flag,
      help = arg$geno_relmatFile$help,
      type = arg$geno_relmatFile$type,
      required = TRUE
    )
    combinedRelMat_parser$add_argument(arg$combine_method$flag,
      help = arg$combine_method$help,
      type = arg$combine_method$type,
      required = TRUE
    )
    combinedRelMat_parser$add_argument(arg$tau$flag,
      help = arg$tau$help,
      type = arg$tau$type
    )
    combinedRelMat_parser$add_argument(arg$omega$flag,
      help = arg$omega$help,
      type = arg$omega$type
    )
    combinedRelMat_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    combinedRelMat_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )


    ## relmat heatmap ----
    relMatHeat_parser <- main_subparsers$add_parser("relmat-heatmap",
      help = "Draw a heatmap of a relationship matrix",
      description = "Draw a heatmap of a relationship matrix",
      formatter_class = formatter_class
    )
    relMatHeat_parser$add_argument(arg$relmatFile$flag,
      help = arg$relmatFile$help,
      type = arg$relmatFile$type,
      required = TRUE
    )
    relMatHeat_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    relMatHeat_parser$add_argument(arg$no_interactive$flag,
      help = arg$no_interactive$help,
      action = "store_true"
    )
    relMatHeat_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## pedNetwork ----
    pedNet_parser <- main_subparsers$add_parser("pedNetwork",
      help = "Draw interactive pedigree network",
      description = "Draw interactive pedigree network",
      formatter_class = formatter_class
    )
    pedNet_parser$add_argument(arg$pedFile$flag,
      help = arg$pedFile$help,
      type = arg$pedFile$type,
      required = TRUE
    )
    pedNet_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    pedNet_parser$add_argument(arg$u_string$flag,
      help = arg$u_string$help,
      type = arg$u_string$type,
      default = arg$u_string$default
    )
    pedNet_parser$add_argument(arg$no_header$flag,
      help = arg$no_header$help,
      action = "store_true"
    )
    pedNet_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )


    ## crossing-simulation ----
    crossSim_parser <- main_subparsers$add_parser("crossing-simulation",
      help = "Simulate crossing between individuals",
      description = paste(
        "Simulate the cross",
        "specified in the crossing",
        "table and write the simulated",
        "genotypes in a phased `.vcf.gz`",
        "file. This command will return",
        "the path of the created file"
      ),
      formatter_class = formatter_class
    )
    crossSim_parser$add_argument(arg$phasedGenoFile$flag,
      help = arg$phasedGenoFile$help,
      type = arg$phasedGenoFile$type,
      required = TRUE
    )
    crossSim_parser$add_argument(arg$crossTableFile$flag,
      help = arg$crossTableFile$help,
      type = arg$crossTableFile$type,
      required = TRUE
    )
    crossSim_parser$add_argument(arg$SNPcoordFile$flag,
      help = arg$SNPcoordFile$help,
      type = arg$SNPcoordFile$type,
      default = NULL
    )
    # required = TRUE)
    crossSim_parser$add_argument(arg$nCross$flag,
      help = arg$nCross$help,
      type = arg$nCross$type,
      default = arg$nCross$default
    )
    crossSim_parser$add_argument(arg$outGenoFile$flag,
      help = arg$outGenoFile$help,
      type = arg$outGenoFile$type,
      required = TRUE
    )
    crossSim_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )


    ## progeny-blup-calculation ----
    prog_blup_parser <- main_subparsers$add_parser(
      "progeny-blup-calculation",
      help = "Estimate progenies' BLUPs of given crosses",
      description = paste(
        "Estimate the BLUPs' expected value and variance of the progenies of a",
        "given crosses specifyed in the crossing table. The results are written in",
        "a `.json` file."
      ),
      formatter_class = formatter_class
    )

    prog_blup_parser$add_argument(arg$phasedGenoFile$flag,
      help = arg$phasedGenoFile$help,
      type = arg$phasedGenoFile$type,
      required = TRUE
    )
    prog_blup_parser$add_argument(arg$crossTableFile$flag,
      help = arg$crossTableFile$help,
      type = arg$crossTableFile$type,
      required = TRUE
    )
    prog_blup_parser$add_argument(arg$SNPcoordFile$flag,
      help = arg$SNPcoordFile$help,
      type = arg$SNPcoordFile$type,
      required = TRUE
    )
    prog_blup_parser$add_argument(arg$markersEffectsFile$flag,
      help = arg$markersEffectsFile$help,
      type = arg$markersEffectsFile$type,
      required = TRUE
    )
    prog_blup_parser$add_argument(arg$outProgBlupFile$flag,
      help = arg$outProgBlupFile$help,
      type = arg$outProgBlupFile$type,
      required = TRUE
    )
    prog_blup_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )


    ## progeny-blup-plot ----
    prog_blup_plot <- main_subparsers$add_parser(
      "progeny-blup-plot",
      help = "Draw a plot of the progenies BLUPs' expected values with error bars",
      description = paste(
        "Draw a plot of the progenies BLUPs' expected values with error bars.",
        "X axis is the crosses, and Y axis the blups. The points are located",
        "at the expected value and the error bar length is the standard deviation."
      ),
      formatter_class = formatter_class
    )

    prog_blup_plot$add_argument(arg$progeniesBlupFile$flag,
      help = arg$progeniesBlupFile$help,
      type = arg$progeniesBlupFile$type,
      required = TRUE
    )
    prog_blup_plot$add_argument(arg$sorting$flag,
      help = arg$sorting$help,
      type = arg$sorting$type,
      default = arg$sorting$default
    )
    prog_blup_plot$add_argument(arg$errorBarInterval$flag,
      help = arg$errorBarInterval$help,
      type = arg$errorBarInterval$type,
      default = arg$errorBarInterval$default
    )
    prog_blup_plot$add_argument(arg$y_axisName$flag,
      help = arg$y_axisName$help,
      type = arg$y_axisName$type,
      default = arg$y_axisName$default
    )
    prog_blup_plot$add_argument(arg$trait$flag,
      help = arg$trait$help,
      type = arg$trait$type,
      default = ""
    )
    prog_blup_plot$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    prog_blup_plot$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## progeny-blup-plot-2-traits ----
    prog_blup_plot_2_traits <- main_subparsers$add_parser(
      "progeny-blup-plot-2-traits",
      help = "Draw a plot of the progenies BLUPs' expected values for 2 traits with prediction ellipses",
      description = paste(
        "Draw a plot of the progenies BLUPs' expected values with prediction ellipses.",
        "The points are located at the expected value and the ellipses",
        "size represent the `--confidenceLevel` prediction interval (default 95%)."
      ),
      formatter_class = formatter_class
    )

    prog_blup_plot_2_traits$add_argument(arg$progeniesBlupFile$flag,
      help = arg$progeniesBlupFile$help,
      type = arg$progeniesBlupFile$type,
      required = TRUE
    )
    prog_blup_plot_2_traits$add_argument(arg$x_trait$flag,
      help = arg$x_trait$help,
      type = arg$x_trait$type,
      required = TRUE
    )
    prog_blup_plot_2_traits$add_argument(arg$y_trait$flag,
      help = arg$y_trait$help,
      type = arg$y_trait$type,
      required = TRUE
    )
    prog_blup_plot_2_traits$add_argument(arg$confidenceLevel$flag,
      help = arg$confidenceLevel$help,
      type = arg$confidenceLevel$type,
      default = arg$confidenceLevel$default
    )
    prog_blup_plot_2_traits$add_argument(arg$x_suffix$flag,
      help = arg$x_suffix$help,
      type = arg$x_suffix$type,
      default = arg$x_suffix$default
    )
    prog_blup_plot_2_traits$add_argument(arg$y_suffix$flag,
      help = arg$y_suffix$help,
      type = arg$y_suffix$type,
      default = arg$y_suffix$default
    )
    prog_blup_plot_2_traits$add_argument(arg$ellipses_npoints$flag,
      help = arg$ellipses_npoints$help,
      type = arg$ellipses_npoints$type,
      default = arg$ellipses_npoints$default
    )
    prog_blup_plot_2_traits$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    prog_blup_plot_2_traits$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## GS model evaluation ----
    gs_model_evaluation_parser <- main_subparsers$add_parser(
      "evaluate-gs-model",
      help = "Evaluate a genomic selection model.",
      description = paste(
        "Evaluate a genomic selection model.",
        "This evaluation is done with a repeated cross-validation."
      ),
      formatter_class = formatter_class
    )

    gs_model_evaluation_parser$add_argument(arg$genoFile$flag,
      help = arg$genoFile$help,
      type = arg$genoFile$type,
      required = TRUE
    )
    gs_model_evaluation_parser$add_argument(arg$phenoFile$flag,
      help = arg$phenoFile$help,
      type = arg$phenoFile$type,
      required = TRUE
    )
    gs_model_evaluation_parser$add_argument(arg$trait$flag,
      help = arg$trait$help,
      type = arg$trait$type,
      required = TRUE
    )
    gs_model_evaluation_parser$add_argument(arg$with_dominance$flag,
      help = arg$with_dominance$help,
      type = arg$with_dominance$type,
      required = TRUE
    )
    gs_model_evaluation_parser$add_argument(arg$n_folds$flag,
      help = arg$n_folds$help,
      type = arg$n_folds$type,
      default = arg$n_folds$default
    )
    gs_model_evaluation_parser$add_argument(arg$n_repetitions$flag,
      help = arg$n_repetitions$help,
      type = arg$n_repetitions$type,
      default = arg$n_repetitions$default
    )
    gs_model_evaluation_parser$add_argument(arg$n_markers$flag,
      help = arg$n_markers$help,
      type = arg$n_markers$type,
      default = arg$n_markers$default
    )
    gs_model_evaluation_parser$add_argument(arg$n_markers_tolerance$flag,
      help = arg$n_markers_tolerance$help,
      type = arg$n_markers_tolerance$type,
      default = arg$n_markers_tolerance$default
    )
    gs_model_evaluation_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    gs_model_evaluation_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## plot-gs-model-evaluation ----
    plot_gs_model_evaluation_parser <- main_subparsers$add_parser(
      "plot-gs-model-evaluation",
      help = "Create an interactive plot of the evaluation resutls",
      description = paste(
        "Create an interactive plot of the evaluation resutls"
      ),
      formatter_class = formatter_class
    )


    plot_gs_model_evaluation_parser$add_argument(arg$evaluation_file$flag,
      help = arg$evaluation_file$help,
      type = arg$evaluation_file$type,
      required = TRUE
    )
    plot_gs_model_evaluation_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    plot_gs_model_evaluation_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## Train GS model ----
    train_gs_model_parser <- main_subparsers$add_parser(
      "train-gs-model",
      help = "Train a genomic selection model.",
      description = paste(
        "Train a genomic selection model.",
        "This model allow to include optionnaly dominance effects.",
        "The result file will contain the estimated markers effects that can",
        "be used to make predictions."
      ),
      formatter_class = formatter_class
    )

    train_gs_model_parser$add_argument(arg$genoFile$flag,
      help = arg$genoFile$help,
      type = arg$genoFile$type,
      required = TRUE
    )
    train_gs_model_parser$add_argument(arg$phenoFile$flag,
      help = arg$phenoFile$help,
      type = arg$phenoFile$type,
      required = TRUE
    )
    train_gs_model_parser$add_argument(arg$trait$flag,
      help = arg$trait$help,
      type = arg$trait$type,
      required = TRUE
    )
    train_gs_model_parser$add_argument(arg$with_dominance$flag,
      help = arg$with_dominance$help,
      type = arg$with_dominance$type,
      required = TRUE
    )
    train_gs_model_parser$add_argument(arg$n_markers$flag,
      help = arg$n_markers$help,
      type = arg$n_markers$type,
      default = arg$n_markers$default
    )
    train_gs_model_parser$add_argument(arg$n_markers_tolerance$flag,
      help = arg$n_markers_tolerance$help,
      type = arg$n_markers_tolerance$type,
      default = arg$n_markers_tolerance$default
    )
    train_gs_model_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    train_gs_model_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## gs-predictions ----
    gs_prediction_parser <- main_subparsers$add_parser(
      "gs-predictions",
      help = "Make phenotypic prediction using a markers effects",
      description = paste(
        "Make phenotypic prediction using markers effects"
      ),
      formatter_class = formatter_class
    )

    gs_prediction_parser$add_argument(arg$genoFile$flag,
      help = arg$genoFile$help,
      type = arg$genoFile$type,
      required = TRUE
    )
    gs_prediction_parser$add_argument(arg$markersEffectsFile$flag,
      help = arg$markersEffectsFile$help,
      type = arg$markersEffectsFile$type,
      required = TRUE
    )
    gs_prediction_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    gs_prediction_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )



    ## generate_rnd_marker_effects ----
    generate_rnd_marker_effects_parser <- main_subparsers$add_parser(
      "generate-rnd-marker-effects",
      help = "Generate random marker effects.",
      description = "Generates additive and dominance marker effects for a set of SNP markers. The generated effects are saved as a JSON file. The command print the distribution of the genetic values for the given genotype file with the generated marker effects.",
      formatter_class = formatter_class
    )

    generate_rnd_marker_effects_parser$add_argument(arg$genoFile$flag,
      help = arg$genoFile$help,
      type = arg$genoFile$type,
      required = TRUE
    )
    generate_rnd_marker_effects_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    generate_rnd_marker_effects_parser$add_argument(arg$rate_add$flag,
      help = arg$rate_add$help,
      type = arg$rate_add$type,
      default = arg$rate_add$default,
      required = FALSE
    )
    generate_rnd_marker_effects_parser$add_argument(arg$rate_dom$flag,
      help = arg$rate_dom$help,
      type = arg$rate_dom$type,
      default = arg$rate_dom$default,
      required = FALSE
    )
    generate_rnd_marker_effects_parser$add_argument(arg$generate_dominance_effects$flag,
      help = arg$generate_dominance_effects$help,
      action = "store_true"
    )
    generate_rnd_marker_effects_parser$add_argument(arg$prop_add$flag,
      help = arg$prop_add$help,
      type = arg$prop_add$type,
      default = arg$prop_add$default,
      required = FALSE
    )
    generate_rnd_marker_effects_parser$add_argument(arg$prop_dom$flag,
      help = arg$prop_dom$help,
      type = arg$prop_dom$type,
      default = arg$prop_dom$default,
      required = FALSE
    )
    generate_rnd_marker_effects_parser$add_argument(arg$no_dual_effects$flag,
      help = arg$no_dual_effects$help,
      action = "store_true"
    )
    generate_rnd_marker_effects_parser$add_argument(arg$rng_seed$flag,
      help = arg$rng_seed$help,
      type = arg$rng_seed$type,
      default = arg$rng_seed$default,
      required = FALSE
    )
    generate_rnd_marker_effects_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )

    ## simulate_phenotype ----
    simulate_phenotype_parser <- main_subparsers$add_parser(
      "simulate-phenotype",
      help = "Simulate Phenotypic Values from Genotypic Data and Marker Effects.",
      description = "Simulates phenotypic values for individuals based on their genotypic data and marker effects. The generated phenotypic values are saved as a CSV file. The command print the distribution of the phenotypic values for the given genotype file.",
      formatter_class = formatter_class
    )
    simulate_phenotype_parser$add_argument(arg$genoFile$flag,
      help = arg$genoFile$help,
      type = arg$genoFile$type,
      required = TRUE
    )
    simulate_phenotype_parser$add_argument(arg$markersEffectsFile_2$flag,
      help = arg$markersEffectsFile_2$help,
      type = arg$markersEffectsFile_2$type,
      required = TRUE
    )
    simulate_phenotype_parser$add_argument(arg$outFile$flag,
      help = arg$outFile$help,
      type = arg$outFile$type,
      required = TRUE
    )
    simulate_phenotype_parser$add_argument(arg$mean$flag,
      help = arg$mean$help,
      type = arg$mean$type,
      default = arg$mean$default,
      required = FALSE
    )
    simulate_phenotype_parser$add_argument(arg$heritability$flag,
      help = arg$heritability$help,
      type = arg$heritability$type,
      default = arg$heritability$default,
      required = FALSE
    )
    simulate_phenotype_parser$add_argument(arg$sd_noise$flag,
      help = arg$sd_noise$help,
      type = arg$sd_noise$type,
      default = arg$sd_noise$default,
      required = FALSE
    )
    simulate_phenotype_parser$add_argument(arg$trait_name$flag,
      help = arg$trait_name$help,
      type = arg$trait_name$type,
      default = arg$trait_name$default,
      required = FALSE
    )
    simulate_phenotype_parser$add_argument(arg$rng_seed$flag,
      help = arg$rng_seed$help,
      type = arg$rng_seed$type,
      default = arg$rng_seed$default,
      required = FALSE
    )
    simulate_phenotype_parser$add_argument(arg$json_error$flag,
      help = arg$json_error$help,
      default = arg$json_error$default,
      action = "store_true",
      required = FALSE
    )





    # Parse the command line arguments
    args <- main_parser$parse_args()

    # show help if no command is provided
    if (is.null(args$command)) {
      main_parser$print_help()
      quit(save = "no", status = 0)
    }

    ERROR_AS_JSON <- args$json_errors

    # log input parameters
    time <- as.character(Sys.time())
    cat(time, "call to R-geno-engine.R:\n", paste0(names(args), " = ", args, "\n"))


    # load engine's R packages dependencies:
    suppressPackageStartupMessages({
      library(R6) # R object oriented
    })

    # source R functions:
    invisible(
      sapply(
        FUN = source,
        X = list.files(file.path(projectRoot, "src"), pattern = ".R$", full.names = T)
      )
    )

    if (args$command == "gwas") {
      # gwas ----
      gwas_results <- run_gwas(
        genoFile = args$genoFile,
        phenoFile = args$phenoFile,
        trait = args$trait,
        test = args$test,
        fixed = args$fixed,
        response = args$response,
        thresh_maf = args$thresh_maf,
        thresh_callrate = args$thresh_callrate,
        n_markers = args$n_markers,
        n_markers_tolerance = args$n_markers_tolerance,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "gwas-manplot") {
      # gwas-manplot ----
      if (identical(args$chr, "")) {
        args$chr <- NA
      }
      p <- draw_manhattanPlot(
        gwasFile = args$gwasFile,
        adj_method = args$adj_method,
        thresh_p = args$thresh_p,
        chr = args$chr,
        interactive = !args$no_interactive,
        filter_pAdj = args$filter_pAdj,
        filter_nPoints = args$filter_nPoints,
        filter_quant = args$filter_quant,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "gwas-adjresults") {
      # gwas-adjresults ----
      gwas_adj <- run_resAdjustment(
        gwasFile = args$gwasFile,
        gwasUrl = args$gwasUrl,
        adj_method = args$adj_method,
        filter_pAdj = args$filter_pAdj,
        filter_nPoints = args$filter_nPoints,
        filter_quant = args$filter_quant,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "ldplot") {
      # ldplot ----
      imgFile <- draw_ldPlot(
        genoFile = args$genoFile,
        from = args$from,
        to = args$to,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "relmat-ped") {
      # relmat-ped ----
      relMat <- calc_pedRelMat(
        pedFile = args$pedFile,
        unknown_string = args$unknown_string,
        header = !args$no_header,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "relmat-geno") {
      # relmat-geno ----
      relMat <- calc_genoRelMat(
        genoFile = args$genoFile,
        n_markers = args$n_markers,
        n_markers_tolerance = args$n_markers_tolerance,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "relmat-combined") {
      # relmat-combined ----
      relMat <- calc_combinedRelMat(
        pedRelMatFile = args$ped_relmatFile,
        genoRelMatFile = args$geno_relmatFile,
        method = args$combine_method,
        tau = args$tau,
        omega = args$omega,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "relmat-heatmap") {
      # relmat-heatmap ----
      p <- draw_relHeatmap(
        relMatFile = args$relmatFile,
        interactive = !args$no_interactive,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "pedNetwork") {
      # pedNetwork ----
      p <- draw_pedNetwork(
        pedFile = args$pedFile,
        unknown_string = args$unknown_string,
        header = !args$no_header,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "crossing-simulation") {
      # crossing-simulation ----
      outFile <- crossingSimulation(
        genoFile = args$genoFile,
        crossTableFile = args$crossTableFile,
        SNPcoordFile = args$SNPcoordFile,
        nCross = args$nCross,
        outFile = args$outFile
      )
      cat(outFile)
      quit(save = "no", status = 0)
    } else if (args$command == "progeny-blup-calculation") {
      # progeny-blup-calculation ----
      progBlups <- calc_progenyBlupEstimation(
        genoFile = args$genoFile,
        crossTableFile = args$crossTableFile,
        SNPcoordFile = args$SNPcoordFile,
        markerEffectsFile = args$markerEffectsFile,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "progeny-blup-plot") {
      # progeny-blup-plot ----
      if (identical(args$trait, "")) {
        args$trait <- NULL
      }
      p <- draw_progBlupsPlot(
        progEstimFile = args$progeniesBlupFile,
        sorting = args$sorting,
        errorBarInterval = args$error_bar_interval,
        y_axisName = args$y_axis_name,
        trait = args$trait,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "progeny-blup-plot-2-traits") {
      # progeny-blup-plot-2-traits ----
      if (identical(args$trait, "")) {
        args$trait <- NULL
      }
      p <- draw_progBlupsPlot_2traits(
        progEstimFile = args$progeniesBlupFile,
        x_trait = args$x_trait,
        y_trait = args$y_trait,
        confidenceLevel = args$confidence_level,
        x_suffix = args$x_suffix,
        y_suffix = args$y_suffix,
        ellipses_npoints = args$ellipses_npoints,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "evaluate-gs-model") {
      # evaluate-gs-model ----
      evaluation <- cross_validation_evaluation_main(
        genoFile = args$geno,
        phenoFile = args$pheno,
        trait = args$trait,
        with_dominance = args$with_dominance,
        n_folds = args$n_folds,
        n_repetitions = args$n_repetitions,
        thresh_maf = 0,
        n_markers = args$n_markers,
        n_markers_tolerance = args$n_markers_tolerance,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "plot-gs-model-evaluation") {
      # plot-gs-model-evaluation ----
      plot <- draw_evaluation_plot(
        evaluationFile = args$evaluation_file,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "train-gs-model") {
      # train-gs-model ----
      model <- train_gs_model_main(
        genoFile = args$geno,
        phenoFile = args$pheno,
        trait = args$trait,
        with_dominance = args$with_dominance,
        thresh_maf = 0,
        n_markers = args$n_markers,
        n_markers_tolerance = args$n_markers_tolerance,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "gs-predictions") {
      # gs-predictions ----
      predictions <- predict_gs_model_main(
        genoFile = args$geno,
        markerEffectsFile = args$markerEffectsFile,
        outFile = args$outFile
      )
      quit(save = "no", status = 0)
    } else if (args$command == "generate-rnd-marker-effects") {
      generate_rnd_marker_effects(
        genoFile = args$genoFile,
        outFile = args$outFile,
        rate_add = args$rate_add,
        rate_dom = args$rate_dom,
        dominance = args$dominance,
        prop_snp_with_effect_add = args$prop_add,
        prop_snp_with_effect_dom = args$prop_dom,
        allow_dual_effects = !args$no_dual_effects,
        rnd_seed = args$rng_seed
      )
      quit(save = "no", status = 0)
    } else if (args$command == "simulate-phenotype") {
      simulate_phenotype(
        genoFile = args$genoFile,
        markerEffectsFile = args$markerEffectsFile,
        outFile = args$outFile,
        mean = args$mean,
        heritability = args$heritability,
        sd_noise = args$sd_noise,
        trait_name = args$trait_name,
        rnd_seed = args$rng_seed
      )
      quit(save = "no", status = 0)
    }



    # Should not arrive here.
    quit(save = "no", status = 1)
  }
)
