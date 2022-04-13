#! /usr/local/bin/Rscript

# Author: Julien Diot juliendiot@ut-biomet.org
# Tue 24 Aug 2021 The University of Tokyo
#
# Description:
# Main script of R-geno-tools-enigne


library(argparse)

# Create parser ----
formatter_class <- "lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=42, width=100)"
main_parser = ArgumentParser(
  prog = './r-geno-tools-engine.R',
  description = "Toolbox for some genomic analysis written in R.",
  formatter_class = formatter_class,
  # argument_default = "True",
  usage = "%(prog)s"
  )
main_subparsers = main_parser$add_subparsers(prog = './r-geno-tools-engine.R',
                                             title = "Available commands",
                                             dest = "command")

# load argument descriptions
source('src/commandArgs.R')
parserList <- c() # to keep names of subcommands

## GWAS gwas ----
gwas_parser = main_subparsers$add_parser("gwas",
                                         help = "Do a gwas analysis",
                                         description = "Run a GWAS analysis",
                                         formatter_class = formatter_class,
                                         argument_default = "True")
gwas_parser$add_argument(arg$genoFile$flag,
                         help = arg$genoFile$help,
                         type = arg$genoFile$type,
                         required = TRUE)
gwas_parser$add_argument(arg$phenoFile$flag,
                         help = arg$phenoFile$help,
                         type = arg$phenoFile$type,
                         required = TRUE)
gwas_parser$add_argument(arg$trait$flag,
                         help = arg$trait$help,
                         type = arg$trait$type,
                         required = TRUE)
gwas_parser$add_argument(arg$outFile$flag,
                         help = arg$outFile$help,
                         type = arg$outFile$type,
                         required = TRUE)
gwas_parser$add_argument(arg$test$flag,
                         help = arg$test$help,
                         type = arg$test$type,
                         default = arg$test$default)
gwas_parser$add_argument(arg$fixed$flag,
                         help = arg$fixed$help,
                         type = arg$fixed$type,
                         default = arg$fixed$default)
gwas_parser$add_argument(arg$response$flag,
                         help = arg$response$help,
                         type = arg$response$type,
                         default = arg$response$default)
gwas_parser$add_argument(arg$thresh_maf$flag,
                         help = arg$thresh_maf$help,
                         type = arg$thresh_maf$type,
                         default = arg$thresh_maf$default)
gwas_parser$add_argument(arg$thresh_callrate$flag,
                         help = arg$thresh_callrate$help,
                         type = arg$thresh_callrate$type,
                         default = arg$thresh_callrate$default)


## GWAS manplot ----
manplot_parser = main_subparsers$add_parser("gwas-manplot",
                                            help = "Draw a Manhattan Plot",
                                            description = "Draw a Manhattan Plot",
                                            formatter_class=formatter_class)
manplot_parser$add_argument(arg$gwasFile$flag,
                            help = arg$gwasFile$help,
                            type = arg$gwasFile$type,
                         required = TRUE)
manplot_parser$add_argument(arg$outFile$flag,
                            help = arg$outFile$help,
                            type = arg$outFile$type,
                         required = TRUE)
manplot_parser$add_argument(arg$adj_method$flag,
                            help = arg$adj_method$help,
                            type = arg$adj_method$type,
                            default = arg$adj_method$default)
manplot_parser$add_argument(arg$thresh_p$flag,
                            help = arg$thresh_p$help,
                            type = arg$thresh_p$type,
                            default = arg$thresh_p$default)
manplot_parser$add_argument(arg$chr$flag,
                            help = arg$chr$help,
                            type = arg$chr$type,
                            default = arg$chr$default)
manplot_parser$add_argument(arg$no_interactive$flag,
                            help = arg$no_interactive$help,
                            action = 'store_true')
manplot_parser$add_argument(arg$filter_pAdj$flag,
                            help = arg$filter_pAdj$help,
                            type = arg$filter_pAdj$type,
                            default = arg$filter_pAdj$default)
manplot_parser$add_argument(arg$filter_nPoints$flag,
                            help = arg$filter_nPoints$help,
                            type = arg$filter_nPoints$type,
                            default = arg$filter_nPoints$default)
manplot_parser$add_argument(arg$filter_quant$flag,
                            help = arg$filter_quant$help,
                            type = arg$filter_quant$type,
                            default = arg$filter_quant$default)


## GWAS adjResults ----
adjResults_parser = main_subparsers$add_parser("gwas-adjresults",
                                               help = "Adjust GWAS p-values",
                                               description = "Adjust GWAS p-values and filter the results",
                                               formatter_class=formatter_class)
adjResults_parser$add_argument(arg$gwasFile$flag,
                               help = arg$gwasFile$help,
                               type = arg$gwasFile$type,
                               required = TRUE)
adjResults_parser$add_argument(arg$outFile$flag,
                               help = arg$outFile$help,
                               type = arg$outFile$type,
                               required = TRUE)
adjResults_parser$add_argument(arg$adj_method$flag,
                               help = arg$adj_method$help,
                               type = arg$adj_method$type,
                               default = arg$adj_method$default)
adjResults_parser$add_argument(arg$filter_pAdj$flag,
                               help = arg$filter_pAdj$help,
                               type = arg$filter_pAdj$type,
                               default = arg$filter_pAdj$default)
adjResults_parser$add_argument(arg$filter_nPoints$flag,
                               help = arg$filter_nPoints$help,
                               type = arg$filter_nPoints$type,
                               default = arg$filter_nPoints$default)
adjResults_parser$add_argument(arg$filter_quant$flag,
                               help = arg$filter_quant$help,
                               type = arg$filter_quant$type,
                               default = arg$filter_quant$default)




## LDplot ----
ldplot_parser = main_subparsers$add_parser("ldplot",
                                           help = "Get plot of the linkage disequilibrium between consecutive markers",
                                           description = "Get plot of the linkage disequilibrium between consecutive markers",
                                           formatter_class=formatter_class)
ldplot_parser$add_argument(arg$genoFile$flag,
                           help = arg$genoFile$help,
                           type = arg$genoFile$type,
                           required = TRUE)
ldplot_parser$add_argument(arg$outFile$flag,
                           help = arg$outFile$help,
                           type = arg$outFile$type,
                           required = TRUE)
ldplot_parser$add_argument(arg$from$flag,
                           help = arg$from$help,
                           type = arg$from$type,
                           required = TRUE)
ldplot_parser$add_argument(arg$to$flag,
                           help = arg$to$help,
                           type = arg$to$type,
                           required = TRUE)




## relmat pedigree ----
pedRelMat_parser = main_subparsers$add_parser("relmat-ped",
                                         help = "Pedigree relationship matrix",
                                         description = "Caclulate a pedigree relationship matrix.\\n\\nThe output file shoud end by either `.csv` or `.json`.",
                                         formatter_class=formatter_class)
pedRelMat_parser$add_argument(arg$pedFile$flag,
                              help = arg$pedFile$help,
                              type = arg$pedFile$type,
                              required = TRUE)
pedRelMat_parser$add_argument(arg$outFile$flag,
                              help = arg$outFile$help,
                              type = arg$outFile$type,
                              required = TRUE)
pedRelMat_parser$add_argument(arg$u_string$flag,
                              help = arg$u_string$help,
                              type = arg$u_string$type,
                              default = arg$u_string$default)
pedRelMat_parser$add_argument(arg$no_header$flag,
                              help = arg$no_header$help,
                              action = 'store_true')


## relmat genomic ----
genoRelMat_parser = main_subparsers$add_parser("relmat-geno",
                                         help = "Genomic relationship matrix",
                                         description = "Caclulate a genomic relationship matrix.\\n\\nThe output file shoud end by either `.csv` or `.json`.",
                                         formatter_class=formatter_class)
genoRelMat_parser$add_argument(arg$genoFile$flag,
                               help = arg$genoFile$help,
                               type = arg$genoFile$type,
                               required = TRUE)
genoRelMat_parser$add_argument(arg$outFile$flag,
                               help = arg$outFile$help,
                               type = arg$outFile$type,
                               required = TRUE)



## relmat combined ----
combinedRelMat_parser = main_subparsers$add_parser(
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
    "\\n \\n \\n "),
  formatter_class=formatter_class)
combinedRelMat_parser$add_argument(arg$ped_relmatFile$flag,
                                   help = arg$ped_relmatFile$help,
                                   type = arg$ped_relmatFile$type,
                                   required = TRUE)
combinedRelMat_parser$add_argument(arg$geno_relmatFile$flag,
                                   help = arg$geno_relmatFile$help,
                                   type = arg$geno_relmatFile$type,
                                   required = TRUE)
combinedRelMat_parser$add_argument(arg$combine_method$flag,
                                   help = arg$combine_method$help,
                                   type = arg$combine_method$type,
                                   required = TRUE)
combinedRelMat_parser$add_argument(arg$tau$flag,
                                   help = arg$tau$help,
                                   type = arg$tau$type)
combinedRelMat_parser$add_argument(arg$omega$flag,
                                   help = arg$omega$help,
                                   type = arg$omega$type)
combinedRelMat_parser$add_argument(arg$outFile$flag,
                                   help = arg$outFile$help,
                                   type = arg$outFile$type,
                                   required = TRUE)



## relmat heatmap ----
relMatHeat_parser = main_subparsers$add_parser("relmat-heatmap",
                                              help = "Draw a heatmap of a relationship matrix",
                                              description = "Draw a heatmap of a relationship matrix",
                                              formatter_class=formatter_class)
relMatHeat_parser$add_argument(arg$relmatFile$flag,
                               help = arg$relmatFile$help,
                               type = arg$relmatFile$type,
                               required = TRUE)
relMatHeat_parser$add_argument(arg$outFile$flag,
                               help = arg$outFile$help,
                               type = arg$outFile$type,
                               required = TRUE)
relMatHeat_parser$add_argument(arg$no_interactive$flag,
                               help = arg$no_interactive$help,
                               action = 'store_true')


## pedNetwork ----
pedNet_parser = main_subparsers$add_parser("pedNetwork",
                                           help = "Draw interactive pedigree network",
                                           description = "Draw interactive pedigree network",
                                           formatter_class=formatter_class)
pedNet_parser$add_argument(arg$pedFile$flag,
                           help = arg$pedFile$help,
                           type = arg$pedFile$type,
                           required = TRUE)
pedNet_parser$add_argument(arg$outFile$flag,
                           help = arg$outFile$help,
                           type = arg$outFile$type,
                           required = TRUE)
pedNet_parser$add_argument(arg$u_string$flag,
                           help = arg$u_string$help,
                           type = arg$u_string$type,
                           default = arg$u_string$default)
pedNet_parser$add_argument(arg$no_header$flag,
                           help = arg$no_header$help,
                           action = 'store_true')


# Parse the command line arguments
args = main_parser$parse_args()

# show help if no command is provided
if (is.null(args$command)) {
  main_parser$print_help()
  quit(save = "no", status = 0)
}

# log input parameters
time <- as.character(Sys.time())
cat(time, 'call to R-geno-engine.R:\n', paste0(names(args), ' = ', args, '\n'))


# load engine's R packages dependencies:
suppressPackageStartupMessages({
  library(R6) # R object oriented
})

# source R functions:
invisible(
  sapply(FUN = source,
         X = list.files("src", pattern = ".R$",full.names = T))
)

if (args$command == "gwas") {
  # gwas ----
  gwas_results <- run_gwas(genoFile = args$genoFile,
                           phenoFile = args$phenoFile,
                           trait = args$trait,
                           test = args$test,
                           fixed = args$fixed,
                           response = args$response,
                           thresh_maf = args$thresh_maf,
                           thresh_callrate = args$thresh_callrate,
                           outFile = args$outFile)
  quit(save = "no", status = 0)

} else if (args$command == "gwas-manplot") {
  # gwas-manplot ----
  if (identical(args$chr, '')) {
    args$chr <- NA
  }
  p <- draw_manhattanPlot(gwasFile = args$gwasFile,
                          adj_method = args$adj_method,
                          thresh_p = args$thresh_p,
                          chr = args$chr,
                          interactive = !args$no_interactive,
                          filter_pAdj = args$filter_pAdj,
                          filter_nPoints = args$filter_nPoints,
                          filter_quant = args$filter_quant,
                          outFile = args$outFile)
  quit(save = "no", status = 0)

} else if (args$command == "gwas-adjresults") {
  # gwas-adjresults ----
  gwas_adj <- run_resAdjustment(gwasFile = args$gwasFile,
                                gwasUrl = args$gwasUrl,
                                adj_method = args$adj_method,
                                filter_pAdj = args$filter_pAdj,
                                filter_nPoints = args$filter_nPoints,
                                filter_quant = args$filter_quant,
                                outFile = args$outFile)
  quit(save = "no", status = 0)

} else if (args$command == "ldplot") {
  # ldplot ----
  imgFile <- draw_ldPlot(genoFile = args$genoFile,
                         from = args$from,
                         to = args$to,
                         outFile = args$outFile)
  quit(save = "no", status = 0)

} else if (args$command == "relmat-ped") {
  # relmat-ped ----
  relMat <- calc_pedRelMAt(pedFile = args$pedFile,
                           unknown_string = args$unknown_string,
                           header = !args$no_header,
                           outFile = args$outFile)
  quit(save = 'no', status = 0)

} else if (args$command == "relmat-geno") {
  # relmat-geno ----
  relMat <- calc_genoRelMAt(genoFile = args$genoFile,
                            outFile = args$outFile)
  quit(save = 'no', status = 0)

} else if (args$command == "relmat-combined") {
  # relmat-combined ----
  relMat <- calc_combinedRelMat(pedRelMatFile = args$ped_relmatFile,
                                genoRelMatFile = args$geno_relmatFile,
                                method = args$combine_method,
                                tau = args$tau,
                                omega = args$omega,
                                outFile = args$outFile)
  quit(save = 'no', status = 0)

} else if (args$command == 'relmat-heatmap') {
  # relmat-heatmap ----
  p <- draw_relHeatmap(relMatFile = args$relmatFile,
                       interactive = !args$no_interactive,
                       outFile = args$outFile)
  quit(save = "no", status = 0)

} else if (args$command == 'pedNetwork') {
  # pedNetwork ----
  p <- draw_pedNetwork(pedFile = args$pedFile,
                       unknown_string = args$unknown_string,
                       header = !args$no_header,
                       outFile = args$outFile)
  quit(save = "no", status = 0)
}

# Should not arrive here.
quit(save = "no", status = 1)
