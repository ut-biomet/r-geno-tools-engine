#! /usr/bin/Rscript

# Author: Julien Diot juliendiot@ut-biomet.org
# Tue 24 Aug 2021 The University of Tokyo
#
# Description:
# Main script for using GWAS-Engine

# Create a parser
suppressPackageStartupMessages(library(argparser, quietly=TRUE))


p <- arg_parser(paste0(
  "Run GWAS-Engine's tools"))

# Add command line arguments
p <- add_argument(p, "fun", help="GWAS-Engine's function you want to run: `gwas`, `manplot`, `ldplot` or `adjresults`", type = "character")

p <- add_argument(p, "--genoFile", help="[`gwas`, `ldplot`] path of the geno data file (`.vcf` or `.vcf.gz` file)", type = "character")
p <- add_argument(p, "--phenoFile", help="[`gwas`] path of the phenotypic data file (`csv` file)", type = "character")
p <- add_argument(p, "--gwasFile", help="[`manplot`, `adjresults`] path of the gwas result data file (json file)", type = "character")
p <- add_argument(p, "--outFile", help="[`gwas`, `manplot`, `ldplot`, `adjresults`]  file where to save the results,", type = "character")


p <- add_argument(p, "--adj_method", help='[`manplot`, `adjresults`] p-value correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none". Default: "bonferroni". (see https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust for more details)', type = "character", default = "bonferroni")
p <- add_argument(p, "--thresh_p", help="[`manplot`, `adjresults`] p value significant threshold (default 0.05)", type = "numeric", default = 0.05)


p <- add_argument(p, "--filter_pAdj", help="[`manplot`] threshold to remove points with pAdj < filter_pAdj from the plot (default no filtering)", type = "numeric", default = 1)
p <- add_argument(p, "--filter_nPoints", help="[`manplot`] threshold to keep only the filter_nPoints with the lowest p-values for the plot (default no filtering)", type = "numeric", default = Inf)
p <- add_argument(p, "--filter_quant", help="[`manplot`] threshold to keep only the filter_quant*100 % of the points with the lowest p-values for the plot (default no filtering)", type = "numeric", default = 1)


# gwas specific parameters:
p <- add_argument(p, "--test", help='[`gwas`] Which test to use. Either `"score"`, `"wald"` or `"lrt"`. For binary phenotypes, test = `"score"` is mandatory. For more information about this parameters see: https://www.rdocumentation.org/packages/gaston/versions/1.4.9/topics/association.test' , type = "character")
p <- add_argument(p, "--fixed", help='[`gwas`] Number of Principal Components to include in the model with fixed effect (for test = `"wald"` or `"lrt"`). Default value is 0. For more information about this parameters see: https://www.rdocumentation.org/packages/gaston/versions/1.4.9/topics/association.test', type = "integer", default = 0)
p <- add_argument(p, "--response", help='[`gwas`] Either "quantitative" or "binary". Is the trait a quantitative or a binary phenotype? Default value is "quantitative"' , type = "character", default = "quantitative")
p <- add_argument(p, "--trait", help="[`gwas`] name of the trait to analyze. Must be a column name of the phenotypic file.", type = "character")

p <- add_argument(p, "--thresh_maf", help="[`gwas`] Threshold for filtering markers. Only markers with minor allele frequency > `thresh_maf` will be kept.", type = "numeric")
p <- add_argument(p, "--thresh_callrate", help="[`gwas`] Threshold for filtering markers. Only markers with a callrate > `thresh_callrate` will be kept.", type = "numeric")


# manplot specific parameters:
p <- add_argument(p, "--chr", help="[`manplot`] name of the chromosome to show (show all by default)", type = "character", default = NA)
p <- add_argument(p, "--interactive", help="[`manplot`] should the plot be interactive: TRUE or FALSE (the default is TRUE)", type = "logical", default = TRUE)


# ldplot specific parameters:
p <- add_argument(p, "--from", help="[`ldplot`] lower bound of the range of SNPs for which the LD is computed (`from` must be lower than `to`)", type = "integer")
p <- add_argument(p, "--to", help="[`ldplot`] upper bound of the range of SNPs for which the LD is computed (the total number of SNP should be lower than 50)", type = "integer")


# Parse the command line arguments
args <- parse_args(p)


# log input parameters
time <- as.character(Sys.time())
inArgs <- args[seq_len(length(args) - 1) + 1]
cat(time, 'call to GWAS-Engine.R:\n', paste0('\t', names(inArgs), ' = ', inArgs, '\n'))


if (!args$fun %in% c('gwas', 'manplot', 'ldplot', 'adjresults')) {
  stop("`fun` should be one of 'gwas', 'manplot', 'ldplot', 'adjresults'")
}

# load R packages:
suppressPackageStartupMessages({
  library(gaston) # for many functions
  library(jsonlite) # manage json format
  library(manhattanly) # manhattan plot using plotly
  library(R6) # R object oriented
})

# source R functions:
invisible(
  sapply(FUN = source,
         X = list.files("src", pattern = ".R$",full.names = T))
)


if (args$fun == "gwas") {
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

} else if (args$fun == "manplot") {
  p <- draw_manhattanPlot(gwasFile = args$gwasFile,
                          adj_method = args$adj_method,
                          thresh_p = args$thresh_p,
                          chr = args$chr,
                          interactive = args$interactive,
                          filter_pAdj = args$filter_pAdj,
                          filter_nPoints = args$filter_nPoints,
                          filter_quant = args$filter_quant,
                          outFile = args$outFile)
  quit(save = "no", status = 0)

} else if (args$fun == "adjresults") {
  gwas_adj <- run_resAdjustment(gwasFile = args$gwasFile,
                                gwasUrl = args$gwasUrl,
                                adj_method = args$adj_method,
                                filter_pAdj = args$filter_pAdj,
                                filter_nPoints = args$filter_nPoints,
                                filter_quant = args$filter_quant,
                                outFile = args$outFile)
  quit(save = "no", status = 0)

} else if (args$fun == "ldplot") {
  imgFile <- draw_ldPlot(genoFile = args$genoFile,
                         from = args$from,
                         to = args$to,
                         outFile = args$outFile)
  quit(save = "no", status = 0)
}

# Should not arrive here.
quit(save = "no", status = 1)
