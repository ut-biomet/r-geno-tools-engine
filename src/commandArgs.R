# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# This file contain the descriptions of all the command line tool's arguments


arg <- list()

# genoFile ----
arg$genoFile$flag <- '--genoFile'
arg$genoFile$help <- 'Path of the geno data file (`.vcf` or `.vcf.gz` file)'
arg$genoFile$type <- 'character'


# phasedGenoFile ----
arg$phasedGenoFile$flag <- '--genoFile'
arg$phasedGenoFile$help <- paste('Path of the phased genetic data file',
                                 '(`.vcf` or `.vcf.gz` file)')
arg$phasedGenoFile$type <- 'character'


# phenoFile ----
arg$phenoFile$flag = '--phenoFile'
arg$phenoFile$help = paste('Path of the phenotypic data file (`csv` file).',
                           'The individuals names should be on the first',
                           'columns, and no duplicate are allowed')
arg$phenoFile$type = 'character'


# pedFile ----
arg$pedFile$flag = '--pedFile'
arg$pedFile$help = paste('Path of the pedigree data file (`csv` file).',
                         'It should have 3 columns with the individual names',
                         'on the 1st, and the parents name on the 2nd and 3rd.',
                         '. No duplicate are allowed')
arg$pedFile$type = 'character'


# gwasFile ----
arg$gwasFile$flag = '--gwasFile'
arg$gwasFile$help = 'Path of the gwas result data file (json file)'
arg$gwasFile$type = 'character'


# relmatFile ----
arg$relmatFile$flag = '--relmatFile'
arg$relmatFile$help = 'path of a file relationship file generated by this engine.'
arg$relmatFile$type = 'character'


# ped-relmatFile ----
arg$ped_relmatFile$flag = '--ped-relmatFile'
arg$ped_relmatFile$help = 'path of a pedigree relationship matrix generated by the the engine.'
arg$ped_relmatFile$type = 'character'


# geno-relmatFile ----
arg$geno_relmatFile$flag = '--geno-relmatFile'
arg$geno_relmatFile$help = 'path of a genomic relationship matrix generated by the the engine.'
arg$geno_relmatFile$type = 'character'


# outFile ----
arg$outFile$flag = '--outFile'
arg$outFile$help = 'Path of the file where to save the results'
arg$outFile$type = 'character'


# adj-method ----
arg$adj_method$flag = '--adj-method'
arg$adj_method$default = 'bonferroni'
arg$adj_method$help = paste0('p-value correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none".\\nDefault: ', arg$adj_method$default, '. (see https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust for more details)')
arg$adj_method$type = 'character'


# combine-method ----
arg$combine_method$flag = '--combine-method'
arg$combine_method$default = 'Legarra'
arg$combine_method$help = paste0('method to use, either "Legarra" or "Martini"',
                                 ' The default is:', arg$combine_method$default)
arg$combine_method$type = 'character'


# thresh-p ----
arg$thresh_p$flag = '--thresh-p'
arg$thresh_p$default = 0.05
arg$thresh_p$help = paste0('p-value significant threshold (default: ',
                           arg$thresh_p$default, ')')
arg$thresh_p$type = 'double'


# filter-pAdj ----
arg$filter_pAdj$flag = '--filter-pAdj'
arg$filter_pAdj$default = 1
arg$filter_pAdj$help = 'Threshold to remove points with pAdj < filter-pAdj from the plot (default no filtering)'
arg$filter_pAdj$type = 'double'


# filter_nPoints ----
arg$filter_nPoints$flag = '--filter-nPoints'
arg$filter_nPoints$default = 1e100
arg$filter_nPoints$help = 'Threshold to keep only the filter-nPoints points with the lowest p-values for the plot (default no filtering)'
arg$filter_nPoints$type = 'double'


# filter_quant ----
arg$filter_quant$flag = '--filter-quant'
arg$filter_quant$default = 1
arg$filter_quant$help = 'Threshold to keep only the filter-quant*100 percent of the points with the lowest p-values for the plot (default no filtering)'
arg$filter_quant$type = 'double'


# test ----
arg$test$flag = '--test'
arg$test$default = 'score'
arg$test$help = paste0('Which test to use. Either `"score"`, `"wald"` or `"lrt"`.\\nFor binary phenotypes, test = `"score"` is mandatory. Default ', arg$test$default, '.\\nFor more information about this parameters see: https://www.rdocumentation.org/packages/gaston/versions/1.4.9/topics/association.test')
arg$test$type = 'character'


# fixed ----
arg$fixed$flag = '--fixed'
arg$fixed$default = 0
arg$fixed$help = paste0('Number of Principal Components to include in the model with fixed effect (for test = `"wald"` or `"lrt"`). Default value is ', arg$fixed$default, '.\\nFor more information about this parameters see: https://www.rdocumentation.org/packages/gaston/versions/1.4.9/topics/association.test')
arg$fixed$type = 'integer'


# response ----
arg$response$flag = '--response'
arg$response$default = 'quantitative'
arg$response$help = paste0('Either "quantitative" or "binary". Is the trait a quantitative or a binary phenotype? Default value is ', arg$response$default)
arg$response$type = 'character'


# trait ----
arg$trait$flag = '--trait'
arg$trait$help = 'Name of the trait to analyze. (Must be a column name of the phenotypic file)'
arg$trait$type = 'character'


# thresh_maf ----
arg$thresh_maf$flag = '--thresh-maf'
arg$thresh_maf$default = 0
arg$thresh_maf$help = 'Threshold for filtering markers. Only markers with minor allele frequency > `thresh-maf` will be kept. Default: no filtering.'
arg$thresh_maf$type = 'double'


# thresh_callrate ----
arg$thresh_callrate$flag = '--thresh-callrate'
arg$thresh_callrate$default = 0
arg$thresh_callrate$help = 'Threshold for filtering markers. Only markers with a callrate > `thresh-callrate` will be kept. Default no filtering.'
arg$thresh_callrate$type = 'double'


# chr ----
arg$chr$flag = '--chr'
arg$chr$default = ""
arg$chr$help = 'Name of the chromosome to show (show all by default)'
arg$chr$type = 'character'


# no-interactive ----
arg$no_interactive$flag = '--no-interactive'
arg$no_interactive$default = FALSE
arg$no_interactive$help = 'Make the output a static png image (if not the output will be an interactive html plot)'


# from ----
arg$from$flag = '--from'
arg$from$help = 'Lower bound of the range of SNPs for which the LD is computed (`from` must be lower than `to`)'
arg$from$type = 'integer'


# to ----
arg$to$flag = '--to'
arg$to$help = 'Upper bound of the range of SNPs for which the LD is computed (the total number of SNP should be lower than 50)'
arg$to$type = 'integer'


# unknown_string ----
arg$u_string$flag = '--unknown-string'
arg$u_string$default = ''
arg$u_string$help = 'A character vector of strings which should be interpreted as "unknown parent". By default: missing value in the file.'
arg$u_string$type = 'character'


# no-header ----
arg$no_header$flag = '--no-header'
arg$no_header$default = FALSE
arg$no_header$help = paste('To specify the input file do not have a header. In any cases,',
                           'the column 1 will be interpreted as the individual id',
                           ', column 2 as the first parent, column 3 as the',
                           'second parent.')


# tau ----
arg$tau$flag = '--tau'
arg$tau$default = NULL
arg$tau$help = "tau parameter of the Martini's method"
arg$tau$type = 'double'


# omega ----
arg$omega$flag = '--omega'
arg$omega$default = NULL
arg$omega$help = "omega parameter of the Martini's method"
arg$omega$type = 'double'


# crossTableFile ----
arg$crossTableFile$flag = '--crossTableFile'
arg$crossTableFile$help = paste("path of the crossing table data file",
                                "(`csv` file of 2 or 3 columns).",
                                "It must contain the names of the variables",
                                "as its first line. The column 1 and 2 will",
                                "be interpreted as the parents ids.",
                                "The optional third column will be interpreted",
                                "as the offspring base name.")
arg$crossTableFile$type = 'character'


# SNPcoordFile ----
arg$SNPcoordFile$flag = '--SNPcoordFile'
arg$SNPcoordFile$help = paste("path of the SNPs coordinates",
                              "file (`csv` file). This `.csv` file can have",
                              "4 named columns:\\n",
                              "- `chr`: chromosome name holding the SNP\\n",
                              "- `physPos`: physical position of the SNP on",
                              "the chromosome\\n",
                              "- `linkMapPos`: linkage map position of the SNP",
                              "on the chromosome in Morgan (mandatory)\\n",
                              "- `SNPid`: ID of the SNP\\n",
                              "Column `physPos` is optional except in some particular case (see below).\\n",
                              "If this column is provided (or contain only missing values), it should",
                              "exactly match the physical positions of the SNP specified in the VCF file.\\n",
                              "If `SNPid` columns is missing or have missing values, the SNPid will be",
                              "automatically imputed using the convention `chr@physPos` therefore columns",
                              "`chr` and `physPos` should not have any missing values in this case.")
arg$SNPcoordFile$type = 'character'


# markersEffects ----
arg$markersEffectsFile$flag = '--markerEffectsFile'
arg$markersEffectsFile$help = paste(
  "path of the marker effects file (`csv`, or `json` file).\\n",
  "For `.csv`, the file file should have 2 named columns:\\n",
  "  - `SNPid`: Marker id\\n",
  "  - `effects`: effect of the corresponding marker\\n",
  "One can specify the intercept using \"--INTERCEPT--\" as SNPid.\\n",
  "For `.json`, the file should have 2 Key-value pairs:\\n",
  "  - `intercept`: a number with the value of the intercept.\\n",
  "  - `coefficient`: a nested object with SNPids as keys and their",
  "corresponding effects as values.\\n",
  "For example :\\n",
  ' {\\n',
  '  "intercept": 100,\\n',
  '  "coefficients": {\\n',
  '     "SNP01": 1.02e-06,\\n',
  '     "SNP02": 0.42,\\n',
  '    "SNP03": 0.0\\n',
  '  }\\n',
  ' })\\n')
arg$markersEffectsFile$type = 'character'


# nCross ----
arg$nCross$flag = '--nCross'
arg$nCross$default = 15
arg$nCross$help = paste("[Optional] Number of cross to simulate for each parent",
                        "pair defined in the crossing table.",
                        "Default is :", arg$nCross$default)
arg$nCross$type = 'integer'


# outGenoFile ----
arg$outGenoFile$flag = '--outFile'
arg$outGenoFile$help = paste("path of the `.vcf.gz` file containing the",
                             "simulated genotypes of the offspring. It must",
                             "end by `.vcf.gz`.")
arg$outGenoFile$type = 'character'


# outProgBlupFile ----
arg$outProgBlupFile$flag = '--outFile'
arg$outProgBlupFile$help = paste(
  "`.json` file path where to save the data. If the file already exists,",
  "it will be overwritten."
)
arg$outProgBlupFile$type = 'character'


# progeniesBlupFile ----
arg$progeniesBlupFile$flag = '--progeniesBlupFile'
arg$progeniesBlupFile$help = paste(
  "path of the progeny BLUP estimation file generated by",
  "r-geno-tools-engine containing the blup estimations of the progenies of",
  "some crosses (`json` file)."
  )
arg$progeniesBlupFile$type = 'character'

# sorting ----
arg$sorting$flag = '--sorting'
arg$sorting$default = "alpha"
arg$sorting$help = paste(
  "method to sort the individuals (X axis) can be:\\n",
  '- "asc": sort the BLUP expected value in ascending order',
  "(from left to right)\\n",
  '- "dec": sort the BLUP expected value in decreasing order',
  "(from left to right)\\n",
  "- any other value will sort the individuals in alphabetical order",
  "(from left to right)\\n",
  "Default is :", arg$sorting$default)
arg$sorting$type = 'character'

# errorBarInterval ----
arg$errorBarInterval$flag = '--error-bar-interval'
arg$errorBarInterval$default = 0.95
arg$errorBarInterval$help = paste(
  "Length of interval of interest represented by the error bars.",
  "Values between 0 and 1, eg. 0.95 give the 95 percent interval.",
  "Default is :", arg$errorBarInterval$default)
arg$errorBarInterval$type = 'double'

# y_axisName ----
arg$y_axisName$flag = '--y-axis-name'
arg$y_axisName$default = "Genetic values"
arg$y_axisName$help = paste(
  "Name of the Y axis of the plot.",
  "Default is :", arg$y_axisName$default)
arg$y_axisName$type = 'character'

# x_trait ----
arg$x_trait$flag = '--x-trait'
arg$x_trait$help = paste("Name of the trait to show on the X axis of the plot.")
arg$x_trait$type = 'character'


# y_trait ----
arg$y_trait$flag = '--y-trait'
arg$y_trait$help = paste("Name of the trait to show on the Y axis of the plot.")
arg$y_trait$type = 'character'

# confidenceLevel ----
arg$confidenceLevel$flag = '--confidence-level'
arg$confidenceLevel$default = 0.95
arg$confidenceLevel$help = paste(
  "Size of the region represented by the ellipses.",
  "Values between 0 and 1, eg. 0.95 give the 95 percent interval.",
  "Default is :", arg$confidenceLevel$default)
arg$confidenceLevel$type = 'double'

# x_suffix ----
arg$x_suffix$flag = '--x-suffix'
arg$x_suffix$default = ""
arg$x_suffix$help = paste(
  "Optional suffix to add to the name of the X axis of the plot.",
  "Default is :", arg$x_suffix$default)
arg$x_suffix$type = 'character'

# y_suffix ----
arg$y_suffix$flag = '--y-suffix'
arg$y_suffix$default = ""
arg$y_suffix$help = paste(
  "Optional suffix to add to the name of the Y axis of the plot.",
  "Default is :", arg$y_suffix$default)
arg$y_suffix$type = 'character'

# ellipses_npoints ----
arg$ellipses_npoints$flag = '--ellipses-npoints'
arg$ellipses_npoints$default = 100
arg$ellipses_npoints$help = paste(
  "Number of points to use to draw the ellipses.",
  "Default is :", arg$ellipses_npoints$default)
arg$ellipses_npoints$type = 'integer'


# json-errors ----
arg$json_error$flag = '--json-errors'
arg$json_error$default = FALSE
arg$json_error$help = 'Write errors as json in stderr (default FALSE)'

# with-dominance ----
arg$with_dominance$flag = '--with-dominance'
arg$with_dominance$help = 'Should the model include dominance effects'
arg$with_dominance$type = 'logical'

# n-folds ----
arg$n_folds$flag = '--n-folds'
arg$n_folds$default = 10
arg$n_folds$help = paste0('Number of folds for each cross-validation (default ',
                          arg$n_folds$default,
                          ')')
arg$n_folds$type = 'integer'

# n-folds ----
arg$n_repetitions$flag = '--n-repetitions'
arg$n_repetitions$default = 5
arg$n_repetitions$help = paste0('Number of cross-validation repetition (default ',
                          arg$n_repetitions$default,
                          ')')
arg$n_repetitions$type = 'integer'

# evaluation-file ----
arg$evaluation_file$flag <- '--evaluation-file'
arg$evaluation_file$help <- 'Path of the model evaluation results'
arg$evaluation_file$type <- 'character'
