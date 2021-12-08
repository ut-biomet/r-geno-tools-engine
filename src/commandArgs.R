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


# phenoFile ----
arg$phenoFile$flag = '--phenoFile'
arg$phenoFile$help = 'Path of the phenotypic data file (`csv` file)'
arg$phenoFile$type = 'character'


# gwasFile ----
arg$gwasFile$flag = '--gwasFile'
arg$gwasFile$help = 'Path of the gwas result data file (json file)'
arg$gwasFile$type = 'character'


# outFile ----
arg$outFile$flag = '--outFile'
arg$outFile$help = 'Path of the file where to save the results'
arg$outFile$type = 'character'


# adj_method ----
arg$adj_method$flag = '--adj_method'
arg$adj_method$default = 'bonferroni'
arg$adj_method$help = paste0('p-value correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none". Default: ', arg$adj_method$default, '. (see https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust for more details)')
arg$adj_method$type = 'character'


# thresh_p ----
arg$thresh_p$flag = '--thresh_p'
arg$thresh_p$default = 0.05
arg$thresh_p$help = paste0('p-value significant threshold (default: ',
                           arg$thresh_p$default, ')')
arg$thresh_p$type = 'double'


# filter_pAdj ----
arg$filter_pAdj$flag = '--filter_pAdj'
arg$filter_pAdj$default = 1
arg$filter_pAdj$help = 'Threshold to remove points with pAdj < filter_pAdj from the plot (default no filtering)'
arg$filter_pAdj$type = 'double'


# filter_nPoints ----
arg$filter_nPoints$flag = '--filter_nPoints'
arg$filter_nPoints$default = "Inf"
arg$filter_nPoints$help = 'Threshold to keep only the filter_nPoints points with the lowest p-values for the plot (default no filtering)'
arg$filter_nPoints$type = 'double'


# filter_quant ----
arg$filter_quant$flag = '--filter_quant'
arg$filter_quant$default = 1
arg$filter_quant$help = 'Threshold to keep only the filter_quant*100 percent of the points with the lowest p-values for the plot (default no filtering)'
arg$filter_quant$type = 'double'


# test ----
arg$test$flag = '--test'
arg$test$default = 'score'
arg$test$help = paste0('Which test to use. Either `"score"`, `"wald"` or `"lrt"`. For binary phenotypes, test = `"score"` is mandatory. Default ', arg$test$default, '. For more information about this parameters see: https://www.rdocumentation.org/packages/gaston/versions/1.4.9/topics/association.test')
arg$test$type = 'character'


# fixed ----
arg$fixed$flag = '--fixed'
arg$fixed$default = 0
arg$fixed$help = paste0('Number of Principal Components to include in the model with fixed effect (for test = `"wald"` or `"lrt"`). Default value is ', arg$fixed$default, '. For more information about this parameters see: https://www.rdocumentation.org/packages/gaston/versions/1.4.9/topics/association.test')
arg$fixed$type = 'integer'


# response ----
arg$response$flag = '--response'
arg$response$default = 'quantitative'
arg$response$help = paste0('Either "quantitative" or "binary". Is the trait a quantitative or a binary phenotype? Default value is ', arg$response$default)
arg$response$type = 'character'


# trait ----
arg$trait$flag = '--trait'
arg$trait$help = 'Name of the trait to analyze. Must be a column name of the phenotypic file.'
arg$trait$type = 'character'


# thresh_maf ----
arg$thresh_maf$flag = '--thresh_maf'
arg$thresh_maf$default = 0
arg$thresh_maf$help = 'Threshold for filtering markers. Only markers with minor allele frequency > `thresh_maf` will be kept. Default: no filtering.'
arg$thresh_maf$type = 'double'


# thresh_callrate ----
arg$thresh_callrate$flag = '--thresh_callrate'
arg$thresh_callrate$default = 0
arg$thresh_callrate$help = 'Threshold for filtering markers. Only markers with a callrate > `thresh_callrate` will be kept. Default no filtering.'
arg$thresh_callrate$type = 'double'


# chr ----
arg$chr$flag = '--chr'
arg$chr$default = ""
arg$chr$help = 'Name of the chromosome to show (show all by default)'
arg$chr$type = 'character'


# interactive ----
arg$interactive$flag = '--interactive'
arg$interactive$default = TRUE
arg$interactive$help = paste0('Should the plot be interactive: TRUE or FALSE (the default is ', arg$interactive$default, ')')
arg$interactive$type = 'logical'


# from ----
arg$from$flag = '--from'
arg$from$help = 'Lower bound of the range of SNPs for which the LD is computed (`from` must be lower than `to`)'
arg$from$type = 'integer'


# to ----
arg$to$flag = '--to'
arg$to$help = 'Upper bound of the range of SNPs for which the LD is computed (the total number of SNP should be lower than 50)'
arg$to$type = 'integer'
