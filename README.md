
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GWAS-Engine

<!-- badges: start -->
<!-- badges: end -->

The goal of GWAS-Engine is to provide simple “input/output” style code
for GWAS analysis in order to be used inside an API.

# How to use GWAS-Engine

## Dependencies

``` r
# required packages
library(gaston) # for many functions
#> Loading required package: Rcpp
#> Loading required package: RcppParallel
#> 
#> Attaching package: 'RcppParallel'
#> The following object is masked from 'package:Rcpp':
#> 
#>     LdFlags
#> Gaston set number of threads to 4. Use setThreadOptions() to modify this.
#> 
#> Attaching package: 'gaston'
#> The following object is masked from 'package:stats':
#> 
#>     sigma
#> The following objects are masked from 'package:base':
#> 
#>     cbind, rbind
library(rjson) # manage json format
library(manhattanly) # manhattan plot using plotly
#> See example usage at http://sahirbhatnagar.com/manhattanly/
```

## Load R scripts

Load all `.R` files in folder `src`.

``` r
invisible(
  sapply(FUN = source,
         X = list.files("src", pattern = ".R$",full.names = T))
)
```

## Load data

Example data are stored in the `data` folder.

``` r
gDta <- readGenoData("data/markers/testMarkerData01.vcf.gz")
#> 2021-06-05 14:55:17 - r-readGenoData(): Check file extention ... 
#> 2021-06-05 14:55:17 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-06-05 14:55:18 - r-readGenoData(): Read geno file DONE 
#> 2021-06-05 14:55:18 - r-readGenoData(): DONE, return output.
pDta <- readPhenoData("data/pheno/testPhenoData01.csv")
#> 2021-06-05 14:55:18 - r-readPhenoData(): Read phenotypic file ... 
#> 2021-06-05 14:55:18 - r-readPhenoData(): Read phenotypic file DONE 
#> 2021-06-05 14:55:18 - r-readPhenoData(): DONE, return output.
```

For online data, they can be loaded using:

``` r
gDta <- readGenoData("data/markers/testMarkerData01.vcf.gz")
pDta <- readPhenoData("data/pheno/testPhenoData01.csv")
```

## GWAS

To run a GWAS analysis you should first process the data using
`prepareData` and then call the `gwas` function.

``` r
data <- prepareData(gDta=gDta, pDta = pDta)
#> 2021-06-05 14:55:18 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set ...
#> 2021-06-05 14:55:19 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set DONE
#> 2021-06-05 14:55:19 - r-prepareData(): reorder matrix ...
#> 2021-06-05 14:55:19 - r-prepareData(): reorder matrix DONE
#> 2021-06-05 14:55:19 - r-prepareData(): remove monomorphic geno ...
#> 2021-06-05 14:55:19 - r-prepareData(): remove monomorphic geno DONE
#> 2021-06-05 14:55:19 - r-prepareData(): calculate genetic relatinoal matrix ...
#> 2021-06-05 14:55:20 - r-prepareData(): calculate genetic relatinoal matrix DONE
#> 2021-06-05 14:55:20 - r-prepareData(): DONE, return output.

gwasRes <- gwas(data = data,
                trait = "Panicle.length",
                test = "score",
                fixed = 0,
                response = "quantitative",
                thresh_maf = 0.05,
                thresh_callrate = 0.95)
#> 2021-06-05 14:55:20 - r-gwas(): Check inputs ...
#> 2021-06-05 14:55:20 - r-gwas(): Check inputs DONE
#> 2021-06-05 14:55:20 - r-gwas(): aggregate data in bed matrix ...
#> 2021-06-05 14:55:20 - r-gwas(): aggregate data DONE
#> 2021-06-05 14:55:20 - r-gwas(): remove samples with missing phenotypic values ...
#> 2021-06-05 14:55:20 - r-gwas(): remove samples with missing phenotypic values DONE
#> 2021-06-05 14:55:20 - r-gwas(): filter SNPs ...
#> 2021-06-05 14:55:20 - r-gwas(): filter SNPs DONE
#> 2021-06-05 14:55:20 - r-gwas(): fit model ...
#> [Iteration 1] theta = 6.11889 3.12119
#> [Iteration 1] log L = -545.19
#> [Iteration 1] AI-REML update
#> [Iteration 1] ||gradient|| = 4.44998
#> [Iteration 2] theta = 2.40545 5.61423
#> [Iteration 2] log L = -536.314
#> [Iteration 2] AI-REML update
#> [Iteration 2] ||gradient|| = 5.68172
#> [Iteration 3] theta = 2.61764  6.5878
#> [Iteration 3] log L = -534.365
#> [Iteration 3] AI-REML update
#> [Iteration 3] ||gradient|| = 0.53883
#> [Iteration 4] theta = 2.60704 6.83002
#> [Iteration 4] log L = -534.325
#> [Iteration 4] AI-REML update
#> [Iteration 4] ||gradient|| = 0.0362724
#> [Iteration 5] theta =  2.5918 6.86787
#> [Iteration 5] log L = -534.324
#> [Iteration 5] AI-REML update
#> [Iteration 5] ||gradient|| = 0.0110436
#> [Iteration 6] theta = 2.58775 6.87669
#> [Iteration 6] log L = -534.324
#> [Iteration 6] AI-REML update
#> [Iteration 6] ||gradient|| = 0.00279298
#> [Iteration 7] theta = 2.58674 6.87888
#> [Iteration 7] log L = -534.324
#> [Iteration 7] AI-REML update
#> [Iteration 7] ||gradient|| = 0.000699933
#> [Iteration 8] theta = 2.58648 6.87942
#> [Iteration 8] log L = -534.324
#> [Iteration 8] AI-REML update
#> [Iteration 8] ||gradient|| = 0.000175022
#> [Iteration 9] theta = 2.58642 6.87956
#> [Iteration 9] log L = -534.324
#> [Iteration 9] AI-REML update
#> [Iteration 9] ||gradient|| = 4.37413e-05
#> [Iteration 10] theta = 2.5864 6.8796
#> [Iteration 10] log L = -534.324
#> [Iteration 10] AI-REML update
#> [Iteration 10] ||gradient|| = 1.09303e-05
#> [Iteration 11] theta = 2.5864 6.8796
#> [Iteration 11] log L = -534.324
#> [Iteration 11] AI-REML update
#> [Iteration 11] ||gradient|| = 2.73122e-06
#> 2021-06-05 14:55:21 - r-gwas(): fit model DONE
#> 2021-06-05 14:55:21 - r-gwas(): DONE, return output.
head(gwasRes)
#>   chr   pos           id A1 A2    freqA2      score         p
#> 1   1  9563  SNP-1.8562.  A  T 0.1382682 1.45574550 0.2276083
#> 2   1 25922 SNP-1.24921.  C  T 0.1278090 0.36309704 0.5467912
#> 3   1 26254 SNP-1.25253.  A  T 0.2955182 0.02865046 0.8655886
#> 4   1 30214 SNP-1.29213.  T  A 0.1344538 0.51074493 0.4748166
#> 5   1 31478 SNP-1.30477.  C  T 0.2262570 0.21535471 0.6426024
#> 6   1 32733 SNP-1.31732.  T  G 0.3011204 0.05018145 0.8227478
```

## Manhattan plot

``` r
#* Manhattan plot (type 2)
#* @tag Plots
#* @param modelS3Path url of the model data file (rds file)
#* @param adj_method either bonferroni or FDR
#* @param thresh_p
#* @param chr names of the chromosomes to show separated using comma. Show all chromosomes if nothing is specified.
#* @serializer htmlwidget
#* @get /manplot
# function(res, modelS3Path, adj_method, thresh_p = 0.05, chr = NA){
```

## LD PLot

``` r
#* LD plot
#* @tag Plots
#* @param geno_url url of the markers data file (.vcf.gz file)
#* @param from (total number of SNP should be < 50)
#* @param to (total number of SNP should be < 50)
#* @serializer png
#* @get /LDplot
# function(res, geno_url, from, to){
```
