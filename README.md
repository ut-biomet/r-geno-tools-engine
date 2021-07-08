
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GWAS-Engine

<!-- badges: start -->
<!-- badges: end -->

The goal of GWAS-Engine is to provide simple “input/output” style code
for GWAS analysis in order to be used inside an API.

# How to use GWAS-Engine

## Dependencies

``` r
# required packages for using GWAS-Engine:
stopifnot(
  compareVersion("1.5.7", as.character(packageVersion("gaston"))) != 1
  )
stopifnot(
  compareVersion("1.7.2", as.character(packageVersion("jsonlite"))) != 1
  )
stopifnot(
  compareVersion("0.3.0", as.character(packageVersion("manhattanly"))) != 1
  )
stopifnot(
  compareVersion("2.5.0", as.character(packageVersion("R6"))) != 1
  )
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
library(jsonlite) # manage json format
library(manhattanly) # manhattan plot using plotly
#> See example usage at http://sahirbhatnagar.com/manhattanly/
library(R6) # R object oriented
```

## Load R scripts

Load all `.R` files in folder `src`.

``` r
invisible(
  sapply(FUN = source,
         X = list.files("src", pattern = ".R$",full.names = T))
)
```

## Main functions

There is the list of some main function that GWAS-engine can run.

### Run GWAS

``` r
gwas_results <- run_gwas(genoFile = "data/geno/testMarkerData01.vcf.gz",
                         phenoFile = "data/pheno/testPhenoData01.csv",
                         genoUrl = NULL,
                         phenoUrl = NULL,
                         trait = "Flowering.time.at.Arkansas",
                         test = "score",
                         fixed = 0,
                         response = "quantitative",
                         thresh_maf = 0.05,
                         thresh_callrate = 0.95,
                         dir = tempdir())
#> 2021-07-08 14:47:59 - r-run_gwas(): Get data ...
#> 2021-07-08 14:47:59 - r-readData(): get geno data ...
#> 2021-07-08 14:47:59 - r-readGenoData(): Check file extention ... 
#> 2021-07-08 14:47:59 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-07-08 14:48:01 - r-readGenoData(): Read geno file DONE 
#> 2021-07-08 14:48:01 - r-readGenoData(): DONE, return output.
#> 2021-07-08 14:48:01 - r-readData(): get geno data DONE
#> 2021-07-08 14:48:01 - r-readData(): get pheno data ...
#> 2021-07-08 14:48:01 - r-readPhenoData(): Read phenotypic file ... 
#> 2021-07-08 14:48:01 - r-readPhenoData(): Read phenotypic file DONE 
#> 2021-07-08 14:48:01 - r-readPhenoData(): DONE, return output.
#> 2021-07-08 14:48:01 - r-readData(): get pheno data DONE
#> 2021-07-08 14:48:01 - r-readData(): prepare data ...
#> 2021-07-08 14:48:01 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set ...
#> 2021-07-08 14:48:01 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set DONE
#> 2021-07-08 14:48:01 - r-prepareData(): reorder matrix ...
#> 2021-07-08 14:48:01 - r-prepareData(): reorder matrix DONE
#> 2021-07-08 14:48:01 - r-prepareData(): remove monomorphic geno ...
#> 2021-07-08 14:48:01 - r-prepareData(): remove monomorphic geno DONE
#> 2021-07-08 14:48:01 - r-prepareData(): calculate genetic relatinoal matrix ...
#> 2021-07-08 14:48:02 - r-prepareData(): calculate genetic relatinoal matrix DONE
#> 2021-07-08 14:48:02 - r-prepareData(): DONE, return output.
#> 2021-07-08 14:48:02 - r-readData(): prepare data DONE
#> 2021-07-08 14:48:02 - r-readData(): DONE, return output.
#> 2021-07-08 14:48:02 - r-run_gwas(): Get data DONE
#> 2021-07-08 14:48:02 - r-run_gwas(): GWAS analysis ...
#> 2021-07-08 14:48:02 - r-gwas(): Check inputs ...
#> 2021-07-08 14:48:02 - r-gwas(): Check inputs DONE
#> 2021-07-08 14:48:02 - r-gwas(): aggregate data in bed matrix ...
#> 2021-07-08 14:48:02 - r-gwas(): aggregate data in bed matrix DONE
#> 2021-07-08 14:48:02 - r-gwas(): remove samples with missing phenotypic values ...
#> 2021-07-08 14:48:02 - r-gwas(): remove samples with missing phenotypic values DONE
#> 2021-07-08 14:48:02 - r-gwas(): filter SNPs ...
#> 2021-07-08 14:48:02 - r-gwas(): filter SNPs DONE
#> 2021-07-08 14:48:02 - r-gwas(): fit model ...
#> [Iteration 1] theta = 78.0648 39.8993
#> [Iteration 1] log L = -1000.94
#> [Iteration 1] AI-REML update
#> [Iteration 1] ||gradient|| = 0.403712
#> [Iteration 2] theta = 19.3246 82.1856
#> [Iteration 2] log L = -991.683
#> [Iteration 2] AI-REML update
#> [Iteration 2] ||gradient|| = 1.19646
#> [Iteration 3] theta = 25.9636 94.4355
#> [Iteration 3] log L = -985.452
#> [Iteration 3] AI-REML update
#> [Iteration 3] ||gradient|| = 0.200588
#> [Iteration 4] theta = 29.0151 95.1047
#> [Iteration 4] log L = -985.125
#> [Iteration 4] AI-REML update
#> [Iteration 4] ||gradient|| = 0.0176334
#> [Iteration 5] theta = 29.5482  94.409
#> [Iteration 5] log L = -985.121
#> [Iteration 5] AI-REML update
#> [Iteration 5] ||gradient|| = 0.00136693
#> [Iteration 6] theta = 29.6081 94.2813
#> [Iteration 6] log L = -985.121
#> [Iteration 6] AI-REML update
#> [Iteration 6] ||gradient|| = 0.000156205
#> [Iteration 7] theta = 29.6153 94.2651
#> [Iteration 7] log L = -985.121
#> [Iteration 7] AI-REML update
#> [Iteration 7] ||gradient|| = 1.8962e-05
#> [Iteration 8] theta = 29.6162 94.2631
#> [Iteration 8] log L = -985.121
#> [Iteration 8] AI-REML update
#> [Iteration 8] ||gradient|| = 2.31736e-06
#> 2021-07-08 14:48:03 - r-gwas(): fit model DONE
#> 2021-07-08 14:48:03 - r-gwas(): DONE, return output.
#> 2021-07-08 14:48:03 - r-run_gwas(): GWAS analysis DONE
#> 2021-07-08 14:48:03 - r-run_gwas(): Save metadata ...
#> 2021-07-08 14:48:03 - r-run_gwas(): Save metadata DONE
#> 2021-07-08 14:48:03 - r-run_gwas(): Save results ...
#> 2021-07-08 14:48:03 - r-saveGWAS(): Check dir ...
#> 2021-07-08 14:48:03 - r-saveGWAS(): Check dir DONE
#> 2021-07-08 14:48:04 - r-run_gwas(): Save results DONE
gwas_results$file
#> [1] "/tmp/RtmpwjRlyt/file3e10b3698ad19.json"
substr(gwas_results$gwasRes, start=1, stop=500)
#> [
#>   {
#>     "chr": "1",
#>     "pos": 9563,
#>     "id": "SNP-1.8562.",
#>     "A1": "A",
#>     "A2": "T",
#>     "freqA2": 0.1359,
#>     "score": 3.6664,
#>     "p": 0.0555
#>   },
#>   {
#>     "chr": "1",
#>     "pos": 25922,
#>     "id": "SNP-1.24921.",
#>     "A1": "C",
#>     "A2": "T",
#>     "freqA2": 0.1254,
#>     "score": 1.3876,
#>     "p": 0.2388
#>   },
#>   {
#>     "chr": "1",
#>     "pos": 26254,
#>     "id": "SNP-1.25253.",
#>     "A1": "A",
#>     "A2": "T",
#>     "freqA2": 0.2935,
#>     "score": 2.9466,
#>     "p": 0.0861
#>   },
#>   {
#>     "chr": "1",
#>     "p
```

### Draw Manhattan Plot

``` r
p <- draw_manhattanPlot(gwasFile = gwas_results$file,
                        gwasUrl = NULL,
                        adj_method = "bonferroni",
                        thresh_p = 0.05,
                        chr = NA)
#> 2021-07-08 14:48:04 - r-draw_manhattanPlot(): Get data ...
#> 2021-07-08 14:48:04 - r-readGWAS(): Read result file ... 
#> 2021-07-08 14:48:04 - r-readGWAS(): Read result file DONE 
#> 2021-07-08 14:48:04 - r-readGWAS(): Convert Json to data.frame ... 
#> 2021-07-08 14:48:05 - r-readGWAS(): Convert Json to data.frame DONE 
#> 2021-07-08 14:48:05 - r-readGWAS(): DONE, return output.
#> 2021-07-08 14:48:05 - r-draw_manhattanPlot(): Get data DONE
#> 2021-07-08 14:48:05 - r-draw_manhattanPlot(): Draw Manhattan Plot ...
#> 2021-07-08 14:48:05 - r-manPlot(): Adjust p-values ...
#> 2021-07-08 14:48:05 - r-adjustPval(): Check adj_method ...
#> 2021-07-08 14:48:05 - r-adjustPval(): Check adj_method DONE
#> 2021-07-08 14:48:05 - r-adjustPval(): Adjust p-values ...
#> 2021-07-08 14:48:05 - r-adjustPval(): Adjust p-values DONE
#> 2021-07-08 14:48:05 - r-adjustPval(): Adjust threshold ...
#> 2021-07-08 14:48:05 - r-adjustPval(): Adjust threshold DONE
#> 2021-07-08 14:48:05 - r-adjustPval(): DONE, return output
#> 2021-07-08 14:48:05 - r-manPlot(): Adjust p-values DONE
#> 2021-07-08 14:48:06 - r-manPlot(): DONE, return output
#> 2021-07-08 14:48:06 - r-draw_manhattanPlot(): Draw Manhattan Plot DONE
```

![screenshot of the plotly graph](README_files/manPlot.png)

### Adjust p-values

``` r
gwas_adj <- run_resAdjustment(gwasFile = gwas_results$file,
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              dir = tempdir())
#> 2021-07-08 14:48:06 - r-run_resAdjustment(): Get data ...
#> 2021-07-08 14:48:06 - r-readGWAS(): Read result file ... 
#> 2021-07-08 14:48:06 - r-readGWAS(): Read result file DONE 
#> 2021-07-08 14:48:06 - r-readGWAS(): Convert Json to data.frame ... 
#> 2021-07-08 14:48:06 - r-readGWAS(): Convert Json to data.frame DONE 
#> 2021-07-08 14:48:06 - r-readGWAS(): DONE, return output.
#> 2021-07-08 14:48:06 - r-run_resAdjustment(): Get data DONE
#> 2021-07-08 14:48:06 - r-run_resAdjustment(): Adjust p-values ...
#> 2021-07-08 14:48:06 - r-adjustPval(): Check adj_method ...
#> 2021-07-08 14:48:06 - r-adjustPval(): Check adj_method DONE
#> 2021-07-08 14:48:06 - r-adjustPval(): Adjust p-values ...
#> 2021-07-08 14:48:06 - r-adjustPval(): Adjust p-values DONE
#> 2021-07-08 14:48:06 - r-adjustPval(): DONE, return output
#> 2021-07-08 14:48:06 - r-run_resAdjustment(): Adjust p-values DONE
#> 2021-07-08 14:48:06 - r-run_resAdjustment(): Save results ...
#> 2021-07-08 14:48:06 - r-saveGWAS(): Check dir ...
#> 2021-07-08 14:48:06 - r-saveGWAS(): Check dir DONE
#> 2021-07-08 14:48:06 - r-run_resAdjustment(): Save results DONE
substr(gwas_adj$gwasAdjusted, start=1, stop=500)
#> {
#>   "p_adj": []
#> }
```

### Draw LD plot

``` r
imgFile <- draw_ldPlot(genoFile = "data/geno/testMarkerData01.vcf.gz",
                       genoUrl = NULL,
                       from = 42,
                       to = 62,
                       dir = tempdir()) 
#> 2021-07-08 14:48:06 - r-draw_ldPlot(): Get data ...
#> 2021-07-08 14:48:06 - r-readGenoData(): Check file extention ... 
#> 2021-07-08 14:48:06 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-07-08 14:48:08 - r-readGenoData(): Read geno file DONE 
#> 2021-07-08 14:48:08 - r-readGenoData(): DONE, return output.
#> 2021-07-08 14:48:08 - r-draw_ldPlot(): Get data DONE
#> 2021-07-08 14:48:08 - r-draw_ldPlot(): Draw LD Plot ...
#> 2021-07-08 14:48:08 - r-LDplot(): Check "from" < "to"...
#> 2021-07-08 14:48:08 - r-LDplot(): Check "from" < "to" DONE
#> 2021-07-08 14:48:08 - r-LDplot(): Check number of SNP < 50...
#> 2021-07-08 14:48:08 - r-LDplot(): Check number of SNP < 50 DONE
#> 2021-07-08 14:48:08 - r-LDplot(): Check dir ...
#> 2021-07-08 14:48:08 - r-LDplot(): Check dir DONE
#> 2021-07-08 14:48:08 - r-LDplot(): Compute LD ...
#> 2021-07-08 14:48:08 - r-LDplot(): Compute LD DONE
#> 2021-07-08 14:48:08 - r-LDplot(): Create LD plot ...
#> 2021-07-08 14:48:08 - r-LDplot(): Create create file: /tmp/RtmpwjRlyt/file3e10b6a7b3839.png
#> 2021-07-08 14:48:08 - r-LDplot(): Create LD plot DONE
#> 2021-07-08 14:48:08 - r-LDplot(): DONE, return output
#> 2021-07-08 14:48:08 - r-draw_ldPlot(): Draw LD Plot DONE
```

``` r
library(png)
img <- readPNG(imgFile)
grid::grid.raster(img)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Load data

Example data are stored in the `data` folder.

`readData` read and prepare the data for GWAS analysis

``` r
data <- readData(genoFile = "data/geno/testMarkerData01.vcf",
                 phenoFile = "data/pheno/testPhenoData01.csv")
#> 2021-07-08 14:48:10 - r-readData(): get geno data ...
#> 2021-07-08 14:48:10 - r-readGenoData(): Check file extention ... 
#> 2021-07-08 14:48:10 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-07-08 14:48:11 - r-readGenoData(): Read geno file DONE 
#> 2021-07-08 14:48:11 - r-readGenoData(): DONE, return output.
#> 2021-07-08 14:48:11 - r-readData(): get geno data DONE
#> 2021-07-08 14:48:11 - r-readData(): get pheno data ...
#> 2021-07-08 14:48:11 - r-readPhenoData(): Read phenotypic file ... 
#> 2021-07-08 14:48:11 - r-readPhenoData(): Read phenotypic file DONE 
#> 2021-07-08 14:48:11 - r-readPhenoData(): DONE, return output.
#> 2021-07-08 14:48:11 - r-readData(): get pheno data DONE
#> 2021-07-08 14:48:11 - r-readData(): prepare data ...
#> 2021-07-08 14:48:11 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set ...
#> 2021-07-08 14:48:12 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set DONE
#> 2021-07-08 14:48:12 - r-prepareData(): reorder matrix ...
#> 2021-07-08 14:48:12 - r-prepareData(): reorder matrix DONE
#> 2021-07-08 14:48:12 - r-prepareData(): remove monomorphic geno ...
#> 2021-07-08 14:48:12 - r-prepareData(): remove monomorphic geno DONE
#> 2021-07-08 14:48:12 - r-prepareData(): calculate genetic relatinoal matrix ...
#> 2021-07-08 14:48:12 - r-prepareData(): calculate genetic relatinoal matrix DONE
#> 2021-07-08 14:48:12 - r-prepareData(): DONE, return output.
#> 2021-07-08 14:48:12 - r-readData(): prepare data DONE
#> 2021-07-08 14:48:12 - r-readData(): DONE, return output.
```

For online data, they can be loaded using:

``` r
data <- downloadData(genoUrl = "url/to/geno_data.vcf",
                     phenoUrl = "url/to/pheno_data.csv")
```

## GWAS

The `gwas` function return a `data.frame`.

``` r
gwasRes <- gwas(data = data,
                trait = "Flowering.time.at.Arkansas",
                test = "score",
                fixed = 0,
                response = "quantitative",
                thresh_maf = 0.05,
                thresh_callrate = 0.95)
#> 2021-07-08 14:48:12 - r-gwas(): Check inputs ...
#> 2021-07-08 14:48:12 - r-gwas(): Check inputs DONE
#> 2021-07-08 14:48:12 - r-gwas(): aggregate data in bed matrix ...
#> 2021-07-08 14:48:12 - r-gwas(): aggregate data in bed matrix DONE
#> 2021-07-08 14:48:12 - r-gwas(): remove samples with missing phenotypic values ...
#> 2021-07-08 14:48:13 - r-gwas(): remove samples with missing phenotypic values DONE
#> 2021-07-08 14:48:13 - r-gwas(): filter SNPs ...
#> 2021-07-08 14:48:13 - r-gwas(): filter SNPs DONE
#> 2021-07-08 14:48:13 - r-gwas(): fit model ...
#> [Iteration 1] theta = 78.0648 39.8994
#> [Iteration 1] log L = -1000.94
#> [Iteration 1] AI-REML update
#> [Iteration 1] ||gradient|| = 0.403711
#> [Iteration 2] theta = 19.3248 82.1855
#> [Iteration 2] log L = -991.683
#> [Iteration 2] AI-REML update
#> [Iteration 2] ||gradient|| = 1.19644
#> [Iteration 3] theta = 25.9638 94.4354
#> [Iteration 3] log L = -985.452
#> [Iteration 3] AI-REML update
#> [Iteration 3] ||gradient|| = 0.200582
#> [Iteration 4] theta = 29.0152 95.1046
#> [Iteration 4] log L = -985.126
#> [Iteration 4] AI-REML update
#> [Iteration 4] ||gradient|| = 0.0176325
#> [Iteration 5] theta = 29.5484 94.4089
#> [Iteration 5] log L = -985.121
#> [Iteration 5] AI-REML update
#> [Iteration 5] ||gradient|| = 0.00136686
#> [Iteration 6] theta = 29.6083 94.2812
#> [Iteration 6] log L = -985.121
#> [Iteration 6] AI-REML update
#> [Iteration 6] ||gradient|| = 0.000156197
#> [Iteration 7] theta = 29.6155  94.265
#> [Iteration 7] log L = -985.121
#> [Iteration 7] AI-REML update
#> [Iteration 7] ||gradient|| = 1.8961e-05
#> [Iteration 8] theta = 29.6164  94.263
#> [Iteration 8] log L = -985.121
#> [Iteration 8] AI-REML update
#> [Iteration 8] ||gradient|| = 2.31723e-06
#> 2021-07-08 14:48:13 - r-gwas(): fit model DONE
#> 2021-07-08 14:48:13 - r-gwas(): DONE, return output.
head(gwasRes)
#>   chr   pos           id A1 A2    freqA2    score          p
#> 1   1  9563  SNP-1.8562.  A  T 0.1358543 3.666395 0.05552015
#> 2   1 25922 SNP-1.24921.  C  T 0.1253521 1.387576 0.23881494
#> 3   1 26254 SNP-1.25253.  A  T 0.2935393 2.946588 0.08605908
#> 4   1 30214 SNP-1.29213.  T  A 0.1320225 2.189340 0.13896890
#> 5   1 31478 SNP-1.30477.  C  T 0.2268908 1.373070 0.24128505
#> 6   1 32733 SNP-1.31732.  T  G 0.2991573 3.335746 0.06778966
```

Save GWAS result in a `.json` file.

``` r
(file <- saveGWAS(gwasRes, metadata = "README"))
#> 2021-07-08 14:48:13 - r-saveGWAS(): Check dir ...
#> 2021-07-08 14:48:13 - r-saveGWAS(): Check dir DONE
#> [1] "/tmp/RtmpwjRlyt/file3e10b31abd8e3.json"
cat(paste(readLines(file, n=20), collapse = "\n"))
#> {
#>   "gwas": [
#>     {
#>       "chr": "1",
#>       "pos": 9563,
#>       "id": "SNP-1.8562.",
#>       "A1": "A",
#>       "A2": "T",
#>       "freqA2": 0.135854341736695,
#>       "score": 3.66639500551692,
#>       "p": 0.0555201477706262
#>     },
#>     {
#>       "chr": "1",
#>       "pos": 25922,
#>       "id": "SNP-1.24921.",
#>       "A1": "C",
#>       "A2": "T",
#>       "freqA2": 0.125352112676056,
#>       "score": 1.38757564950292,
```

## Manhattan plot

This function generates a “plotly” graph.

``` r
p <- manPlot(gwas = gwasRes,
             adj_method = "bonferroni",
             thresh_p = 0.05,
             chr = NA,
             title = "Readme Example")
#> 2021-07-08 14:48:14 - r-manPlot(): Adjust p-values ...
#> 2021-07-08 14:48:14 - r-adjustPval(): Check adj_method ...
#> 2021-07-08 14:48:14 - r-adjustPval(): Check adj_method DONE
#> 2021-07-08 14:48:14 - r-adjustPval(): Adjust p-values ...
#> 2021-07-08 14:48:14 - r-adjustPval(): Adjust p-values DONE
#> 2021-07-08 14:48:14 - r-adjustPval(): Adjust threshold ...
#> 2021-07-08 14:48:14 - r-adjustPval(): Adjust threshold DONE
#> 2021-07-08 14:48:14 - r-adjustPval(): DONE, return output
#> 2021-07-08 14:48:14 - r-manPlot(): Adjust p-values DONE
#> 2021-07-08 14:48:14 - r-manPlot(): DONE, return output
```

![screenshot of the plotly graph](README_files/manPlot2.png)

## LD PLot

Compute *r*<sup>2</sup> Linkage Disequilibrium (LD) between given SNPs
and save a plot in a temporary PNG file.

``` r
gDta <- readGenoData("data/geno/testMarkerData01.vcf.gz")
#> 2021-07-08 14:48:14 - r-readGenoData(): Check file extention ... 
#> 2021-07-08 14:48:14 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-07-08 14:48:15 - r-readGenoData(): Read geno file DONE 
#> 2021-07-08 14:48:15 - r-readGenoData(): DONE, return output.
(imgFile <- LDplot(geno = gDta,
                   from = 1,
                   to = 20,
                   dir = tempdir()))
#> 2021-07-08 14:48:15 - r-LDplot(): Check "from" < "to"...
#> 2021-07-08 14:48:15 - r-LDplot(): Check "from" < "to" DONE
#> 2021-07-08 14:48:15 - r-LDplot(): Check number of SNP < 50...
#> 2021-07-08 14:48:15 - r-LDplot(): Check number of SNP < 50 DONE
#> 2021-07-08 14:48:15 - r-LDplot(): Check dir ...
#> 2021-07-08 14:48:15 - r-LDplot(): Check dir DONE
#> 2021-07-08 14:48:15 - r-LDplot(): Compute LD ...
#> 2021-07-08 14:48:15 - r-LDplot(): Compute LD DONE
#> 2021-07-08 14:48:15 - r-LDplot(): Create LD plot ...
#> 2021-07-08 14:48:15 - r-LDplot(): Create create file: /tmp/RtmpwjRlyt/file3e10b7f14c9fc.png
#> 2021-07-08 14:48:16 - r-LDplot(): Create LD plot DONE
#> 2021-07-08 14:48:16 - r-LDplot(): DONE, return output
#> [1] "/tmp/RtmpwjRlyt/file3e10b7f14c9fc.png"
```

``` r
img <- readPNG(imgFile)
grid::grid.raster(img)
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Functions documentation

All the functions of this engine use `roxygen2` for documentation. You
can generate a markdown documentation located in `./doc/README.md` using
the function `writeDoc` (defined in `./src/utils.R`):

``` r
writeDoc(srcDir = "src",
         docDir = "doc")
```

# Tests

The R package `testthat` is needed to run the unit tests. To run the
unit tests of this engine you can use the command:

``` sh
Rscript ./tests/testthat.R
```
