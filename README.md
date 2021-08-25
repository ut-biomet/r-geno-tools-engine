
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GWAS-Engine

<!-- badges: start -->

<!-- badges: end -->

The goal of GWAS-Engine is to provide simple “input/output” style code
for GWAS analysis in order to be used inside an API.

# How to use GWAS-Engine

## Dependencies

Package dependencies are managed using
[Renv](https://rstudio.github.io/renv/articles/renv.html). To install
dependencies do in R console: `renv::restore()`.

``` r
# required packages for using GWAS-Engine:
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

``` sh
Rscript ./gwas-engine.R gwas --genoFile "data/geno/testMarkerData01.vcf.gz" --phenoFile "data/pheno/testPhenoData01.csv" --trait "Flowering.time.at.Arkansas" --test "score" --fixed 0 --response "quantitative" --thresh_maf 0.05 --thresh_callrate 0.95 --outFile readmeTemp/gwasRes.json
#> Loading required package: Rcpp
#> Loading required package: RcppParallel
#> 
#> Attaching package: 'RcppParallel'
#> 
#> The following object is masked from 'package:Rcpp':
#> 
#>     LdFlags
#> 
#> Gaston set number of threads to 4. Use setThreadOptions() to modify this.
#> 
#> Attaching package: 'gaston'
#> 
#> The following object is masked from 'package:stats':
#> 
#>     sigma
#> 
#> The following objects are masked from 'package:base':
#> 
#>     cbind, rbind
#> 
#> See example usage at http://sahirbhatnagar.com/manhattanly/
#> 2021-08-25 15:06:39 - r-run_gwas(): Get data ...
#> 2021-08-25 15:06:39 - r-readData(): get geno data ...
#> 2021-08-25 15:06:39 - r-readGenoData(): Check file extention ... 
#> 2021-08-25 15:06:39 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-08-25 15:06:40 - r-readGenoData(): Read geno file DONE 
#> 2021-08-25 15:06:40 - r-readGenoData(): DONE, return output.
#> 2021-08-25 15:06:40 - r-readData(): get geno data DONE
#> 2021-08-25 15:06:40 - r-readData(): get pheno data ...
#> 2021-08-25 15:06:40 - r-readPhenoData(): Read phenotypic file ... 
#> 2021-08-25 15:06:40 - r-readPhenoData(): Read phenotypic file DONE 
#> 2021-08-25 15:06:40 - r-readPhenoData(): DONE, return output.
#> 2021-08-25 15:06:40 - r-readData(): get pheno data DONE
#> 2021-08-25 15:06:40 - r-readData(): prepare data ...
#> 2021-08-25 15:06:40 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set ...
#> 2021-08-25 15:06:40 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set DONE
#> 2021-08-25 15:06:40 - r-prepareData(): reorder matrix ...
#> 2021-08-25 15:06:40 - r-prepareData(): reorder matrix DONE
#> 2021-08-25 15:06:40 - r-prepareData(): remove monomorphic geno ...
#> 2021-08-25 15:06:40 - r-prepareData(): remove monomorphic geno DONE
#> 2021-08-25 15:06:40 - r-prepareData(): calculate genetic relatinoal matrix ...
#> 2021-08-25 15:06:40 - r-prepareData(): calculate genetic relatinoal matrix DONE
#> 2021-08-25 15:06:40 - r-prepareData(): DONE, return output.
#> 2021-08-25 15:06:40 - r-readData(): prepare data DONE
#> 2021-08-25 15:06:40 - r-readData(): DONE, return output.
#> 2021-08-25 15:06:40 - r-run_gwas(): Get data DONE
#> 2021-08-25 15:06:40 - r-run_gwas(): GWAS analysis ...
#> 2021-08-25 15:06:40 - r-gwas(): Check inputs ...
#> 2021-08-25 15:06:40 - r-gwas(): Check inputs DONE
#> 2021-08-25 15:06:40 - r-gwas(): aggregate data in bed matrix ...
#> 2021-08-25 15:06:40 - r-gwas(): aggregate data in bed matrix DONE
#> 2021-08-25 15:06:40 - r-gwas(): remove individuals with missing phenotypic values ...
#> 2021-08-25 15:06:41 - r-gwas(): remove samples with missing phenotypic values DONE
#> 2021-08-25 15:06:41 - r-gwas(): filter SNPs ...
#> 2021-08-25 15:06:41 - r-gwas(): filter SNPs DONE
#> 2021-08-25 15:06:41 - r-gwas(): fit model ...
#> [Iteration 1] theta = 78.0648 39.8994
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
#> [Iteration 3] ||gradient|| = 0.200589
#> [Iteration 4] theta = 29.0151 95.1047
#> [Iteration 4] log L = -985.125
#> [Iteration 4] AI-REML update
#> [Iteration 4] ||gradient|| = 0.0176334
#> [Iteration 5] theta = 29.5483  94.409
#> [Iteration 5] log L = -985.121
#> [Iteration 5] AI-REML update
#> [Iteration 5] ||gradient|| = 0.00136693
#> [Iteration 6] theta = 29.6082 94.2812
#> [Iteration 6] log L = -985.12
#> [Iteration 6] AI-REML update
#> [Iteration 6] ||gradient|| = 0.000156204
#> [Iteration 7] theta = 29.6154  94.265
#> [Iteration 7] log L = -985.12
#> [Iteration 7] AI-REML update
#> [Iteration 7] ||gradient|| = 1.89618e-05
#> [Iteration 8] theta = 29.6163 94.2631
#> [Iteration 8] log L = -985.12
#> [Iteration 8] AI-REML update
#> [Iteration 8] ||gradient|| = 2.31732e-06
#> 2021-08-25 15:06:41 - r-gwas(): fit model DONE
#> 2021-08-25 15:06:41 - r-gwas(): DONE, return output.
#> 2021-08-25 15:06:41 - r-run_gwas(): GWAS analysis DONE
#> 2021-08-25 15:06:41 - r-run_gwas(): Save metadata ...
#> 2021-08-25 15:06:41 - r-run_gwas(): Save metadata DONE
#> 2021-08-25 15:06:41 - r-run_gwas(): Save results ...
#> 2021-08-25 15:06:41 - r-saveGWAS(): Check file ...
#> 2021-08-25 15:06:41 - r-saveGWAS(): Check file DONE
#> 2021-08-25 15:06:41 - r-run_gwas(): Save results DONE
```

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
                         outFile = tempfile(fileext = ".json"))
#> 2021-08-25 15:06:41 - r-run_gwas(): Get data ...
#> 2021-08-25 15:06:41 - r-readData(): get geno data ...
#> 2021-08-25 15:06:41 - r-readGenoData(): Check file extention ... 
#> 2021-08-25 15:06:41 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-08-25 15:06:42 - r-readGenoData(): Read geno file DONE 
#> 2021-08-25 15:06:42 - r-readGenoData(): DONE, return output.
#> 2021-08-25 15:06:42 - r-readData(): get geno data DONE
#> 2021-08-25 15:06:42 - r-readData(): get pheno data ...
#> 2021-08-25 15:06:42 - r-readPhenoData(): Read phenotypic file ... 
#> 2021-08-25 15:06:42 - r-readPhenoData(): Read phenotypic file DONE 
#> 2021-08-25 15:06:42 - r-readPhenoData(): DONE, return output.
#> 2021-08-25 15:06:42 - r-readData(): get pheno data DONE
#> 2021-08-25 15:06:42 - r-readData(): prepare data ...
#> 2021-08-25 15:06:42 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set ...
#> 2021-08-25 15:06:42 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set DONE
#> 2021-08-25 15:06:42 - r-prepareData(): reorder matrix ...
#> 2021-08-25 15:06:42 - r-prepareData(): reorder matrix DONE
#> 2021-08-25 15:06:42 - r-prepareData(): remove monomorphic geno ...
#> 2021-08-25 15:06:42 - r-prepareData(): remove monomorphic geno DONE
#> 2021-08-25 15:06:42 - r-prepareData(): calculate genetic relatinoal matrix ...
#> 2021-08-25 15:06:43 - r-prepareData(): calculate genetic relatinoal matrix DONE
#> 2021-08-25 15:06:43 - r-prepareData(): DONE, return output.
#> 2021-08-25 15:06:43 - r-readData(): prepare data DONE
#> 2021-08-25 15:06:43 - r-readData(): DONE, return output.
#> 2021-08-25 15:06:43 - r-run_gwas(): Get data DONE
#> 2021-08-25 15:06:43 - r-run_gwas(): GWAS analysis ...
#> 2021-08-25 15:06:43 - r-gwas(): Check inputs ...
#> 2021-08-25 15:06:43 - r-gwas(): Check inputs DONE
#> 2021-08-25 15:06:43 - r-gwas(): aggregate data in bed matrix ...
#> 2021-08-25 15:06:43 - r-gwas(): aggregate data in bed matrix DONE
#> 2021-08-25 15:06:43 - r-gwas(): remove individuals with missing phenotypic values ...
#> 2021-08-25 15:06:43 - r-gwas(): remove samples with missing phenotypic values DONE
#> 2021-08-25 15:06:43 - r-gwas(): filter SNPs ...
#> 2021-08-25 15:06:43 - r-gwas(): filter SNPs DONE
#> 2021-08-25 15:06:43 - r-gwas(): fit model ...
#> [Iteration 1] theta = 78.0648 39.8993
#> [Iteration 1] log L = -1000.94
#> [Iteration 1] AI-REML update
#> [Iteration 1] ||gradient|| = 0.403711
#> [Iteration 2] theta = 19.3246 82.1855
#> [Iteration 2] log L = -991.683
#> [Iteration 2] AI-REML update
#> [Iteration 2] ||gradient|| = 1.19646
#> [Iteration 3] theta = 25.9636 94.4353
#> [Iteration 3] log L = -985.452
#> [Iteration 3] AI-REML update
#> [Iteration 3] ||gradient|| = 0.200589
#> [Iteration 4] theta = 29.0151 95.1044
#> [Iteration 4] log L = -985.126
#> [Iteration 4] AI-REML update
#> [Iteration 4] ||gradient|| = 0.0176335
#> [Iteration 5] theta = 29.5483 94.4087
#> [Iteration 5] log L = -985.121
#> [Iteration 5] AI-REML update
#> [Iteration 5] ||gradient|| = 0.00136694
#> [Iteration 6] theta = 29.6082  94.281
#> [Iteration 6] log L = -985.121
#> [Iteration 6] AI-REML update
#> [Iteration 6] ||gradient|| = 0.000156205
#> [Iteration 7] theta = 29.6154 94.2648
#> [Iteration 7] log L = -985.121
#> [Iteration 7] AI-REML update
#> [Iteration 7] ||gradient|| = 1.89618e-05
#> [Iteration 8] theta = 29.6163 94.2628
#> [Iteration 8] log L = -985.121
#> [Iteration 8] AI-REML update
#> [Iteration 8] ||gradient|| = 2.31731e-06
#> 2021-08-25 15:06:43 - r-gwas(): fit model DONE
#> 2021-08-25 15:06:43 - r-gwas(): DONE, return output.
#> 2021-08-25 15:06:43 - r-run_gwas(): GWAS analysis DONE
#> 2021-08-25 15:06:43 - r-run_gwas(): Save metadata ...
#> 2021-08-25 15:06:43 - r-run_gwas(): Save metadata DONE
#> 2021-08-25 15:06:43 - r-run_gwas(): Save results ...
#> 2021-08-25 15:06:43 - r-saveGWAS(): Check file ...
#> 2021-08-25 15:06:43 - r-saveGWAS(): Check file DONE
#> 2021-08-25 15:06:44 - r-run_gwas(): Save results DONE
gwas_results$file
#> [1] "/tmp/RtmpAxU90K/file3b7122da07697.json"
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

``` sh
Rscript ./gwas-engine.R manplot --gwasFile readmeTemp/gwasRes.json --adj_method "bonferroni" --thresh_p 0.05 --outFile readmeTemp/manPlot.html
#> Loading required package: Rcpp
#> Loading required package: RcppParallel
#> 
#> Attaching package: 'RcppParallel'
#> 
#> The following object is masked from 'package:Rcpp':
#> 
#>     LdFlags
#> 
#> Gaston set number of threads to 4. Use setThreadOptions() to modify this.
#> 
#> Attaching package: 'gaston'
#> 
#> The following object is masked from 'package:stats':
#> 
#>     sigma
#> 
#> The following objects are masked from 'package:base':
#> 
#>     cbind, rbind
#> 
#> See example usage at http://sahirbhatnagar.com/manhattanly/
#> 2021-08-25 15:06:44 - r-draw_manhattanPlot(): Get data ...
#> 2021-08-25 15:06:44 - r-readGWAS(): Read result file ... 
#> 2021-08-25 15:06:44 - r-readGWAS(): Read result file DONE 
#> 2021-08-25 15:06:44 - r-readGWAS(): Convert Json to data.frame ... 
#> 2021-08-25 15:06:45 - r-readGWAS(): Convert Json to data.frame DONE 
#> 2021-08-25 15:06:45 - r-readGWAS(): DONE, return output.
#> 2021-08-25 15:06:45 - r-draw_manhattanPlot(): Get data DONE
#> 2021-08-25 15:06:45 - r-draw_manhattanPlot(): Draw Manhattan Plot ...
#> 2021-08-25 15:06:45 - r-manPlot(): Check chromosome name ...
#> 2021-08-25 15:06:45 - r-manPlot(): Check chromosome name DONE
#> 2021-08-25 15:06:45 - r-manPlot(): Adjust p-values ...
#> 2021-08-25 15:06:45 - r-adjustPval(): Check adj_method ...
#> 2021-08-25 15:06:45 - r-adjustPval(): Check adj_method DONE
#> 2021-08-25 15:06:45 - r-adjustPval(): Check p values ...
#> 2021-08-25 15:06:45 - r-adjustPval(): Check p values DONE
#> 2021-08-25 15:06:45 - r-adjustPval(): Adjust p-values ...
#> 2021-08-25 15:06:45 - r-adjustPval(): Adjust p-values DONE
#> 2021-08-25 15:06:45 - r-adjustPval(): Adjust threshold ...
#> 2021-08-25 15:06:45 - r-adjustPval(): Adjust threshold DONE
#> 2021-08-25 15:06:45 - r-adjustPval(): DONE, return output
#> 2021-08-25 15:06:45 - r-manPlot(): Adjust p-values DONE
#> 2021-08-25 15:06:45 - r-manPlot(): DONE, return output
#> 2021-08-25 15:06:45 - r-draw_manhattanPlot(): Draw Manhattan Plot DONE
#> 2021-08-25 15:06:45 - r-draw_manhattanPlot(): Save results ...
#> 2021-08-25 15:06:50 - r-draw_manhattanPlot(): Save results DONE
```

``` r
p <- draw_manhattanPlot(gwasFile = gwas_results$file,
                        gwasUrl = NULL,
                        adj_method = "bonferroni",
                        thresh_p = 0.05,
                        chr = NA,
                        outFile = tempfile(fileext = ".html"))
#> 2021-08-25 15:06:50 - r-draw_manhattanPlot(): Get data ...
#> 2021-08-25 15:06:50 - r-readGWAS(): Read result file ... 
#> 2021-08-25 15:06:50 - r-readGWAS(): Read result file DONE 
#> 2021-08-25 15:06:50 - r-readGWAS(): Convert Json to data.frame ... 
#> 2021-08-25 15:06:50 - r-readGWAS(): Convert Json to data.frame DONE 
#> 2021-08-25 15:06:50 - r-readGWAS(): DONE, return output.
#> 2021-08-25 15:06:50 - r-draw_manhattanPlot(): Get data DONE
#> 2021-08-25 15:06:50 - r-draw_manhattanPlot(): Draw Manhattan Plot ...
#> 2021-08-25 15:06:50 - r-manPlot(): Check chromosome name ...
#> 2021-08-25 15:06:50 - r-manPlot(): Check chromosome name DONE
#> 2021-08-25 15:06:50 - r-manPlot(): Adjust p-values ...
#> 2021-08-25 15:06:50 - r-adjustPval(): Check adj_method ...
#> 2021-08-25 15:06:50 - r-adjustPval(): Check adj_method DONE
#> 2021-08-25 15:06:50 - r-adjustPval(): Check p values ...
#> 2021-08-25 15:06:50 - r-adjustPval(): Check p values DONE
#> 2021-08-25 15:06:50 - r-adjustPval(): Adjust p-values ...
#> 2021-08-25 15:06:50 - r-adjustPval(): Adjust p-values DONE
#> 2021-08-25 15:06:50 - r-adjustPval(): Adjust threshold ...
#> 2021-08-25 15:06:50 - r-adjustPval(): Adjust threshold DONE
#> 2021-08-25 15:06:50 - r-adjustPval(): DONE, return output
#> 2021-08-25 15:06:50 - r-manPlot(): Adjust p-values DONE
#> 2021-08-25 15:06:50 - r-manPlot(): DONE, return output
#> 2021-08-25 15:06:50 - r-draw_manhattanPlot(): Draw Manhattan Plot DONE
#> 2021-08-25 15:06:50 - r-draw_manhattanPlot(): Save results ...
#> 2021-08-25 15:06:54 - r-draw_manhattanPlot(): Save results DONE
```

![screenshot of the plotly graph](README_files/manPlot.png)

### Adjust p-values

``` sh
Rscript ./gwas-engine.R adjresults --gwasFile readmeTemp/gwasRes.json --adj_method "bonferroni" --outFile readmeTemp/adjRes.json
#> Loading required package: Rcpp
#> Loading required package: RcppParallel
#> 
#> Attaching package: 'RcppParallel'
#> 
#> The following object is masked from 'package:Rcpp':
#> 
#>     LdFlags
#> 
#> Gaston set number of threads to 4. Use setThreadOptions() to modify this.
#> 
#> Attaching package: 'gaston'
#> 
#> The following object is masked from 'package:stats':
#> 
#>     sigma
#> 
#> The following objects are masked from 'package:base':
#> 
#>     cbind, rbind
#> 
#> See example usage at http://sahirbhatnagar.com/manhattanly/
#> 2021-08-25 15:06:55 - r-run_resAdjustment(): Get data ...
#> 2021-08-25 15:06:55 - r-readGWAS(): Read result file ... 
#> 2021-08-25 15:06:55 - r-readGWAS(): Read result file DONE 
#> 2021-08-25 15:06:55 - r-readGWAS(): Convert Json to data.frame ... 
#> 2021-08-25 15:06:55 - r-readGWAS(): Convert Json to data.frame DONE 
#> 2021-08-25 15:06:55 - r-readGWAS(): DONE, return output.
#> 2021-08-25 15:06:55 - r-run_resAdjustment(): Get data DONE
#> 2021-08-25 15:06:55 - r-run_resAdjustment(): Adjust p-values ...
#> 2021-08-25 15:06:55 - r-adjustPval(): Check adj_method ...
#> 2021-08-25 15:06:55 - r-adjustPval(): Check adj_method DONE
#> 2021-08-25 15:06:55 - r-adjustPval(): Check p values ...
#> 2021-08-25 15:06:55 - r-adjustPval(): Check p values DONE
#> 2021-08-25 15:06:55 - r-adjustPval(): Adjust p-values ...
#> 2021-08-25 15:06:55 - r-adjustPval(): Adjust p-values DONE
#> 2021-08-25 15:06:55 - r-adjustPval(): DONE, return output
#> 2021-08-25 15:06:55 - r-run_resAdjustment(): Adjust p-values DONE
#> 2021-08-25 15:06:55 - r-run_resAdjustment(): Save results ...
#> 2021-08-25 15:06:55 - r-saveGWAS(): Check file ...
#> 2021-08-25 15:06:55 - r-saveGWAS(): Check file DONE
#> 2021-08-25 15:06:55 - r-run_resAdjustment(): Save results DONE
```

``` r
gwas_adj <- run_resAdjustment(gwasFile = gwas_results$file,
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              outFile = tempfile(fileext = ".json"))
#> 2021-08-25 15:06:55 - r-run_resAdjustment(): Get data ...
#> 2021-08-25 15:06:55 - r-readGWAS(): Read result file ... 
#> 2021-08-25 15:06:56 - r-readGWAS(): Read result file DONE 
#> 2021-08-25 15:06:56 - r-readGWAS(): Convert Json to data.frame ... 
#> 2021-08-25 15:06:56 - r-readGWAS(): Convert Json to data.frame DONE 
#> 2021-08-25 15:06:56 - r-readGWAS(): DONE, return output.
#> 2021-08-25 15:06:56 - r-run_resAdjustment(): Get data DONE
#> 2021-08-25 15:06:56 - r-run_resAdjustment(): Adjust p-values ...
#> 2021-08-25 15:06:56 - r-adjustPval(): Check adj_method ...
#> 2021-08-25 15:06:56 - r-adjustPval(): Check adj_method DONE
#> 2021-08-25 15:06:56 - r-adjustPval(): Check p values ...
#> 2021-08-25 15:06:56 - r-adjustPval(): Check p values DONE
#> 2021-08-25 15:06:56 - r-adjustPval(): Adjust p-values ...
#> 2021-08-25 15:06:56 - r-adjustPval(): Adjust p-values DONE
#> 2021-08-25 15:06:56 - r-adjustPval(): DONE, return output
#> 2021-08-25 15:06:56 - r-run_resAdjustment(): Adjust p-values DONE
#> 2021-08-25 15:06:56 - r-run_resAdjustment(): Save results ...
#> 2021-08-25 15:06:56 - r-saveGWAS(): Check file ...
#> 2021-08-25 15:06:56 - r-saveGWAS(): Check file DONE
#> 2021-08-25 15:06:56 - r-run_resAdjustment(): Save results DONE
substr(gwas_adj$gwasAdjusted, start=1, stop=500)
#> [
#>   {
#>     "chr": "1",
#>     "pos": 9563,
#>     "id": "SNP-1.8562.",
#>     "A1": "A",
#>     "A2": "T",
#>     "freqA2": 0.1359,
#>     "score": 3.6664,
#>     "p": 0.0555,
#>     "p_adj": 1
#>   },
#>   {
#>     "chr": "1",
#>     "pos": 25922,
#>     "id": "SNP-1.24921.",
#>     "A1": "C",
#>     "A2": "T",
#>     "freqA2": 0.1254,
#>     "score": 1.3876,
#>     "p": 0.2388,
#>     "p_adj": 1
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
```

### Draw LD plot

``` sh
Rscript ./gwas-engine.R ldplot --genoFile "data/geno/testMarkerData01.vcf.gz" --from 42 --to 62 --outFile readmeTemp/ldplot.png
#> Loading required package: Rcpp
#> Loading required package: RcppParallel
#> 
#> Attaching package: 'RcppParallel'
#> 
#> The following object is masked from 'package:Rcpp':
#> 
#>     LdFlags
#> 
#> Gaston set number of threads to 4. Use setThreadOptions() to modify this.
#> 
#> Attaching package: 'gaston'
#> 
#> The following object is masked from 'package:stats':
#> 
#>     sigma
#> 
#> The following objects are masked from 'package:base':
#> 
#>     cbind, rbind
#> 
#> See example usage at http://sahirbhatnagar.com/manhattanly/
#> 2021-08-25 15:06:56 - r-draw_ldPlot(): Get data ...
#> 2021-08-25 15:06:57 - r-readGenoData(): Check file extention ... 
#> 2021-08-25 15:06:57 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-08-25 15:06:57 - r-readGenoData(): Read geno file DONE 
#> 2021-08-25 15:06:57 - r-readGenoData(): DONE, return output.
#> 2021-08-25 15:06:57 - r-draw_ldPlot(): Get data DONE
#> 2021-08-25 15:06:57 - r-draw_ldPlot(): Draw LD Plot ...
#> 2021-08-25 15:06:57 - r-LDplot(): Check "from" and "to" format ...
#> 2021-08-25 15:06:57 - r-LDplot(): Check "from" and "to" format DONE
#> 2021-08-25 15:06:57 - r-LDplot(): Check "from" < "to"...
#> 2021-08-25 15:06:57 - r-LDplot(): Check "from" < "to" DONE
#> 2021-08-25 15:06:57 - r-LDplot(): Check number of SNP < 50...
#> 2021-08-25 15:06:57 - r-LDplot(): Check number of SNP < 50 DONE
#> 2021-08-25 15:06:57 - r-LDplot(): Check file ...
#> 2021-08-25 15:06:57 - r-LDplot(): Check file DONE
#> 2021-08-25 15:06:57 - r-LDplot(): Compute LD ...
#> 2021-08-25 15:06:57 - r-LDplot(): Compute LD DONE
#> 2021-08-25 15:06:57 - r-LDplot(): Create LD plot ...
#> 2021-08-25 15:06:57 - r-LDplot(): Create create file: readmeTemp/ldplot.png
#> 2021-08-25 15:06:58 - r-LDplot(): Create LD plot DONE
#> 2021-08-25 15:06:58 - r-LDplot(): DONE, return output
#> 2021-08-25 15:06:58 - r-draw_ldPlot(): Draw LD Plot DONE
```

``` r
imgFile <- draw_ldPlot(genoFile = "data/geno/testMarkerData01.vcf.gz",
                       genoUrl = NULL,
                       from = 42,
                       to = 62,
                       outFile = tempfile(fileext = ".png")) 
#> 2021-08-25 15:06:58 - r-draw_ldPlot(): Get data ...
#> 2021-08-25 15:06:58 - r-readGenoData(): Check file extention ... 
#> 2021-08-25 15:06:58 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-08-25 15:06:58 - r-readGenoData(): Read geno file DONE 
#> 2021-08-25 15:06:58 - r-readGenoData(): DONE, return output.
#> 2021-08-25 15:06:58 - r-draw_ldPlot(): Get data DONE
#> 2021-08-25 15:06:58 - r-draw_ldPlot(): Draw LD Plot ...
#> 2021-08-25 15:06:58 - r-LDplot(): Check "from" and "to" format ...
#> 2021-08-25 15:06:58 - r-LDplot(): Check "from" and "to" format DONE
#> 2021-08-25 15:06:58 - r-LDplot(): Check "from" < "to"...
#> 2021-08-25 15:06:58 - r-LDplot(): Check "from" < "to" DONE
#> 2021-08-25 15:06:58 - r-LDplot(): Check number of SNP < 50...
#> 2021-08-25 15:06:58 - r-LDplot(): Check number of SNP < 50 DONE
#> 2021-08-25 15:06:58 - r-LDplot(): Check file ...
#> 2021-08-25 15:06:58 - r-LDplot(): Check file DONE
#> 2021-08-25 15:06:58 - r-LDplot(): Compute LD ...
#> 2021-08-25 15:06:58 - r-LDplot(): Compute LD DONE
#> 2021-08-25 15:06:58 - r-LDplot(): Create LD plot ...
#> 2021-08-25 15:06:58 - r-LDplot(): Create create file: /tmp/RtmpAxU90K/file3b7122741cba4.png
#> 2021-08-25 15:06:59 - r-LDplot(): Create LD plot DONE
#> 2021-08-25 15:06:59 - r-LDplot(): DONE, return output
#> 2021-08-25 15:06:59 - r-draw_ldPlot(): Draw LD Plot DONE
```

``` r
library(png)
img <- readPNG(imgFile)
grid::grid.raster(img)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Load data

Example data are stored in the `data` folder.

`readData` read and prepare the data for GWAS analysis

``` r
data <- readData(genoFile = "data/geno/testMarkerData01.vcf",
                 phenoFile = "data/pheno/testPhenoData01.csv")
#> 2021-08-25 15:07:00 - r-readData(): get geno data ...
#> 2021-08-25 15:07:00 - r-readGenoData(): Check file extention ... 
#> 2021-08-25 15:07:00 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-08-25 15:07:00 - r-readGenoData(): Read geno file DONE 
#> 2021-08-25 15:07:00 - r-readGenoData(): DONE, return output.
#> 2021-08-25 15:07:00 - r-readData(): get geno data DONE
#> 2021-08-25 15:07:00 - r-readData(): get pheno data ...
#> 2021-08-25 15:07:00 - r-readPhenoData(): Read phenotypic file ... 
#> 2021-08-25 15:07:00 - r-readPhenoData(): Read phenotypic file DONE 
#> 2021-08-25 15:07:00 - r-readPhenoData(): DONE, return output.
#> 2021-08-25 15:07:00 - r-readData(): get pheno data DONE
#> 2021-08-25 15:07:00 - r-readData(): prepare data ...
#> 2021-08-25 15:07:00 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set ...
#> 2021-08-25 15:07:00 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set DONE
#> 2021-08-25 15:07:00 - r-prepareData(): reorder matrix ...
#> 2021-08-25 15:07:00 - r-prepareData(): reorder matrix DONE
#> 2021-08-25 15:07:00 - r-prepareData(): remove monomorphic geno ...
#> 2021-08-25 15:07:01 - r-prepareData(): remove monomorphic geno DONE
#> 2021-08-25 15:07:01 - r-prepareData(): calculate genetic relatinoal matrix ...
#> 2021-08-25 15:07:01 - r-prepareData(): calculate genetic relatinoal matrix DONE
#> 2021-08-25 15:07:01 - r-prepareData(): DONE, return output.
#> 2021-08-25 15:07:01 - r-readData(): prepare data DONE
#> 2021-08-25 15:07:01 - r-readData(): DONE, return output.
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
#> 2021-08-25 15:07:01 - r-gwas(): Check inputs ...
#> 2021-08-25 15:07:01 - r-gwas(): Check inputs DONE
#> 2021-08-25 15:07:01 - r-gwas(): aggregate data in bed matrix ...
#> 2021-08-25 15:07:01 - r-gwas(): aggregate data in bed matrix DONE
#> 2021-08-25 15:07:01 - r-gwas(): remove individuals with missing phenotypic values ...
#> 2021-08-25 15:07:01 - r-gwas(): remove samples with missing phenotypic values DONE
#> 2021-08-25 15:07:01 - r-gwas(): filter SNPs ...
#> 2021-08-25 15:07:01 - r-gwas(): filter SNPs DONE
#> 2021-08-25 15:07:01 - r-gwas(): fit model ...
#> [Iteration 1] theta = 78.0648 39.8993
#> [Iteration 1] log L = -1000.94
#> [Iteration 1] AI-REML update
#> [Iteration 1] ||gradient|| = 0.403711
#> [Iteration 2] theta = 19.3245 82.1855
#> [Iteration 2] log L = -991.683
#> [Iteration 2] AI-REML update
#> [Iteration 2] ||gradient|| = 1.19647
#> [Iteration 3] theta = 25.9636 94.4354
#> [Iteration 3] log L = -985.452
#> [Iteration 3] AI-REML update
#> [Iteration 3] ||gradient|| = 0.20059
#> [Iteration 4] theta = 29.0151 95.1045
#> [Iteration 4] log L = -985.126
#> [Iteration 4] AI-REML update
#> [Iteration 4] ||gradient|| = 0.0176337
#> [Iteration 5] theta = 29.5482 94.4088
#> [Iteration 5] log L = -985.121
#> [Iteration 5] AI-REML update
#> [Iteration 5] ||gradient|| = 0.00136696
#> [Iteration 6] theta = 29.6081  94.281
#> [Iteration 6] log L = -985.121
#> [Iteration 6] AI-REML update
#> [Iteration 6] ||gradient|| = 0.000156208
#> [Iteration 7] theta = 29.6154 94.2648
#> [Iteration 7] log L = -985.121
#> [Iteration 7] AI-REML update
#> [Iteration 7] ||gradient|| = 1.89622e-05
#> [Iteration 8] theta = 29.6162 94.2629
#> [Iteration 8] log L = -985.121
#> [Iteration 8] AI-REML update
#> [Iteration 8] ||gradient|| = 2.31737e-06
#> 2021-08-25 15:07:01 - r-gwas(): fit model DONE
#> 2021-08-25 15:07:01 - r-gwas(): DONE, return output.
head(gwasRes)
#>   chr   pos           id A1 A2    freqA2
#> 1   1  9563  SNP-1.8562.  A  T 0.1358543
#> 2   1 25922 SNP-1.24921.  C  T 0.1253521
#> 3   1 26254 SNP-1.25253.  A  T 0.2935393
#> 4   1 30214 SNP-1.29213.  T  A 0.1320225
#> 5   1 31478 SNP-1.30477.  C  T 0.2268908
#> 6   1 32733 SNP-1.31732.  T  G 0.2991573
#>      score          p
#> 1 3.666425 0.05551914
#> 2 1.387572 0.23881556
#> 3 2.946578 0.08605963
#> 4 2.189352 0.13896780
#> 5 1.373079 0.24128359
#> 6 3.335736 0.06779009
```

Save GWAS result in a `.json` file.

``` r
(file <- saveGWAS(gwasRes, metadata = "README"))
#> 2021-08-25 15:07:01 - r-saveGWAS(): Check dir ...
#> 2021-08-25 15:07:01 - r-saveGWAS(): Check dir DONE
#> [1] "/tmp/RtmpAxU90K/file3b712565bff11.json"
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
#>       "score": 3.66642523198947,
#>       "p": 0.055519140778278
#>     },
#>     {
#>       "chr": "1",
#>       "pos": 25922,
#>       "id": "SNP-1.24921.",
#>       "A1": "C",
#>       "A2": "T",
#>       "freqA2": 0.125352112676056,
#>       "score": 1.38757200435625,
```

## Manhattan plot

This function generates a “plotly” graph.

``` r
p <- manPlot(gwas = gwasRes,
             adj_method = "bonferroni",
             thresh_p = 0.05,
             chr = NA,
             title = "Readme Example")
#> 2021-08-25 15:07:01 - r-manPlot(): Check chromosome name ...
#> 2021-08-25 15:07:01 - r-manPlot(): Check chromosome name DONE
#> 2021-08-25 15:07:01 - r-manPlot(): Adjust p-values ...
#> 2021-08-25 15:07:01 - r-adjustPval(): Check adj_method ...
#> 2021-08-25 15:07:01 - r-adjustPval(): Check adj_method DONE
#> 2021-08-25 15:07:01 - r-adjustPval(): Check p values ...
#> 2021-08-25 15:07:01 - r-adjustPval(): Check p values DONE
#> 2021-08-25 15:07:01 - r-adjustPval(): Adjust p-values ...
#> 2021-08-25 15:07:01 - r-adjustPval(): Adjust p-values DONE
#> 2021-08-25 15:07:01 - r-adjustPval(): Adjust threshold ...
#> 2021-08-25 15:07:01 - r-adjustPval(): Adjust threshold DONE
#> 2021-08-25 15:07:01 - r-adjustPval(): DONE, return output
#> 2021-08-25 15:07:01 - r-manPlot(): Adjust p-values DONE
#> 2021-08-25 15:07:02 - r-manPlot(): DONE, return output
```

![screenshot of the plotly graph](README_files/manPlot2.png)

## LD PLot

Compute \(r^2\) Linkage Disequilibrium (LD) between given SNPs and save
a plot in a temporary PNG file.

``` r
gDta <- readGenoData("data/geno/testMarkerData01.vcf.gz")
#> 2021-08-25 15:07:02 - r-readGenoData(): Check file extention ... 
#> 2021-08-25 15:07:02 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-08-25 15:07:02 - r-readGenoData(): Read geno file DONE 
#> 2021-08-25 15:07:02 - r-readGenoData(): DONE, return output.
(imgFile <- LDplot(geno = gDta,
                   from = 1,
                   to = 20,
                   file = tempfile(fileext = ".png")))
#> 2021-08-25 15:07:02 - r-LDplot(): Check "from" and "to" format ...
#> 2021-08-25 15:07:02 - r-LDplot(): Check "from" and "to" format DONE
#> 2021-08-25 15:07:02 - r-LDplot(): Check "from" < "to"...
#> 2021-08-25 15:07:02 - r-LDplot(): Check "from" < "to" DONE
#> 2021-08-25 15:07:02 - r-LDplot(): Check number of SNP < 50...
#> 2021-08-25 15:07:02 - r-LDplot(): Check number of SNP < 50 DONE
#> 2021-08-25 15:07:02 - r-LDplot(): Check file ...
#> 2021-08-25 15:07:02 - r-LDplot(): Check file DONE
#> 2021-08-25 15:07:02 - r-LDplot(): Compute LD ...
#> 2021-08-25 15:07:02 - r-LDplot(): Compute LD DONE
#> 2021-08-25 15:07:02 - r-LDplot(): Create LD plot ...
#> 2021-08-25 15:07:02 - r-LDplot(): Create create file: /tmp/RtmpAxU90K/file3b7124b2c31d0.png
#> 2021-08-25 15:07:03 - r-LDplot(): Create LD plot DONE
#> 2021-08-25 15:07:03 - r-LDplot(): DONE, return output
#> [1] "/tmp/RtmpAxU90K/file3b7124b2c31d0.png"
```

``` r
img <- readPNG(imgFile)
grid::grid.raster(img)
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

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

# Data references

The genotypic and phenotypic data used as example come from:

Keyan Zhao, Chih-Wei Tung, Georgia C. Eizenga, Mark H. Wright, M. Liakat
Ali, Adam H. Price, Gareth J. Norton, M. Rafiqul Islam, Andy Reynolds,
Jason Mezey, Anna M. McClung, Carlos D. Bustamante & Susan R. McCouch
(2011). [Genome-wide association mapping reveals a rich genetic
architecture of complex traits in *Oryza
sativa*.](http://www.nature.com/ncomms/journal/v2/n9/full/ncomms1467.html)
Nat Comm 2:467 | DOI: 10.1038/ncomms1467, Published Online 13 Sep 2011.
