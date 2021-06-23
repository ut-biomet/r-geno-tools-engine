
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

## Load data

Example data are stored in the `data` folder.

``` r
gDta <- readGenoData("data/markers/testMarkerData01.vcf.gz")
#> 2021-06-23 11:16:58 - r-readGenoData(): Check file extention ... 
#> 2021-06-23 11:16:58 - r-readGenoData(): Read geno file ... 
#> ped stats and snps stats have been set. 
#> 'p' has been set. 
#> 'mu' and 'sigma' have been set.
#> 2021-06-23 11:17:00 - r-readGenoData(): Read geno file DONE 
#> 2021-06-23 11:17:00 - r-readGenoData(): DONE, return output.
pDta <- readPhenoData("data/pheno/testPhenoData01.csv")
#> 2021-06-23 11:17:00 - r-readPhenoData(): Read phenotypic file ... 
#> 2021-06-23 11:17:00 - r-readPhenoData(): Read phenotypic file DONE 
#> 2021-06-23 11:17:00 - r-readPhenoData(): DONE, return output.
```

For online data, they can be loaded using:

``` r
gDta <- downloadGenoData("url/to/geno_data.vcf")
pDta <- downloadPhenoData("url/to/pheno_data.csv")
dta <- downloadData("url/to/geno_data.vcf", "url/to/pheno_data.csv")
```

## GWAS

To run a GWAS analysis you should first process the data using
`prepareData` and then call the `gwas` function.

The `gwas` function return a `data.frame`.

``` r
data <- prepareData(gDta=gDta, pDta = pDta)
#> 2021-06-23 11:17:00 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set ...
#> 2021-06-23 11:17:00 - r-prepareData(): Remove from geno data individuals that are not in phenotypic data-set DONE
#> 2021-06-23 11:17:00 - r-prepareData(): reorder matrix ...
#> 2021-06-23 11:17:00 - r-prepareData(): reorder matrix DONE
#> 2021-06-23 11:17:00 - r-prepareData(): remove monomorphic geno ...
#> 2021-06-23 11:17:01 - r-prepareData(): remove monomorphic geno DONE
#> 2021-06-23 11:17:01 - r-prepareData(): calculate genetic relatinoal matrix ...
#> 2021-06-23 11:17:01 - r-prepareData(): calculate genetic relatinoal matrix DONE
#> 2021-06-23 11:17:01 - r-prepareData(): DONE, return output.

gwasRes <- gwas(data = data,
                trait = "Panicle.length",
                test = "score",
                fixed = 0,
                response = "quantitative",
                thresh_maf = 0.05,
                thresh_callrate = 0.95)
#> 2021-06-23 11:17:01 - r-gwas(): Check inputs ...
#> 2021-06-23 11:17:01 - r-gwas(): Check inputs DONE
#> 2021-06-23 11:17:01 - r-gwas(): aggregate data in bed matrix ...
#> 2021-06-23 11:17:01 - r-gwas(): aggregate data in bed matrix DONE
#> 2021-06-23 11:17:01 - r-gwas(): remove samples with missing phenotypic values ...
#> 2021-06-23 11:17:02 - r-gwas(): remove samples with missing phenotypic values DONE
#> 2021-06-23 11:17:02 - r-gwas(): filter SNPs ...
#> 2021-06-23 11:17:02 - r-gwas(): filter SNPs DONE
#> 2021-06-23 11:17:02 - r-gwas(): fit model ...
#> [Iteration 1] theta = 6.11889 3.12119
#> [Iteration 1] log L = -545.19
#> [Iteration 1] AI-REML update
#> [Iteration 1] ||gradient|| = 4.44998
#> [Iteration 2] theta = 2.40544 5.61424
#> [Iteration 2] log L = -536.314
#> [Iteration 2] AI-REML update
#> [Iteration 2] ||gradient|| = 5.68174
#> [Iteration 3] theta = 2.61763 6.58782
#> [Iteration 3] log L = -534.365
#> [Iteration 3] AI-REML update
#> [Iteration 3] ||gradient|| = 0.538832
#> [Iteration 4] theta = 2.60703 6.83004
#> [Iteration 4] log L = -534.325
#> [Iteration 4] AI-REML update
#> [Iteration 4] ||gradient|| = 0.0362738
#> [Iteration 5] theta = 2.59179 6.86789
#> [Iteration 5] log L = -534.324
#> [Iteration 5] AI-REML update
#> [Iteration 5] ||gradient|| = 0.0110441
#> [Iteration 6] theta = 2.58774 6.87671
#> [Iteration 6] log L = -534.324
#> [Iteration 6] AI-REML update
#> [Iteration 6] ||gradient|| = 0.00279312
#> [Iteration 7] theta = 2.58672  6.8789
#> [Iteration 7] log L = -534.324
#> [Iteration 7] AI-REML update
#> [Iteration 7] ||gradient|| = 0.000699973
#> [Iteration 8] theta = 2.58647 6.87945
#> [Iteration 8] log L = -534.324
#> [Iteration 8] AI-REML update
#> [Iteration 8] ||gradient|| = 0.000175033
#> [Iteration 9] theta = 2.58641 6.87958
#> [Iteration 9] log L = -534.324
#> [Iteration 9] AI-REML update
#> [Iteration 9] ||gradient|| = 4.37443e-05
#> [Iteration 10] theta = 2.58639 6.87962
#> [Iteration 10] log L = -534.324
#> [Iteration 10] AI-REML update
#> [Iteration 10] ||gradient|| = 1.09311e-05
#> [Iteration 11] theta = 2.58639 6.87963
#> [Iteration 11] log L = -534.324
#> [Iteration 11] AI-REML update
#> [Iteration 11] ||gradient|| = 2.73144e-06
#> 2021-06-23 11:17:02 - r-gwas(): fit model DONE
#> 2021-06-23 11:17:02 - r-gwas(): DONE, return output.
head(gwasRes)
#>   chr   pos           id A1 A2    freqA2      score         p
#> 1   1  9563  SNP-1.8562.  A  T 0.1382682 1.45576079 0.2276058
#> 2   1 25922 SNP-1.24921.  C  T 0.1278090 0.36310219 0.5467884
#> 3   1 26254 SNP-1.25253.  A  T 0.2955182 0.02865155 0.8655861
#> 4   1 30214 SNP-1.29213.  T  A 0.1344538 0.51075282 0.4748131
#> 5   1 31478 SNP-1.30477.  C  T 0.2262570 0.21535690 0.6426007
#> 6   1 32733 SNP-1.31732.  T  G 0.3011204 0.05018114 0.8227484
```

Save GWAS result in a `.json` file.

``` r
(file <- saveGWAS(gwasRes))
#> [1] "/tmp/Rtmpb6ocsA/file1f66943907fd1.json"
cat(paste(readLines(file, n=20), collapse = "\n"))
#> [
#>   {
#>     "chr": "1",
#>     "pos": 9563,
#>     "id": "SNP-1.8562.",
#>     "A1": "A",
#>     "A2": "T",
#>     "freqA2": 0.1383,
#>     "score": 1.4558,
#>     "p": 0.2276
#>   },
#>   {
#>     "chr": "1",
#>     "pos": 25922,
#>     "id": "SNP-1.24921.",
#>     "A1": "C",
#>     "A2": "T",
#>     "freqA2": 0.1278,
#>     "score": 0.3631,
#>     "p": 0.5468
```

## Manhattan plot

This function generates a “plotly” graph.

``` r
manPlot(gwas = gwasRes,
        adj_method = "bonferroni",
        thresh_p = 0.05,
        chr = NA,
        title = "Readme Example")
```

![screenshot of the plotly graph](README_files/manPlot.png)

## LD PLot

Compute *r*<sup>2</sup> Linkage Disequilibrium (LD) between given SNPs
and save a plot in a temporary PNG file.

``` r
(imgFile <- LDplot(gDta, 1, 20, write = TRUE))
#> 2021-06-23 11:17:03 - r-LDplot(): Check "from" < "to"...
#> 2021-06-23 11:17:03 - r-LDplot(): Check "from" < "to" DONE
#> 2021-06-23 11:17:03 - r-LDplot(): Check number of SNP < 50...
#> 2021-06-23 11:17:03 - r-LDplot(): Check number of SNP < 50 DONE
#> 2021-06-23 11:17:03 - r-LDplot(): Compute LD ...
#> 2021-06-23 11:17:03 - r-LDplot(): Compute LD DONE
#> 2021-06-23 11:17:03 - r-LDplot(): Create LD plot ...
#> 2021-06-23 11:17:03 - r-LDplot(): Create create file: /tmp/Rtmpb6ocsA/file1f66924ef4de.png
#> 2021-06-23 11:17:03 - r-LDplot(): Create LD plot DONE
#> 2021-06-23 11:17:03 - r-LDplot(): DONE, return output
#> [1] "/tmp/Rtmpb6ocsA/file1f66924ef4de.png"
```

``` r
library(png)
img <- readPNG(imgFile)
grid::grid.raster(img)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
