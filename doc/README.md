# `gwas`

perform GWAS analysis


## Description

perform GWAS analysis


## Usage

```r
gwas(
  data,
  trait,
  test,
  fixed = 0,
  response = "quantitative",
  thresh_maf,
  thresh_callrate
)
```


## Arguments

Argument      |Description
------------- |----------------
`data`     |     List return by `prepareDta` function
`trait`     |     Chraracter of length 1, name of the trait to analyze. Could be a column name of the phenotypic file
`test`     |     Which test to use. Either `"score"`,  `"wald"` or `"lrt"`. For binary phenotypes, test = `"score"` is mandatory. For more information about this parameters see: `??gaston::association.test`
`fixed`     |     Number of Principal Components to include in the model with fixed effect (for test = `"wald"` or `"lrt"`). Default value is 0. For more information about this parameters see: `??gaston::association.test`
`response`     |     Character of length 1, Either "quantitative" or "binary". Is the trait a quantitative or a binary phenotype? Default value is "quantitative"
`thresh_maf`     |     Threshold for filtering markers. Only markers with minor allele frequency > `thresh_maf` will be kept.
`thresh_callrate`     |     Threshold for filtering markers. Only markers with a callrate > `thresh_callrate` will be kept.


## Details

For the calculation, the genetic relationship matrix need to be calculated. This is done based on the genetic matrix standardized by the the genetic mean "mu" and the genetic variance "sigma", after having filtering according to the `thresh_maf` and `thresh_callrate`.


## Value

`data.frame`


# `adjustPval`

Adjust P-values for Multiple Comparisons


## Description

Adjust P-values for Multiple Comparisons


## Usage

```r
adjustPval(p, adj_method, thresh_p = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`adj_method`     |     correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
`thresh_p`     |     optional value of the p value significant threshold (default NULL)
`vector`     |     of p-values


## Details

The method "hommel" is not implemented because it is too long to calculate.


## Value

list of two elements: "p_adj" vector of adjusted p values, "thresh_adj" the adjusted threshold (if thresh_p is preovided, NULL if not)


# `run_gwas`

Run GWAS analysis


## Description

Run GWAS analysis


## Usage

```r
run_gwas(
  genoFile = NULL,
  phenoFile = NULL,
  genoUrl = NULL,
  phenoUrl = NULL,
  trait,
  test,
  fixed = 0,
  response = "quantitative",
  thresh_maf,
  thresh_callrate,
  outFile = tempfile(fileext = ".json")
)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     path of the geno data file (`.vcf` or `.vcf.gz` file)
`phenoFile`     |     path of the phenotypic data file (`csv` file)
`genoUrl`     |     url of the geno data file (`.vcf` or `.vcf.gz` file)
`phenoUrl`     |     url of the phenotypic data file (`csv` file)
`trait`     |     Chraracter of length 1, name of the trait to analyze. Must be a column name of the phenotypic file.
`test`     |     Which test to use. Either `"score"`,  `"wald"` or `"lrt"`. For binary phenotypes, test = `"score"` is mandatory. For more information about this parameters see: `??gaston::association.test`
`fixed`     |     Number of Principal Components to include in the model with fixed effect (for test = `"wald"` or `"lrt"`). Default value is 0. For more information about this parameters see: `??gaston::association.test`
`response`     |     Character of length 1, Either "quantitative" or "binary". Is the trait a quantitative or a binary phenotype? Default value is "quantitative"
`thresh_maf`     |     Threshold for filtering markers. Only markers with minor allele frequency > `thresh_maf` will be kept.
`thresh_callrate`     |     Threshold for filtering markers. Only markers with a callrate > `thresh_callrate` will be kept.
`outFile`     |     path of the output file. If `NULL`, the output will not be written in any file. By default write in an tempoary `.json` file.


## Value

list with 3 elements `gwasRes` for the results of the gwas analysis in json, `metadata` a list of metadata of these analysis and `file` path of the json file containing the results (id `dir` is not `NULL`)


# `draw_manhattanPlot`

Draw a Manhattan Plot


## Description

Draw a Manhattan Plot


## Usage

```r
draw_manhattanPlot(
  gwasFile = NULL,
  gwasUrl = NULL,
  adj_method = "bonferroni",
  thresh_p = 0.05,
  chr = NA,
  interactive = TRUE,
  filter_pAdj = 1,
  filter_nPoints = Inf,
  filter_quant = 1,
  outFile = tempfile()
)
```


## Arguments

Argument      |Description
------------- |----------------
`gwasFile`     |     path of the gwas result data file (json file)
`gwasUrl`     |     url of the gwas result data file (json file)
`adj_method`     |     correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
`thresh_p`     |     p value significant threshold (default 0.05)
`chr`     |     name of the chromosome to show (show all if NA)
`interactive`     |     [bool] should the plot be interactive (the default)
`filter_pAdj`     |     [numeric] threshold to remove points with pAdj > filter_pAdj from the plot (default no filtering)
`filter_nPoints`     |     [numeric] threshold to keep only the filter_nPoints with the lowest p-values for the plot (default no filtering)
`filter_quant`     |     [numeric] threshold to keep only the filter_quant*100 %  of the points with the lowest p-values for the plot (default no filtering)
`outFile`     |     path of the file containing the plot. If `NULL`, the output will not be written in any file. By default write in an tempoary file.


## Details

If several filtering rules are given, the filtering process apply
 the filtering process sequentially (this lead to having the same result
 that if only the strongest rules were given).
 Moreover, the number of points kept for the plot will be display
 in the plot title.


## Value

plotly graph if interactive is TRUE, or NULL if not.


# `draw_ldPlot`

Draw an LD Plot


## Description

Draw an LD Plot


## Usage

```r
draw_ldPlot(
  genoFile = NULL,
  genoUrl = NULL,
  from,
  to,
  outFile = tempfile(fileext = ".png")
)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     path of the geno data file (`.vcf` or `.vcf.gz` file)
`from`     |     lower bound of the range of SNPs for which the LD is computed
`to`     |     upper bound of the range of SNPs for which the LD is computed
`outFile`     |     path of the png file to save the plot. If `NULL`, the image file will not be created. By default write in an tempoary `.png` file.
`phenoFile`     |     path of the phenotypic data file (`csv` file)


## Value

path of the created file (or NULL if `file` is NULL)


# `run_resAdjustment`

Adjust GWAS p-values


## Description

Adjust GWAS p-values


## Usage

```r
run_resAdjustment(
  gwasFile = NULL,
  gwasUrl = NULL,
  adj_method = "bonferroni",
  outFile = tempfile(fileext = ".json")
)
```


## Arguments

Argument      |Description
------------- |----------------
`gwasFile`     |     path of the gwas result data file (json file)
`gwasUrl`     |     url of the gwas result data file (json file)
`adj_method`     |     correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
`outFile`     |     path of the output file. If `NULL`, the output will not be written in any file. By default write in an tempoary `.json` file.


## Value

list with 3 elements `gwasAdjusted` for the results of the gwas analysis in json with adjusted p-values, `metadata` a list of metadata of the gwas analysis in json with adjusted p-values, and `file` path of the json file containing the results (if `dir` is not `NULL`)


# `manPlot`

create manhatan plot


## Description

create manhatan plot


## Usage

```r
manPlot(
  gwas,
  adj_method,
  thresh_p = 0.05,
  chr = NA,
  title = "Manhattan Plot",
  filter_pAdj = 1,
  filter_nPoints = Inf,
  filter_quant = 1,
  interactive = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gwas`     |     [data.frame] output of the gwas function
`adj_method`     |     correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
`thresh_p`     |     p value significant threshold (default 0.05)
`chr`     |     [char] name of the chromosome to show (show all if NA)
`title`     |     [char] Title of the plot. Default is "Manhattan Plot"
`filter_pAdj`     |     [numeric] threshold to remove points with pAdj > filter_pAdj from the plot (default no filtering)
`filter_nPoints`     |     [numeric] threshold to keep only the filter_nPoints with the lowest p-values for the plot (default no filtering)
`filter_quant`     |     [numeric] threshold to keep only the filter_quant*100 %  of the points with the lowest p-values for the plot (default no filtering)
`interactive`     |     [bool] should the plot be interactive (the default) or not


## Details

If several filtering rules are given, the filtering process apply
 the filtering process sequentially (this lead to having the same result
 that if only the strongest rules were given).
 Moreover, the number of points kept for the plot will be display
 in the plot title.


## Value

plotly graph if interactive is TRUE, or NULL if not.


# `LDplot`

writeLDplot Compute r2 Linkage Disequilibrium (LD) between given SNPs and return a plot


## Description

writeLDplot Compute r2 Linkage Disequilibrium (LD) between given SNPs and return a plot


## Usage

```r
LDplot(geno, from, to, file = tempfile(fileext = ".png"))
```


## Arguments

Argument      |Description
------------- |----------------
`geno`     |     [bed.matrix] geno data return by function `readGenoData` or `downloadGenoData`.
`from`     |     lower bound of the range of SNPs for which the LD is computed
`to`     |     upper bound of the range of SNPs for which the LD is computed
`file`     |     path of the png file to save the plot. If `NULL`, the image file will not be created. By default write in an tempoary `.png` file.


## Details

`from` should be lower than `to`, and the maximum ranger size is 50.
 (In order to get a readable image). If write is `TRUE`, the function will write the plot in a png file, else it will plot it.


## Value

null if `dir` is NULL, else the path of the png file.


# `downloadGenoData`

Download geno data


## Description

Download geno data


## Usage

```r
downloadGenoData(url)
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the geno data file (.vcf.gz file)


## Value

`gaston::bed.matrix`


# `downloadPhenoData`

Download phenotypic data


## Description

Download phenotypic data


## Usage

```r
downloadPhenoData(url)
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the phenotypic data file (csv file)


## Value

`data.frame`


# `downloadData`

Download and prepare data for GWAS analysis


## Description

Download and prepare data for GWAS analysis


## Usage

```r
downloadData(genoUrl, phenoUrl)
```


## Arguments

Argument      |Description
------------- |----------------
`genoUrl`     |     url of the geno data file (.vcf.gz file)
`phenoUrl`     |     url of the phenotypic data file (csv file)


## Value

List


# `downloadGWAS`

Download a gwas reults


## Description

Download a gwas reults


## Usage

```r
downloadGWAS(url)
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the result data file (json file)


## Value

`data.frame`


# `readGenoData`

Read geno data from a file


## Description

Read geno data from a file


## Usage

```r
readGenoData(file)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     VCF file path (ext `.vcf` or `.vcf.gz`)


## Value

`gaston::bed.matrix`


# `readPhenoData`

Read phenotypic data file


## Description

Read phenotypic data file


## Usage

```r
readPhenoData(file, row.names = 1, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     file path
`row.names`     |     [default 1] a single number giving the column of the table which contains the row names
`...`     |     Further arguments to be passed to `read.csv`


## Value

`data.frame`


# `readData`

Read and prepare data for GWAS result


## Description

Read and prepare data for GWAS result


## Usage

```r
readData(genoFile, phenoFile)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     path of the geno data file (`.vcf` or `.vcf.gz` file)
`phenoFile`     |     path of the phenotypic data file (`csv` file)


## Value

List


# `readGWAS`

Read GWAS analysis result file (`.json`)


## Description

Read GWAS analysis result file (`.json`)


## Usage

```r
readGWAS(file)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the json file generated by the function `gwas` containing GWAS result


## Value

`list` of 2 elements `gwas` (data.frame) and `metadata` (list)


# `saveGWAS`

saveGWAS save gwas result in a temporary file


## Description

saveGWAS save gwas result in a temporary file


## Usage

```r
saveGWAS(gwasRes, metadata, dir = NULL, file = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`gwasRes`     |     data.frame return by `gwas` function
`metadata`     |     list of metadata of the gwas results
`dir`     |     if `filename` is NULL, directory where to save the data, by default it is a temporary directory
`file`     |     file path where to save the data. If the file already exists, it will be overwritten. Default NULL


## Value

path of the created filed


# `prepareData`

Filter individuals and remove monomorphic markers


## Description

Filter individuals and remove monomorphic markers


## Usage

```r
prepareData(gDta, pDta)
```


## Arguments

Argument      |Description
------------- |----------------
`gDta`     |     output of `downloadGenoData` or `readGenoData` functions
`pDta`     |     output of `downloadPhenoData` or `readPhenoData` functions


## Details

The function remove the monomorphic markers and


## Value

List of 2 elements: `genoData` (a bed matrix), `phenoData` (a data.frame)


# `logger`

R6 class use to log messages in this engine's function


## Description

R6 class use to log messages in this engine's function
 
 R6 class use to log messages in this engine's function


## Examples

```r
## ------------------------------------------------
## Method `logger$new`
## ------------------------------------------------

mylogger <- logger$new(context = NULL)
```


# `writeDoc`

Write engine documentation


## Description

Write engine documentation


## Usage

```r
writeDoc(srcDir = "./src", docDir = "./doc")
```


## Arguments

Argument      |Description
------------- |----------------
`srcDir`     |     path of R sources folder (default "./src")
`docDir`     |     path of documentation folder (default "./doc")


## Details

Write engine's functions documentation in a README.md file located in `docDir`.


