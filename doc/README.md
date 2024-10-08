# `checkAndFilterSNPcoord`

Check compatibility between 2 snps coordinates data set
 and keep only genotypes SNPs


## Description

If the `user_SNPcoord` data miss some SNPs defined in the `.vcf`
 file, an error is raised. If the `user_SNPcoord` data have additional SNPs,
 those SNPs will be removed.
 If the data between those two data-set are inconsistent, an error is raised.


## Usage

```r
checkAndFilterSNPcoord(user_SNPcoord, vcf_SNPcoord)
```


## Arguments

Argument      |Description
------------- |----------------
`user_SNPcoord`     |     SNPs coordinates data coming from the user (`.csv` file)
`vcf_SNPcoord`     |     SNPs coordinates data coming from the `.vcf` file


## Value

The filtered `user_SNPcoord` data.frame


# `checkIndNamesConsistency`

Check individuals in the crossing table are in the haplotype data


## Description

Check individuals in the crossing table are in the haplotype data


## Usage

```r
checkIndNamesConsistency(crossTable, haplo)
```


## Arguments

Argument      |Description
------------- |----------------
`crossTable`     |     the crossing table
`haplo`     |     the haplotype data given by the function `readPhasedGeno`


## Value

NULL, raise error if missing individuals are detected.


# `initializeSimulation`

Initialise the simulation


## Description

Initialise the simulation


## Usage

```r
initializeSimulation(haplotypes, SNPcoord)
```


## Arguments

Argument      |Description
------------- |----------------
`haplotypes`     |     haplotypes of the parents, data.frame with genotype values in row, and individuals'haplotype in columns. The columns name should be `individualName_1` and `individualName_2` for the first/second haplotype of the individual named `individualName`. (the list item named `haplo` return by the function `readPhasedGeno`)
`SNPcoord`     |     snp coordinates, data.frame of 4 columns: - `chr`: Name of the chromosome holding the SNP - `physPos`: SNP physical position on the chromosome - `linkMapPos`: SNP linkage map position on the chromosome **in morgan** - `SNPid`: SNP's IDs


## Value

`breedSimulatR`'s population object


# `engineError`

Main function to raise an `engineError` (ie. expected error)


## Description

This function is similar to the R's `stop` function and should be used
 instead to ensure we raise expected errors.


## Usage

```r
engineError(message, extra = list(), n_skip_caller = 1)
```


## Arguments

Argument      |Description
------------- |----------------
`message`     |     error message
`extra`     |     list of extra information that will be included in the error
`n_skip_caller`     |     (int, default 1) This is to catch the function where the error happens. 0 will show this function, 1 will show the function calling this one and so on.


# `bad_argument`

Helper function to raise an `engineError` for a "bad argument"


## Description

To be used inside functions to check their arguments.
 the error message will be:
 "`arg` must be `must_be` not `not`"
 (adapted from https://adv-r.hadley.nz/conditions.html#signalling)


## Usage

```r
bad_argument(
  arg,
  must_be,
  not = NULL,
  errType = "value",
  class = "engineError",
  extra = NULL,
  n_skip_caller = 2
)
```


## Arguments

Argument      |Description
------------- |----------------
`arg`     |     (character) tested argument name
`must_be`     |     (character) expected value or type
`not`     |     (any) provided value
`errType`     |     "type" is the only recognised value. Use it to report an error about the type of the argument
`class`     |     use "engineError" (default) to raise an expected "engineError"
`extra`     |     (list, default NULL) extra information to pass to engineError
`n_skip_caller`     |     (int, default 2) see `engineError`


## Examples

```r
x = "a"
bad_argument("x", must_be = 42, not = x)
# will stop with the error msg: "`x` must be 42 not `a`"
bad_argument("x", must_be = "numeric", not = x, "type")
# will stop with the error msg: "`x` must be numeric not `character`"
```


# `errorCode`

list of error codes as function so that unexpected error code raise and error


## Description

list of error codes as function so that unexpected error code raise and error


## Usage

```r
errorCode(code)
```


# `filterGWAS`

Filter gwas results


## Description

Filter gwas results


## Usage

```r
filterGWAS(gwas, filter_pAdj = 1, filter_nPoints = Inf, filter_quant = 1)
```


## Arguments

Argument      |Description
------------- |----------------
`gwas`     |     [data.frame] output of the gwas function
`filter_pAdj`     |     [numeric] threshold to remove points with pAdj < filter_pAdj from the plot (default no filtering)
`filter_nPoints`     |     [numeric] threshold to keep only the filter_nPoints with the lowest p-values for the plot (default no filtering)
`filter_quant`     |     [numeric] threshold to keep only the filter_quant*100 %  of the points with the lowest p-values for the plot (default no filtering)


# `fit_with_gaston`

Fit a GS model using `gaston` package


## Description

Fit a GS model using `gaston` package


## Usage

```r
fit_with_gaston(pheno, K)
```


## Arguments

Argument      |Description
------------- |----------------
`pheno`     |     **1 dimensional** vector of phenotypic values
`K`     |     list of 1 or 2 variance covariance matrices to consider for the random effect.


## Details

This is a wrapper around `gaston::lmm.aireml()` see
 `?gaston::lmm.aireml()` for more information


## Value

list of 3 elements
 - `logL`: Value of log-likelihood of the model (cf. `?gaston::lmm.aireml()`)
 - `blups`: `length(pheno)` by `length(K)` matrix of the blups of the
 individuals (one random effects per column).
 - `intercept`: value of the intercept of the model


# `fit_with_rainbowr`

Fit a GS model using `RAINBOWR` package


## Description

Fit a GS model using `RAINBOWR` package


## Usage

```r
fit_with_rainbowr(pheno, K)
```


## Arguments

Argument      |Description
------------- |----------------
`pheno`     |     **1 dimensional** vector of phenotypic values
`K`     |     list of 1 or 2 variance covariance matrices used for the random effect.


## Details

This is a wrapper around `RAINBOWR::EM3.general()` see
 `?RAINBOWR::EM3.general()` for more information


## Value

list of 3 elements
 - `logL`: Value of log-likelihood of the model (cf. `?RAINBOWR::EM3.general()`)
 - `blups`: `length(pheno)` by `length(K)` matrix of the blups of the
 individuals (one random effects per column).
 - `intercept`: value of the intercept of the model


# `calc_mark_eff`

Calculate markers effects from the BLUPS, the genetic matrix and the
 variance covariance matrix of the blups (ie. genetic relationship matrices)


## Description

Calculate markers effects from the BLUPS, the genetic matrix and the
 variance covariance matrix of the blups (ie. genetic relationship matrices)


## Usage

```r
calc_mark_eff(geno_mat, rel_mat, blups)
```


## Arguments

Argument      |Description
------------- |----------------
`geno_mat`     |     `n_inds` by `n_marker` genomic matrix
`rel_mat`     |     `n_inds` by `n_inds` relationship matrix associated to `geno_mat`
`blups`     |     `n_inds` vector of the individuals blups effect.


## Details

Calculation is made by
 
$$(\frac{X}{n_{mark}})(\frac{X}{n_{mark}})'Z^{-1}B$$
 
 
 with:
 - $X$ the genetic matrix
 - $n_ mark $ the number of markers (ie. number of columns of $X$)
 - $Z$ the relationship matrix of $X$
 - $B$ the blups of the individuals


## Value

vector of the estimated markers effects


# `train_gs_model`

Train GS model optionally with dominance


## Description

Train GS model optionally with dominance


## Usage

```r
train_gs_model(pheno, geno, with_dominance = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
`pheno`     |     1 column data.frame of the phenotype
`geno`     |     bed.matrix returned by `readGenoData()`
`with_dominance`     |     (default = FALSE) control if the model should include dominance effects


## Details

IMPORTANT ! It will return the markers effects for the genotypes encoded:
 - **in allele dose (0, 1, 2) for the additive effects**
 - **as (0, 1, 0) for the dominance effects**


## Value

list of 2 elements:
 - `intercept`: the intercept of the model
 - `eff`: data.frame of 2 columns `additive` and `dominance` with
 their corresponding markers effects


# `predict_gs_model`

Make GS predictions


## Description

Make GS predictions


## Usage

```r
predict_gs_model(geno, estim_mark_eff)
```


## Arguments

Argument      |Description
------------- |----------------
`geno`     |     bed.matrix returned by `readGenoData()` on wich we want to make prediction
`estim_mark_eff`     |     (list, output of `train_gs_model()`) estimated markers effects and intercept.


## Value

data.frame of one column containing the prediction


# `cross_validation_evaluation`

Repeated K-Folds cross-validation of a GS model


## Description

Repeated K-Folds cross-validation of a GS model


## Usage

```r
cross_validation_evaluation(
  pheno,
  geno,
  with_dominance = TRUE,
  n_folds = 10,
  n_repetitions = 5
)
```


## Arguments

Argument      |Description
------------- |----------------
`pheno`     |     1 column data.frame of the phenotype
`geno`     |     bed.matrix returned by `readGenoData()`
`with_dominance`     |     (default = FALSE) control if the model should include dominance effects
`n_folds`     |     (default `10`) number of folds for each repetition
`n_repetitions`     |     (default `5`) number of repetition


## Value

list of 2 elements:
 - `predictions`: data.frame of the predicted values during the
 cross-validation. Available columns:
 - `ind` individual id
 - `actual` phenotype value in the training data
 - `predicted` predicted values
 - `fold` cross-validation fold id
 - `repetition` cross-validation repetition id
 - `metrics`: data.frame of the calculated model metrics during the
 cross-validation. Available columns are:
 - `fold` cross-validation fold id
 - `repetition` cross-validation repetition id
 - `...` values of the metrics returned by `get_model_metrics()` (one column
 per returned list's elements with the same name)


# `get_model_metrics`

Calculate model metrics


## Description

Calculate model metrics


## Usage

```r
get_model_metrics(actual, predictions)
```


## Arguments

Argument      |Description
------------- |----------------
`actual`     |     vector of actual values
`predictions`     |     vector of model's predicted values


## Value

list of metrics values:
 - `rmse` root mean square error
 - `corel_pearson` Pearson's corelation coefficient
 - `corel_spearman` Spearman's corelation coefficient
 - `r2` R squared


# `check_dominance_model_is_applicable`

Check if a dominance model is applicable with the provided genetic data


## Description

Check if a dominance model is applicable with the provided genetic data


## Usage

```r
check_dominance_model_is_applicable(dominance, homozygous_threshold = 0.95)
```


## Arguments

Argument      |Description
------------- |----------------
`dominance`     |     list returned by
`homozygous_threshold`     |     A threshold used to identify individuals or SNPs with a homozygosity proportion exceeding this value. These will be counted and reported in the error message if applicable. This parameter does not influence the function's behavior, only the error message it can raise.


## Details

The dominance model need the dominance relationship matrix to be invertible
 in order to be able to calculate the dominance effects. This function will
 return an error if it is not the case.
 The error will contain information about the homozygousity of the markers and
 individuals as if the data have too many homozygous individuals/markers it
 will probably not be suited for a dominance model. The ``


## Value

`TRUE` or raise an `engineError`


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


# `Logger`

R6 class use to log messages in this engine's function


## Description

R6 class use to log messages in this engine's function
 
 R6 class use to log messages in this engine's function


## Examples

```r
## ------------------------------------------------
## Method `Logger$new`
## ------------------------------------------------

mylogger <- Logger$new(context = NULL)
```


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
  fixed = NULL,
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
`phenoFile`     |     path of the phenotypic data file (`csv` file). Individuals' name should be the first column of the file and no duplication is allowed.
`genoUrl`     |     url of the geno data file (`.vcf` or `.vcf.gz` file)
`phenoUrl`     |     url of the phenotypic data file (`csv` file) Individuals' name should be the first column of the file and no duplication is allowed.
`trait`     |     Chraracter of length 1, name of the trait to analyze. Must be a column name of the phenotypic file.
`test`     |     Which test to use. Either `"score"`,  `"wald"` or `"lrt"`. For binary phenotypes, test = `"score"` is mandatory. For more information about this parameters see: `??gaston::association.test`
`fixed`     |     Number of Principal Components to include in the model with fixed effect (for test = `"wald"` or `"lrt"`). Default value is 0. For more information about this parameters see: `??gaston::association.test`
`response`     |     Character of length 1, Either "quantitative" or "binary". Is the trait a quantitative or a binary phenotype? Default value is "quantitative"
`thresh_maf`     |     Threshold for filtering markers. Only markers with minor allele frequency > `thresh_maf` will be kept.
`thresh_callrate`     |     Threshold for filtering markers. Only markers with a callrate > `thresh_callrate` will be kept.
`outFile`     |     path of the output file. If `NULL`, the output will not be written in any file. By default write in an tempoary `.json` file.


## Value

list with 3 elements `gwasRes` for the results of the gwas analysis in json, `metadata` a list of metadata of these analysis and `file` path of the json file containing the results


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
`filter_pAdj`     |     [numeric] threshold to remove points with pAdj < filter_pAdj from the plot (default no filtering)
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
  filter_pAdj = 1,
  filter_nPoints = Inf,
  filter_quant = 1,
  outFile = tempfile(fileext = ".json")
)
```


## Arguments

Argument      |Description
------------- |----------------
`gwasFile`     |     path of the gwas result data file (json file)
`gwasUrl`     |     url of the gwas result data file (json file)
`adj_method`     |     correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
`filter_pAdj`     |     [numeric] threshold to remove points with pAdj < filter_pAdj from the plot (default no filtering)
`filter_nPoints`     |     [numeric] threshold to keep only the filter_nPoints with the lowest p-values for the plot (default no filtering)
`filter_quant`     |     [numeric] threshold to keep only the filter_quant*100 %  of the points with the lowest p-values for the plot (default no filtering)
`outFile`     |     path of the output file. If `NULL`, the output will not be written in any file. By default write in an tempoary `.json` file.


## Value

list with 3 elements `gwasAdjusted` for the results of the gwas analysis in json with adjusted p-values, `metadata` a list of metadata of the gwas analysis in json with adjusted p-values, and `file` path of the json file containing the results (if `dir` is not `NULL`)


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
  n_max = 50,
  outFile = tempfile(fileext = ".png")
)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     path of the geno data file (`.vcf` or `.vcf.gz` file)
`from`     |     lower bound of the range of SNPs for which the LD is computed
`to`     |     upper bound of the range of SNPs for which the LD is computed
`n_max`     |     maximum number of marker to show (to avoid unreadable plot) created. By default write in an tempoary `.png` file.
`outFile`     |     path of the png file to save the plot. If `NULL`, the image file will not be


## Value

path of the created file (or NULL if `file` is NULL)


# `calc_pedRelMat`

Calculate pedigree relationship matrix


## Description

Calculate pedigree relationship matrix


## Usage

```r
calc_pedRelMat(
  pedFile = NULL,
  pedUrl = NULL,
  unknown_string = "",
  header = TRUE,
  outFile = tempfile(fileext = ".csv"),
  outFormat = tools::file_ext(outFile)
)
```


## Arguments

Argument      |Description
------------- |----------------
`pedFile`     |     path of the pedigree data file (`csv` file).
`pedUrl`     |     url of the pedigree data file (`csv` file).
`unknown_string`     |     [default: ""] a character vector of strings which are to be interpreted as "unknown parent". By default: missing value in the file.
`header`     |     [default: TRUE] a logical value indicating whether the file contains the names of the variables as its first line. The default value is TRUE. In any cases, the column 1 will be interpreted as the individual id, column 2 as the first parent, column 3 as the second parent.
`outFile`     |     path of the output file. If `NULL`, the output will not be written in any file. By default write in an tempoary `.json` file.
`outFormat`     |     Format of the output file, either `csv` or `json`. by default it will use the file extension of `outfile`.


## Details

For `csv` output, the file will include some metadata lines (starting by a `#`
 symbol), a header and the row ids in its first column.


## Value

list with 3 elements `relMat` the relationship matrix, `metadata` a
 list of metadata of these analysis (pedigree fingerprint,
 number of individuals, creation time) and `file` path
 of the file containing the results.


# `calc_genoRelMat`

Calculate genomic relationship matrix


## Description

Calculate genomic relationship matrix


## Usage

```r
calc_genoRelMat(
  genoFile = NULL,
  genoUrl = NULL,
  outFile = tempfile(fileext = ".csv"),
  outFormat = tools::file_ext(outFile)
)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     path of the geno data file (`.vcf` or `.vcf.gz` file)
`genoUrl`     |     url of the geno data file (`.vcf` or `.vcf.gz` file)
`outFile`     |     path of the output file. If `NULL`, the output will not be written in any file. By default write in an tempoary `.json` file.
`outFormat`     |     Format of the output file, either `csv` or `json`. by default it will use the file extension of `outfile`.


## Details

For `csv` output, the file will include some metadata lines (starting by a `#`
 symbol), a header and the row ids in its first column.


## Value

list with 3 elements `relMat` the relationship matrix, `metadata` a
 list of metadata of these analysis (pedigree fingerprint,
 number of individuals, creation time) and `file` path
 of the file containing the results.


# `calc_combinedRelMat`

Combined (pedigree + genomic) Relationship Matrix


## Description

Correct a pedigree relationship matrix by using genomic relationship matrix.


## Usage

```r
calc_combinedRelMat(
  pedRelMatFile = NULL,
  pedRelMatUrl = NULL,
  genoRelMatFile = NULL,
  genoRelMatUrl = NULL,
  method = "Legarra",
  tau = NULL,
  omega = NULL,
  outFile = tempfile(fileext = ".csv"),
  outFormat = tools::file_ext(outFile)
)
```


## Arguments

Argument      |Description
------------- |----------------
`pedRelMatFile`     |     path of a pedigree relationship matrix generated by the the engine.
`pedRelMatUrl`     |     url of a pedigree relationship matrix generated by the the engine.
`genoRelMatFile`     |     path of a genomic relationship matrix generated by the the engine.
`genoRelMatUrl`     |     url of a genomic relationship matrix generated by the the engine.
`method`     |     method to use, either "Legarra" or "Martini"
`tau`     |     tau parameter of the Martini's method
`omega`     |     omega parameter of the Martini's method
`outFile`     |     path of the output file. If `NULL`, the output will not be written in any file. By default write in an tempoary `.json` file.
`outFormat`     |     Format of the output file, either `csv` or `json`. by default it will use the file extension of `outfile`.


## Details

This method correct the pedigree matrix with the genomic relationship matrix.
 Therefore, individuals in the genomic relationship matrix not in the pedigree
 relationship matrix will be ignored.
 
 Using the Martini's method with `tau=1`, and `omega=1` is equivalent of
 Legarra's method.
 
 For `csv` output, the file will include some metadata lines (starting by a `#`
 symbol), a header and the row ids in its first column.


## Value

list with 3 elements `relMat` the relationship matrix, `metadata` a
 list of metadata of these analysis (pedigree fingerprint,
 number of individuals, creation time) and `file` path
 of the file containing the results.


## References

*Martini, JW, et al. 2018 The effect of the H-1 scaling factors tau and omega on the structure of H in the single-step procedure. Genetics Selection Evolution 50(1), 16*
 
 *Legarra, A, et al. 2009 A relationship matrix including full pedigree and genomic information. Journal of Dairy Science 92, 4656–4663*


# `draw_relHeatmap`

Draw a heatmap of a relationship matrix


## Description

Draw a heatmap of a relationship matrix


## Usage

```r
draw_relHeatmap(
  relMatFile = NULL,
  relMatUrl = NULL,
  format = NULL,
  interactive = TRUE,
  outFile = tempfile()
)
```


## Arguments

Argument      |Description
------------- |----------------
`relMatFile`     |     path of a file generated by the function `saveRelMat`
`relMatUrl`     |     url of a file generated by the function `saveRelMat`
`interactive`     |     [bool] should the plot be interactive (the default) or not
`outFile`     |     path of the file containing the plot. If `NULL`, the output will not be written in any file. By default write in an tempoary file.


## Value

plotly graph if interactive is TRUE, or NULL if not.


# `draw_pedNetwork`

Draw interactive pedigree network


## Description

Draw interactive pedigree network


## Usage

```r
draw_pedNetwork(
  pedFile = NULL,
  pedUrl = NULL,
  unknown_string = "",
  header = TRUE,
  outFile = tempfile(fileext = ".html")
)
```


## Arguments

Argument      |Description
------------- |----------------
`pedFile`     |     path of the pedigree data file (`csv` file).
`pedUrl`     |     url of the pedigree data file (`csv` file).
`unknown_string`     |     [default: ""] a character vector of strings which are to be interpreted as "unknown parent". By default: missing value in the file.
`header`     |     [default: TRUE] a logical value indicating whether the file contains the names of the variables as its first line. The default value is TRUE. In any cases, the column 1 will be interpreted as the individual id, column 2 as the first parent, column 3 as the second parent.
`outFile`     |     path of the file containing the plot. If `NULL`, the output will not be written in any file. By default write in an temporary file.


## Value

plotly graph if interactive is TRUE, or NULL if not.


# `crossingSimulation`

Simulate the genotypes of offspring given the parent genotypes.


## Description

Simulate the genotypes of offspring given the parent genotypes.


## Usage

```r
crossingSimulation(
  genoFile = NULL,
  genoUrl = NULL,
  crossTableFile = NULL,
  crossTableUrl = NULL,
  SNPcoordFile = NULL,
  SNPcoordUrl = NULL,
  nCross = 30,
  outFile = tempfile(fileext = ".vcf.gz")
)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     phased VCF file path (ext `.vcf` or `.vcf.gz`)
`genoUrl`     |     url of a phased VCF file path (ext `.vcf` or `.vcf.gz`)
`crossTableFile`     |     path of the crossing table data file (`csv` file of 2 or 3 columns). It must contain the names of the variables as its first line. The column 1 and 2 will be interpreted as the parents ids. The optional third column will be interpreted as the offspring base name.
`crossTableUrl`     |     URL of a crossing table file
`SNPcoordFile`     |     path of the SNPs coordinates file (`csv` file). This `.csv` file should have 4 named columns: - `chr`: Chromosome holding the SNP (mandatory) - `physPos`: SNP physical position on the chromosome - `linkMapPos`: SNP linkage map position on the chromosome in Morgan (mandatory) - `SNPid`: SNP's IDs
`SNPcoordUrl`     |     URL of a SNP coordinate file
`nCross`     |     number of cross to simulate for each parent pair defined in the crossing table.
`outFile`     |     path of the `.vcf.gz` file containing the simulated genotypes of the offspring. It must end by `.vcf.gz`. By default write in an temporary file.  For SNPcoordFile/Url, column `physPos` is optional except in some particular case (see below). If this column is provided (or contain only missing values), it should exactly match the physical positions of the SNP specified in the VCF file. If `SNPid` columns is missing or have missing values, the SNPid will be automatically imputed using the convention `chr@physPos` therefore columns `chr` and `physPos` should not have any missing values in this case.


## Value

path of the `.vcf.gz` file containing the simulated genotypes
 of the offspring.


# `calc_progenyBlupEstimation`

Calculate progenies BLUPs variance and expected values based on parents'
 genotype and markers effects


## Description

Calculate progenies BLUPs variance and expected values based on parents'
 genotype and markers effects


## Usage

```r
calc_progenyBlupEstimation(
  genoFile = NULL,
  genoUrl = NULL,
  crossTableFile = NULL,
  crossTableUrl = NULL,
  SNPcoordFile = NULL,
  SNPcoordUrl = NULL,
  markerEffectsFile = NULL,
  markerEffectsUrl = NULL,
  outFile = tempfile(fileext = ".json")
)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     phased VCF file path (ext `.vcf` or `.vcf.gz`)
`genoUrl`     |     url of a phased VCF file path (ext `.vcf` or `.vcf.gz`)
`crossTableFile`     |     path of the crossing table data file (`csv` file of 2 or 3 columns). It must contain the names of the variables as its first line. The column 1 and 2 will be interpreted as the parents ids. The optional third column will be interpreted as the offspring base name.
`crossTableUrl`     |     URL of a crossing table file
`SNPcoordFile`     |     path of the SNPs coordinates file (`csv` file). This `.csv` file should have 4 named columns: - `chr`: Chromosome holding the SNP - `physPos`: SNP physical position on the chromosome - `linkMapPos`: SNP linkage map position on the chromosome in Morgan - `SNPid`: SNP's IDs
`SNPcoordUrl`     |     URL of a SNP coordinate file
`markerEffectsFile`     |     path of the marker effects file (`csv` or `json` file).
`markerEffectsUrl`     |     URL of a marker effect file
`outFile`     |     `.json` file path where to save the data. If the file already exists, it will be overwritten.  For SNPcoordFile/Url, column `physPos` is optional except in some particular case (see below). If this column is provided (or contain only missing values), it should exactly match the physical positions of the SNP specified in the VCF file. If `SNPid` columns is missing or have missing values, the SNPid will be automatically imputed using the convention `chr@physPos` therefore columns `chr` and `physPos` should not have any missing values in this case.


## Value

data.frame containing the calculations results


# `draw_progBlupsPlot`

Draw a plot of the progenies BLUPs' expected values with error bars


## Description

X axis is the crosses, and Y axis the blups. The points are located at the
 expected value and the error bar length is the standard deviation.


## Usage

```r
draw_progBlupsPlot(
  progEstimFile = NULL,
  progEstimUrl = NULL,
  errorBarInterval = 0.95,
  y_axisName = "Genetic values",
  sorting = "alpha",
  trait = NULL,
  outFile = tempfile(fileext = ".html")
)
```


## Arguments

Argument      |Description
------------- |----------------
`progEstimFile`     |     path of the progeny BLUP estimation file generated by r-geno-tools-engine containing the blup estimations of the progenies of some crosses (`json` file).
`progEstimUrl`     |     URL of a progeny BLUP estimation file
`errorBarInterval`     |     length of XX% interval of interest represented by the error bars (default=0.95)
`y_axisName`     |     The Y axis name (default = "genetic values")
`sorting`     |     method to sort the individuals (X axis) can be: - "asc": sort the BLUP expected value in ascending order (from left to right) - "dec": sort the BLUP expected value in decreasing order (from left to right) - any other value will sort the individuals in alphabetical order (from left to right)
`outFile`     |     outFile path of the file containing the plot. If `NULL`, the output will not be written in any file. By default write in an temporary file.


## Value

plotly graph


# `draw_progBlupsPlot_2traits`

Draw a plotly graph of blups data for 2 traits


## Description

The points are located at the expected value and the ellipses
 size represent the `confidenceLevel` prediction interval.


## Usage

```r
draw_progBlupsPlot_2traits(
  progEstimFile = NULL,
  progEstimUrl = NULL,
  x_trait,
  y_trait,
  confidenceLevel = 0.95,
  x_suffix = "",
  y_suffix = "",
  ellipses_npoints = 100,
  outFile = tempfile(fileext = ".html")
)
```


## Arguments

Argument      |Description
------------- |----------------
`progEstimFile`     |     path of the progeny BLUP estimation file generated by r-geno-tools-engine containing the blup estimations of the progenies of some crosses (`json` file).
`progEstimUrl`     |     URL of a progeny BLUP estimation file
`x_trait`     |     name of the trait to show on the x axis
`y_trait`     |     name of the trait to show on the y axis
`confidenceLevel`     |     level of the prediction ellipses (default 0.95, ie 95 %  ellypses)
`x_suffix`     |     suffix to add to the x axis's name
`y_suffix`     |     suffix to add to the y axis's name
`ellipses_npoints`     |     number of points used to draw the ellipses (default 100)
`outFile`     |     outFile path of the file containing the plot. If `NULL`, the output will not be written in any file. By default write in an temporary file.


## Value

plotly graph


# `train_gs_model_main`

Train GS model


## Description

Train GS model


## Usage

```r
train_gs_model_main(
  genoFile = NULL,
  phenoFile = NULL,
  genoUrl = NULL,
  phenoUrl = NULL,
  trait,
  with_dominance,
  thresh_maf,
  outFile = tempfile(fileext = ".json")
)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     path of the geno data file (`.vcf` or `.vcf.gz` file)
`phenoFile`     |     path of the phenotypic data file (`csv` file). Individuals' name should be the first column of the file and no duplication is allowed.
`genoUrl`     |     url of the geno data file (`.vcf` or `.vcf.gz` file)
`phenoUrl`     |     url of the phenotypic data file (`csv` file) Individuals' name should be the first column of the file and no duplication is allowed.
`trait`     |     Chraracter of length 1, name of the trait to analyze. Must be a column name of the phenotypic file.
`with_dominance`     |     should the model include dominance effects
`thresh_maf`     |     threshold to keep only markers with minor allele frequency greater than `thresh_maf`.
`outFile`     |     paht of the `.json` file where to save the model's estimated markers effects.


# `predict_gs_model_main`

Make phenotypic prediction using a markers effects


## Description

Make phenotypic prediction using a markers effects


## Usage

```r
predict_gs_model_main(
  genoFile = NULL,
  genoUrl = NULL,
  markerEffectsFile = NULL,
  markerEffectsUrl = NULL,
  outFile = tempfile(fileext = ".csv")
)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     path of the geno data file (`.vcf` or `.vcf.gz` file)
`genoUrl`     |     url of the geno data file (`.vcf` or `.vcf.gz` file)
`markerEffectsFile`     |     path of the marker effects file (`csv` or `json` file).
`markerEffectsUrl`     |     URL of a marker effect file
`outFile`     |     `.csv` file path where to save the predictions. If the file already exists, it will be overwritten.


# `cross_validation_evaluation_main`

Evaluate a model with repeated cross validation


## Description

Evaluate a model with repeated cross validation


## Usage

```r
cross_validation_evaluation_main(
  genoFile = NULL,
  phenoFile = NULL,
  genoUrl = NULL,
  phenoUrl = NULL,
  trait,
  n_folds = 10,
  n_repetitions = 5,
  with_dominance,
  thresh_maf,
  outFile = tempfile(fileext = ".json")
)
```


## Arguments

Argument      |Description
------------- |----------------
`genoFile`     |     path of the geno data file (`.vcf` or `.vcf.gz` file)
`phenoFile`     |     path of the phenotypic data file (`csv` file). Individuals' name should be the first column of the file and no duplication is allowed.
`genoUrl`     |     url of the geno data file (`.vcf` or `.vcf.gz` file)
`phenoUrl`     |     url of the phenotypic data file (`csv` file) Individuals' name should be the first column of the file and no duplication is allowed.
`trait`     |     Chraracter of length 1, name of the trait to analyze. Must be a column name of the phenotypic file.
`n_folds`     |     number of fold for each cross-validation
`n_repetitions`     |     number of cross-validation repetitions
`with_dominance`     |     should the model include dominance effects
`thresh_maf`     |     threshold to keep only markers with minor allele frequency greater than `thresh_maf`.
`outFile`     |     paht of the `.json` file where to save the evaluation results


# `draw_evaluation_plot`

Draw the evaluation result plot


## Description

Draw the evaluation result plot


## Usage

```r
draw_evaluation_plot(
  evaluationFile = NULL,
  outFile = tempfile(fileext = ".html")
)
```


## Arguments

Argument      |Description
------------- |----------------
`evaluationFile`     |     path of the evaluation result file generated by `save_GS_evaluation` function
`outFile`     |     path of the file containing the plot. If `NULL`, the output will not be written in any file. By default write in an temporary file.


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
`filter_pAdj`     |     [numeric] threshold to remove points with pAdj < filter_pAdj from the plot (default no filtering)
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


# `relMatHeatmap`

Heatmap of a relationship matrix


## Description

Heatmap of a relationship matrix


## Usage

```r
relMatHeatmap(relMat, interactive = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
`relMat`     |     relationship matrix return by `pedRelMat`
`interactive`     |     [bool] should the plot be interactive (the default) or not


## Value

plotly graph if interactive is TRUE, or NULL if not.


# `pedNetwork`

Draw an interactive Pedigree network


## Description

Draw an interactive Pedigree network


## Usage

```r
pedNetwork(ped)
```


## Arguments

Argument      |Description
------------- |----------------
`ped`     |     List return by `readPedData` function


## Value

a `forceNetwork` object (`htmlwidget`)


# `plotBlup_1trait`

Draw a plotly graph of blups data for 1 trait


## Description

X axis is the crosses, and Y axis the blups. The points are located at the
 expected value and the error bar length is the standard deviation.


## Usage

```r
plotBlup_1trait(
  blupDta,
  sorting = "alpha",
  y_axisName = NULL,
  errorBarInterval = 0.95,
  trait
)
```


## Arguments

Argument      |Description
------------- |----------------
`blupDta`     |     list of data.frame of 4 columns: "ind1", "ind2", "blup_exp", "blup_var"
`sorting`     |     method to sort the individuals (X axis) can be: - "asc": sort the BLUP expected value in ascending order (from left to right) - "dec": sort the BLUP expected value in decreasing order (from left to right) - any other value will sort the individuals in alphabetical order (from left to right)
`y_axisName`     |     Name of the Y axis (default = `trait`)
`errorBarInterval`     |     length of XX% interval of interest represented by the error bars (default=0.95)
`trait`     |     name of the trait to plot. This should be a name of the blupDta list. (optional if only one trait in `blupDta`, it will be set to the name of this trait)


## Value

plotly graph


# `plotBlup_2traits`

Draw a plotly graph of blups data for 2 traits


## Description

The points are located at the expected value and the ellipses
 size represent the `confidenceLevel` prediction interval.


## Usage

```r
plotBlup_2traits(
  blupDta,
  x_trait,
  y_trait,
  confidenceLevel = 0.95,
  x_suffix = "",
  y_suffix = "",
  ellipses_npoints = 100
)
```


## Arguments

Argument      |Description
------------- |----------------
`blupDta`     |     list returned by calc_progenyBlupEstimation
`x_trait`     |     name of the trait to show on the x axis
`y_trait`     |     name of the trait to show on the y axis
`confidenceLevel`     |     level of the prediction ellipses (default 0.95, ie 95 %  ellypses)
`x_suffix`     |     suffix to add to the x axis's name
`y_suffix`     |     suffix to add to the y axis's name
`ellipses_npoints`     |     number of points used to draw the ellipses (default 100)


## Value

plotly graph


# `evaluation_plot`

Draw a plotly graph of a GS model cross-validation evaluation


## Description

This plots is composed of several subplots:
 - Observed vs Predicted scatter plot for all the cross-validation folds
 - Horizontal box plots of models metrics calculated during the
 cross-validation


## Usage

```r
evaluation_plot(evaluation_results)
```


## Arguments

Argument      |Description
------------- |----------------
`evaluation_results`     |     list returned by `cross_validation_evaluation()`


## Value

plotly graph


# `calcRecombRate`

Calculate the recombination rate matrix for each couple of SNP


## Description

The recombination rate is alculated using the "haldane inverse" function


## Usage

```r
calcRecombRate(SNPcoord)
```


## Arguments

Argument      |Description
------------- |----------------
`SNPcoord`     |     SNP coordinate data.frame return by `readSNPcoord`


## Value

named list of matrices. List names are the chromosomes' names.
 Matrices' row and columns names are the SNP ids
 (of the corresponding chromosome). Matrices values are the recombination rate
 between the corresponding SNPs.


# `calcProgenyGenetCovar`

Calculate the genetic variance-covariance matrix of the progenies of 2 given
 parents


## Description

Calculate the genetic variance-covariance matrix of the progenies of 2 given
 parents


## Usage

```r
calcProgenyGenetCovar(SNPcoord, r, haplo, p1.id, p2.id)
```


## Arguments

Argument      |Description
------------- |----------------
`SNPcoord`     |     SNP coordinate data.frame return by `readSNPcoord`
`r`     |     recombination rate matrices return by `calcRecombRate`
`haplo`     |     haplotypes of individuals ("haplotypes" element of the list return by `readPhasedGeno` function)
`p1.id`     |     id of the first parent
`p2.id`     |     id of the second parent


## Value

named list of matrices. List names are the chromosomes' names.
 Matrices' row and columns names are the SNP ids
 (of the corresponding chromosome). Matrices values are the genetic covariance
 between the corresponding SNPs for the progeny of the given parents.


# `calcProgenyBlupVariance`

Calculate the BULP variance of the progeny


## Description

Calculate the BULP variance of the progeny


## Usage

```r
calcProgenyBlupVariance(SNPcoord, markerEffects, geneticCovar)
```


## Arguments

Argument      |Description
------------- |----------------
`SNPcoord`     |     SNP coordinate data.frame return by `readSNPcoord`
`geneticCovar`     |     list of the genetic variance covariance matrices return by `calcProgenyGenetCovar`
`makrerEffects`     |     (output of `extract_additive_effects` function) list of 2 elements: `intercept`: named vector of the intercepts `SNPeffects`: data.frame of the additive effects, 1 columns per phenotype with the marker ids as row names.


## Value

numeric


# `calcProgenyBlupExpected`

Calculate the BULP expected values of the progeny


## Description

Calculate the BULP expected values of the progeny


## Usage

```r
calcProgenyBlupExpected(SNPcoord, haplo, p1.id, p2.id, markerEffects)
```


## Arguments

Argument      |Description
------------- |----------------
`SNPcoord`     |     SNP coordinate data.frame return by `readSNPcoord`
`haplo`     |     haplotypes of individuals ("haplotypes" element of the list return by `readPhasedGeno` function)
`p1.id`     |     id of the first parent
`p2.id`     |     id of the second parent
`makrerEffects`     |     (output of `extract_additive_effects` function) list of 2 elements: `intercept`: named vector of the intercepts `SNPeffects`: data.frame of the additive effects, 1 columns per phenotype with the marker ids as row names.


## Value

list for each trait with a list with
 - `sum` the global expected value for the trait (taking in account the intercept)
 - `by_chr` list of the expected value for each chromosome (NOT taking in account the intercept)


# `calcProgenyBlupCovariance`

Calculate the variance covariance matrix of the progeny's blups for several
 traits


## Description

Calculate the variance covariance matrix of the progeny's blups for several
 traits


## Usage

```r
calcProgenyBlupCovariance(
  SNPcoord,
  r,
  haplo,
  p1.id,
  p2.id,
  markerEffects,
  blupExpectedValues
)
```


## Arguments

Argument      |Description
------------- |----------------
`SNPcoord`     |     SNP coordinate data.frame return by `readSNPcoord`
`r`     |     recombination rate matrices return by `calcRecombRate`
`haplo`     |     haplotypes of individuals ("haplotypes" element of the list return by `readPhasedGeno` function)
`p1.id`     |     id of the first parent
`p2.id`     |     id of the second parent
`blupExpectedValues`     |     output of `calcProgenyBlupExpected`
`makrerEffects`     |     (output of `extract_additive_effects` function) list of 2 elements: `intercept`: named vector of the intercepts `SNPeffects`: data.frame of the additive effects, 1 columns per phenotype with the marker ids as row names.


## Value

matrix


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


# `downloadPhasedGeno`

Download phased geno data


## Description

Download phased geno data


## Usage

```r
downloadPhasedGeno(url)
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the geno data file (.vcf.gz file)


## Value

list of 2: `haplotypes` a matrix of the individuals haplotypes
 and `SNPcoord`, data frame of the SNP coordinates.


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


## Details

The individuals' names must be on the first column. No duplication
 is allowed.


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


# `downloadPedData`

Download pedigree data


## Description

Download pedigree data


## Usage

```r
downloadPedData(url, unknown_string = "", header = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the pedigree data file (`csv` file).
`unknown_string`     |     [default: ""] a character vector of strings which are to be interpreted as "unknown parent". By default: missing value in the file.
`header`     |     [default: TRUE] a logical value indicating whether the file contains the names of the variables as its first line. The default value is TRUE. In any cases, the column 1 will be interpreted as the individual id, column 2 as the first parent, column 3 as the second parent.


## Value

List of 2: `data` pedigree data, `graph` "igraph" object of the pedigree graph.


# `downloadRelMat`

Download relationship matrix


## Description

Download relationship matrix


## Usage

```r
downloadRelMat(url, format = tools::file_ext(url))
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the result data file (csv or json file)
`format`     |     format of the input file. Either "csv" or "json" (optional, by default it will use "json").


## Value

matrix


# `downloadCrossTable`

Download crossing table


## Description

Download crossing table


## Usage

```r
downloadCrossTable(url, header = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the crossing table data file (`csv` file of 2 or 3 columns).
`header`     |     [default: TRUE] a logical value indicating whether the file contains the names of the variables as its first line. The default value is TRUE. In any cases, the column 1 and 2 will be interpreted as the parents id. The optional third column will be interpreted as the offspring base name.


## Value

data.frame with the crossing table information.


# `downloadSNPcoord`

download SNP coordinates `.csv` file


## Description

download SNP coordinates `.csv` file


## Usage

```r
downloadSNPcoord(url)
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the SNPs coordinates file (`csv` file). This `.csv` file can have 4 named columns: - `chr`: Chromosome holding the SNP - `physPos`: SNP physical position on the chromosome - `linkMapPos`: SNP linkage map position on the chromosome in Morgan - `SNPid`: SNP's IDs  If `SNPid` columns is missing or have missing values, the SNPid will be automatically imputed using the convention `chr@physPos` therefore columns `chr` and `physPos` should not have any missing values


## Value

data.frame of 4 columns: 'chr', 'physPos', 'linkMapPos', 'SNPid'


# `downloadMarkerEffects`

Download marker effects file


## Description

Download marker effects file


## Usage

```r
downloadMarkerEffects(url)
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the marker effects file (`csv`, or `json` file).


## Value

data.frame of 1 columns named `effects` with the marker ids as
 row names.


# `downloadProgBlupEstim`

Download progeny BLUP estimation file


## Description

Download progeny BLUP estimation file


## Usage

```r
downloadProgBlupEstim(url)
```


## Arguments

Argument      |Description
------------- |----------------
`url`     |     url of the progeny BLUP estimation file generated by r-geno-tools-engine containing the blup estimations of the progenies of some crosses (`json` file).


## Value

data.frame of 4 columns named "ind1", "ind2", "blup_var", "blup_exp"


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


# `readPhasedGeno`

Read phased genetic data from a file


## Description

Read phased genetic data from a file


## Usage

```r
readPhasedGeno(file)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     phased VCF file path (ext `.vcf` or `.vcf.gz`)


## Value

list of 2: `haplotypes` a matrix of the individuals haplotypes
 and `SNPcoord`, data frame of the SNP coordinates.


# `readPhenoData`

Read phenotypic data file


## Description

Read phenotypic data file


## Usage

```r
readPhenoData(file, ind.names = 1, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     file path
`ind.names`     |     [default 1] a single number giving the column of the table which contains the individuals' names.
`...`     |     Further arguments to be passed to `read.csv`


## Details

Any duplication in the phenotypic file is forbidden.


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


# `readPedData`

Read and prepare pedigree data


## Description

Read and prepare pedigree data


## Usage

```r
readPedData(file, unknown_string = "", header = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the pedigree data file (`csv` file).
`unknown_string`     |     [default: ""] a character vector of strings which are to be interpreted as "unknown parent". By default: missing value in the file.
`header`     |     [default: TRUE] a logical value indicating whether the file contains the names of the variables as its first line. The default value is TRUE. In any cases, the column 1 will be interpreted as the individual id, column 2 as the first parent, column 3 as the second parent.


## Details

We consider here only allo-fecundation or auto-fecundation. For
 auto-fecundation, use the parental individual id in both column 2 and 3.
 Doubles haploids can not be interpreted, please avoid them in the file.
 
 Please be sure that all individuals id in columns 2 and 3 are defined in the
 column 1. If columns 2 and/or 3 contain id of individuals
 that are not in the first column, a warning will be raised and these
 individuals will be added to the pedigree with unknown
 parents as founder individuals.


## Value

List of 2: `data` pedigree data, `graph` "igraph" object of the pedigree graph.


# `readRelMat`

Read a relationship matrix file


## Description

Read a relationship matrix file


## Usage

```r
readRelMat(file, format = tools::file_ext(file))
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the file generated by the function `saveRelMat` containing relationship matrix
`format`     |     format of the input file. Either "csv" or "json" (optional, by default it will use the `file` extension).


## Details

The metadata of the file are not kept.


## Value

matrix


# `readCrossTable`

Read crossing table


## Description

Read crossing table


## Usage

```r
readCrossTable(file, header = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the crossing table data file (`csv` file of 2 or 3 columns).
`header`     |     [default: TRUE] a logical value indicating whether the file contains the names of the variables as its first line. The default value is TRUE. In any cases, the column 1 and 2 will be interpreted as the parents id. The optional third column will be interpreted as the offspring base name.


## Value

data.frame with the crossing table information.


# `readSNPcoord`

Read SNP coordinates `.csv` file


## Description

Read SNP coordinates `.csv` file


## Usage

```r
readSNPcoord(file)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the SNPs coordinates file (`csv` file). This `.csv` file can have 4 named columns: - `chr`: Chromosome holding the SNP - `physPos`: SNP physical position on the chromosome - `linkMapPos`: SNP linkage map position on the chromosome in Morgan (mandatory) - `SNPid`: SNP's IDs  Column `physPos` is optional except in some particular case (see below)  If `SNPid` columns is missing or have missing values, the SNPid will be automatically imputed using the convention `chr@physPos` therefore columns `chr` and `physPos` should not have any missing values


## Value

data.frame of 4 columns: 'chr', 'linkMapPos', 'SNPid'


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


# `readMarkerEffects`

Read marker effects file


## Description

Read marker effects file


## Usage

```r
readMarkerEffects(file)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the marker effects file (`csv`, or `json` file)


## Details

For `.csv`, the file file should
 have 2 named columns:
 - `SNPid`: Marker id
 - `effects`: effect of the corresponding marker
 one can specify the intercept using a "--INTERCEPT--" as SNPid.
 
 For `.json`, the file should
 have 2 Key-value pairs:
 - `intercept`: a number with the value of the intercept.
 - `coefficient`: a nested object with SNPids as keys and their corresponding
 effects as values.
 For example :
 list("\n", "  \"intercept\": 100,\n", "  \"coefficients\": ", list("\n", "    \"SNP01\": 1.02e-06,\n", "    \"SNP02\": 0.42,\n", "    \"SNP03\": 0.0\n", "  "), "\n")


## Value

list of 2 elements:
 `intercept`: the value of the intercept,
 `effects`: data.frame of 1 columns named `SNPeffects` with the marker ids as
 row names.


# `readMarkerEffects_csv`

Read marker effects CSV file


## Description

Read marker effects CSV file


## Usage

```r
readMarkerEffects_csv(file)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the marker effects file (`csv` file). This `.csv` file should have 2 named columns: - `SNPid`: Marker id - `effects`: effect of the corresponding marker one can specify the intercept using a "--INTERCEPT--" as SNPid.


## Value

list of 2 elements:
 `intercept`: the value of the intercept,
 `effects`: data.frame of 1 columns named `SNPeffects` with the marker ids as
 row names.


# `readMarkerEffects_json`

Read marker effects JSON file


## Description

Read marker effects JSON file


## Usage

```r
readMarkerEffects_json(file)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the marker effects file (`json` file). This `.json` file should have 2 Key-value pairs: - `intercept`: a number with the value of the intercept. - `coefficient`: a nested object with SNPids as keys and their corresponding effects as values. For example : list("\n", "  \"intercept\": 100,\n", "  \"coefficients\": ", list("\n", "    \"SNP01\": 1.02e-06,\n", "    \"SNP02\": 0.42,\n", "    \"SNP03\": 0.0\n", "  "), "\n")


## Value

list of 2 elements:
 `intercept`: the value of the intercept,
 `effects`: data.frame of 1 columns named `SNPeffects` with the marker ids as
 row names.


# `readProgBlupEstim`

Read progeny BLUP estimation file


## Description

Read progeny BLUP estimation file


## Usage

```r
readProgBlupEstim(file)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the progeny BLUP estimation file generated by r-geno-tools-engine containing the blup estimations of the progenies of some crosses (`json` file).


## Value

data.frame of 4 columns named "ind1", "ind2", "blup_var", "blup_exp"


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

path of the created file


# `saveRelMat`

Save relationship matrix in file


## Description

Save relationship matrix in file


## Usage

```r
saveRelMat(
  relMat,
  metadata = NULL,
  dir = NULL,
  file = NULL,
  format = tools::file_ext(file)
)
```


## Arguments

Argument      |Description
------------- |----------------
`relMat`     |     relationship matrix created with `pedRelMat`
`metadata`     |     list of metadata of the relationship matrix (optional).
`dir`     |     if `file` is NULL, directory where to save the data, by default it is a temporary directory
`file`     |     file path where to save the data. If the file already exists, it will be overwritten. Default NULL (it will create a new "csv" file)
`format`     |     format of the output file. Either "csv" or "json". (optional, by default it will use the `file` extension, if file is NULL, "csv").


## Value

path of the created file


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


# `saveVcf`

Save phased genotypes of simulatied population to vcf.gz file


## Description

Save phased genotypes of simulatied population to vcf.gz file


## Usage

```r
saveVcf(file, pop, SNPcoord)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     file path where to save the data. If the file already exists, it will be overwritten.
`pop`     |     simulated population (`breedSimulatR`'s population)
`SNPcoord`     |     snp coordinate of the genotypes. (data.frame with `chr`, `physPos`, and `SNPid` columns)


# `save_dataFrame_as_json`

Save a R data.frame as json file


## Description

Save a R data.frame as json file


## Usage

```r
save_dataFrame_as_json(df, file)
```


## Arguments

Argument      |Description
------------- |----------------
`df`     |     data.frame
`file`     |     file path where to save the data. If the file already exists, it will be overwritten.


## Value

path of the created file


# `save_blupVarExp_as_json`

Save a R blupVarExp as json file


## Description

Save a R blupVarExp as json file


## Usage

```r
save_blupVarExp_as_json(blupVarExp, file)
```


## Arguments

Argument      |Description
------------- |----------------
`blupVarExp`     |     list (cf. `calc_progenyBlupEstimation()`)
`file`     |     file path where to save the data. If the file already exists, it will be overwritten.


## Value

path of the created file


# `save_plotly`

Save a plotly graph


## Description

Save a plotly graph


## Usage

```r
save_plotly(plot, file)
```


## Arguments

Argument      |Description
------------- |----------------
`plot`     |     graph generated with plotly
`file`     |     file path where to save the html graph. If the file already exists, it will be overwritten.


## Value

path of the created file


# `save_GS_model`

Save a GS model


## Description

Save a GS model


## Usage

```r
save_GS_model(model, trait_name, file)
```


## Arguments

Argument      |Description
------------- |----------------
`model`     |     GS model generated with
`trait_name`     |     name of the phenotypic trait predicted by this model
`file`     |     file path where to save the model. If the file already exists, it will be overwritten.


## Value

path of the created file


# `extract_additive_effects`

Format `markerEffects` to be used by `calcProgenyBlupCovariance` function.


## Description

Format `markerEffects` to be used by `calcProgenyBlupCovariance` function.


## Usage

```r
extract_additive_effects(markerEffects)
```


## Arguments

Argument      |Description
------------- |----------------
`markerEffects`     |     list


## Value

list of 2 elements:
 `intercept`: named vector of the intercepts
 `SNPeffects`: data.frame of the additive effects, 1 columns per phenotype with
 the marker ids as row names.


# `save_GS_evaluation`

Save evaluation result in .json file


## Description

Save evaluation result in .json file


## Usage

```r
save_GS_evaluation(evaluation, trait_name, file)
```


## Arguments

Argument      |Description
------------- |----------------
`evaluation`     |     output of `cross_validation_evaluation` function
`trait_name`     |     trait name
`file`     |     file path where to save the evaluations resutls. If the file already exists, it will be overwritten.


# `read_GS_evaluation`

read evaluation results file


## Description

read evaluation results file


## Usage

```r
read_GS_evaluation(file)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     path of the evaluation result file generated by `save_GS_evaluation` function


# `pedRelMat`

Pedigree Relationship Matrix calculation


## Description

Pedigree Relationship Matrix calculation


## Usage

```r
pedRelMat(ped)
```


## Arguments

Argument      |Description
------------- |----------------
`ped`     |     List return by `readPedData` function


## Value

matrix


## Author

Hiroyoshi Iwata, Julien Diot


# `calc_additive_geno`

Additive Genomic Matrix


## Description

Additive Genomic Matrix


## Usage

```r
calc_additive_geno(geno, standardized = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
`geno`     |     `gaston::bed.matrix` return by `readGenoData` function
`standardized`     |     `boolean` (default TRUE) control if the returned genetic matrix should be standardized.


## Details

The standardization is made with: (x - 2*p) / sqrt(2*p*(1 - p)) with x the
 genetic values in alleles dose and p the allelic frequency.


## Value

matrix


# `calc_additive_rel_mat`

Additive Genomic Relationship Matrix calculation


## Description

Additive Genomic Relationship Matrix calculation


## Usage

```r
calc_additive_rel_mat(geno, standardized = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
`geno`     |     `gaston::bed.matrix` return by `readGenoData` function
`standardized`     |     `boolean` (default TRUE) control if the calculation should be done on the standardized additive genetic matrix.


## Value

list of 2 elements:
 - `geno_mat`: the genetic matrix used for calculating the relationship matrix
 - `rel_mat`: the relationship matrix


# `calc_dominance_geno`

Dominance Genomic Matrix


## Description

Dominance Genomic Matrix


## Usage

```r
calc_dominance_geno(geno, standardized = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
`geno`     |     `gaston::bed.matrix` return by `readGenoData` function
`standardized`     |     `boolean` (default TRUE) control if the returned genetic matrix should be standardized.


## Details

The standardization is made with the matrix with entries
 `p/(1 - p)`, `-1`, `(1-p)/p` according to the values `0`, `1`, `2` in the
 genetic matrix, with p the allele frequency (cf. `?gaston::DM`)


## Value

matrix


# `calc_dominance_rel_mat`

Dominance Genomic Relationship Matrix calculation


## Description

Dominance Genomic Relationship Matrix calculation


## Usage

```r
calc_dominance_rel_mat(geno, standardized = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
`geno`     |     `gaston::bed.matrix` return by `readGenoData` function
`standardized`     |     `boolean` (default TRUE) control if the calculation should be done on the standardized dominance genetic matrix.


## Value

list of 2 elements:
 - `geno_mat`: the genetic matrix used for calculating the relationship matrix
 - `rel_mat`: the relationship matrix


# `combinedRelMat`

Combined (pedigree + genomic) Relationship Matrix


## Description

Correct a pedigree relationship matrix by using genomic relationship matrix.


## Usage

```r
combinedRelMat(ped_rm, geno_rm, method = "Legarra", tau = NULL, omega = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`ped_rm`     |     pedigree relationship matrix (matrix from function pedRelMat)
`geno_rm`     |     genomic relationship matrix (matrix from function genoRelMat)
`method`     |     method to use, either "Legarra" or "Martini"
`tau`     |     tau parameter of the Martini's method
`omega`     |     omega parameter of the Martini's method


## Value

matrix


## Author

Hiroyoshi Iwata, Julien Diot


# `supThisWarning`

Do not show a specific warning message


## Description

Do not show a specific warning message


## Usage

```r
supThisWarning(expr, warnMessage)
```


## Arguments

Argument      |Description
------------- |----------------
`expr`     |     expression to evaluate.
`warnMessage`     |     warning message to catch


## Details

This is based on the `base::suppressWarnings` function.


