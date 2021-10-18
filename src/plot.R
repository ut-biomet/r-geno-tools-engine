#' create manhatan plot
#' @param gwas [data.frame] output of the gwas function
#' @param adj_method correction method: "holm", "hochberg",
#' "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
#' @param thresh_p p value significant threshold (default 0.05)
#' @param chr [char] name of the chromosome to show (show all if NA)
#' @param title [char] Title of the plot. Default is "Manhattan Plot"
#' @param interactive [bool] should the plot be interactive (the default)
#'   or not.
#'
#' @return plotly graph if interactive is TRUE, or a
manPlot <- function(gwas,
                    adj_method,
                    thresh_p = 0.05,
                    chr = NA,
                    title = "Manhattan Plot",
                    interactive = TRUE) {

  logger <- logger$new("r-manPlot()")


  # Check chromosome name
  logger$log("Check chromosome name ...")
  if (!is.na(chr)) {
    if (!chr %in% unique(gwas$chr)) {
      stop("`chr` should match the chromosome names of `gwas`")
    }
  }
  logger$log("Check chromosome name DONE")

  # P-Values adjustment
  logger$log("Adjust p-values ...")
  adj <- adjustPval(gwas$p, adj_method, thresh_p)
  gwas$p_adj <- adj$p_adj
  thresh_pAdj <- adj$thresh_adj
  logger$log("Adjust p-values DONE")

  # filter according to "chr"
  if (!is.na(chr)) {
    gwas <- gwas[as.character(gwas$chr) %in% chr,]
  }

  # plot functions require chr as numeric values:
  chrlabels <- unique(gwas$chr)
  gwas$chr <- as.numeric(factor(gwas$chr,
                                levels = chrlabels))

  # manage duplicate in SNP's ID
  logger$log("Check duplicated SNP ID ...")
  if (anyDuplicated(gwas$id) != 0) {
    warning("Duplicated in SNP's ID detected, replacing SNP ID by: CHR@POS.")
    gwas$id <- paste0(gwas$chr, "@", gwas$pos)
    if (anyDuplicated(gwas$id) != 0) {
      warning("After eplacing SNP ID by: CHR@POS,",
              "there is still some duplicates. No SNP will be highlighted.")
      highlightSinif <- FALSE
    } else  highlightSinif <- TRUE
  } else  highlightSinif <- TRUE
  logger$log("Check duplicated SNP ID DONE")

  # get significant SNP
  if (highlightSinif) {
    logger$log("Extract significant SNP ...")
    significantSNP <- gwas[gwas$p_adj <= thresh_p, "id"]
    if (length(significantSNP) == 0) {
      significantSNP <- NULL
    }
    logger$log("Extract significant SNP DONE")
  } else significantSNP <- NULL

  logger$log("Draw plot ...")
  if (interactive) {
  p <- manhattanly::manhattanly(
    data.frame(CHR = gwas$chr,
               BP = gwas$pos,
               SNP = gwas$id,
               P = gwas$p),
    snp = "SNP",
    labelChr = chrlabels,
    highlight = significantSNP,
    genomewideline = -log10(thresh_pAdj),
    suggestiveline = FALSE,
    title = title
  )
  logger$log("Draw plot DONE")
  logger$log("DONE, return output")
  } else {
    if (length(chrlabels) == 1) {
      # qqman::manhattan raise error if there is only one chromosome in the data
      # and chrlabs is defined
      chrlabels <- NULL
    }
    qqman::manhattan(x = gwas,
                     chr = 'chr',
                     bp = 'pos',
                     p = 'p',
                     snp = 'id',
                     chrlabs = chrlabels,
                     suggestiveline = FALSE,
                     genomewideline = -log10(thresh_pAdj),
                     annotatePval = -log10(thresh_pAdj),
                     main = title)
    p <- NULL
  }
  p
}





#'writeLDplot Compute r2 Linkage Disequilibrium (LD) between given SNPs and return a plot
#'
#' @param geno [bed.matrix] geno data return by function `readGenoData` or
#' `downloadGenoData`.
#' @param from lower bound of the range of SNPs for which the LD is computed
#' @param to upper bound of the range of SNPs for which the LD is computed
#' @param file path of the png file to save the plot. If `NULL`, the image file will not be
#' created. By default write in an tempoary `.png` file.
#'
#' @details `from` should be lower than `to`, and the maximum ranger size is 50.
#' (In order to get a readable image). If write is `TRUE`, the function will write the plot in a png file, else it will plot it.
#' @return null if `dir` is NULL, else the path of the png file.
LDplot <- function(geno, from, to, file = tempfile(fileext = ".png")) {
  logger <- logger$new("r-LDplot()")

  # Checks:
  logger$log('Check "from" and "to" format ...')
  if (length(from) != 1) {
    stop(
         "`from` should be an positive integer of length 1 ",
         "It's length is \"",
         length(from),
         "\""
    )
  }
  if (!is.numeric(from) || from < 0 || from%%1!=0) {
    stop("`from` should be an positive integer of length 1 ",
         "It is \"",
         class(from),
         "\"")
  }
  if (from < 0 || from%%1!=0) {
    stop("`from` should be an positive integer of length 1 ",
         "It's value is \"",
         from,
         "\"")
  }

  if (length(to) != 1) {
    stop(
         "`to` should be an positive integer of length 1 ",
         "It's length is \"",
         length(to),
         "\""
    )
  }
  if (!is.numeric(to) || to < 0 || to%%1!=0) {
    stop("`to` should be an positive integer of length 1 ",
         "It is \"",
         class(to),
         "\"")
  }
  if (to < 0 || to%%1!=0) {
    stop("`to` should be an positive integer of length 1 ",
         "It's value is \"",
         to,
         "\"")
  }
  logger$log('Check "from" and "to" format DONE')

  logger$log('Check "from" < "to"...')
  if (from >= to) {
    logger$log('Error: "from" greater than "to".')
    stop('"from" should be lower than "to"')
  }
  logger$log('Check "from" < "to" DONE')


  logger$log('Check number of SNP < 50...')
  if (to - from > 50) {
    logger$log('Error: number of SNP is > 50.')
    stop('range size ("to" - "from") should be lower or equal than 50')
  }
  logger$log('Check number of SNP < 50 DONE')

  if (!is.null(file)) {
    logger$log('Check file ...')
    if (length(file) != 1) {
      logger$log('Error: only one file name should be provided')
      stop('Error: only one file name should be provided')
    }
    if (file.exists(file)) {
      logger$log('Warning: "file" directory already exists. This file will be overwritten.')
    } else {
      file.create(file)
    }
    logger$log('Check file DONE')
  }

  # COMPUTE LD:
  logger$log("Compute LD ...")
  ld <- gaston::LD(geno, c(from, to), measure = "r2")
  logger$log("Compute LD DONE")
  logger$log("Create LD plot ...")
  if (!is.null(file)) {
    png(file, width = 2400, height = 1600)
    logger$log("Create create file:", file)
  }
  gaston::LD.plot(ld, snp.positions = geno@snps$pos[from:to])
  if (!is.null(file)) {
    dev.off()
  }
  logger$log("Create LD plot DONE")

  logger$log("DONE, return output")
  if (!is.null(file)) {
    return(file)
  } else {return(NULL)}
}
