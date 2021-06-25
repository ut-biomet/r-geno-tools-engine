#' create manhatan plot
#' @param gwas [data.frame] output of the gwas function
#' @param adj_method correction method: "holm", "hochberg",
#' "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
#' @param thresh_p significant threshold (default 0.05)
#' @param chr [char] filter to show only some chromosomes (show all if NA)
#' @param title [char] Title of the plot. Default is "Manhattan Plot"
#'
#' @return plotly graph
manPlot <- function(gwas,
                    adj_method,
                    thresh_p = 0.05,
                    chr = NA,
                    title = "Manhattan Plot") {

  logger <- logger$new("r-manPlot()")


  # P-Values adjustment
  adj <- adjustPval(gwas$p, adj_method, thresh_p)
  gwas$p_adj <- adj$p_adj
  thresh_pAdj <- adj$thresh_adj
  logger$log("Adjust p-values threshold DONE")

  # filter according to "chr"
  if (!is.na(chr)) {
    gwas <- gwas[as.character(gwas$chr) %in% chr,]
  }

  # get significant SNP
  significantSNP <- gwas[gwas$p_adj <= thresh_p, "id"]
  if (length(significantSNP) == 0) {
    significantSNP <- NULL
  }

  p <- manhattanly::manhattanly(
    data.frame(CHR = as.numeric(factor(gwas$chr,
                                       levels = unique(gwas$chr))),
               BP = gwas$pos,
               SNP = gwas$id,
               P = gwas$p),
    snp = "SNP",
    labelChr = unique(gwas$chr),
    highlight = significantSNP,
    genomewideline = -log10(thresh_pAdj),
    suggestiveline = FALSE,
    title = title
  )
  logger$log("DONE, return output")
  p
}





#'writeLDplot Compute r2 Linkage Disequilibrium (LD) between given SNPs and return a plot
#'
#' @param geno [bed.matrix] geno data return by function `readGenoData` or
#' `downloadGenoData`.
#' @param from lower bound of the range of SNPs for which the LD is computed
#' @param to upper bound of the range of SNPs for which the LD is computed
#' @param write should the function write the plot in a png file.
#'  Default: `FALSE`
#' @param dir if write is TRUE, path of the directory where to save the file, by default dir is set to a temporary directory.
#'  the image will be written in a temporary dir. Default: `NULL`
#'
#' @details `from` should be lower than `to`, and the maximum ranger size is 50.
#' (In order to get a readable image). If write is `TRUE`, the function will write the plot in a png file, else it will plot it.
#' @return null if `write` is false, else the path of the png file.
#' @export
#'
#' @examples
LDplot <- function(geno, from, to, write = FALSE, dir = tempdir()) {
  logger <- logger$new("r-LDplot()")

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

  logger$log('Check write ...')
  if (!is.logical(write)) {
    logger$log('Error: "write" should be a boolean')
    stop('Error: "write" should be a boolean')
  }
  logger$log('Check write DONE')

  if (write) {
    logger$log('Check dir ...')
    if (!dir.exists(dir)) {
      logger$log('Error: "dir" directory should exists')
      stop('Error: "dir" directory should exists')
    }
    logger$log('Check dir DONE')
  }



  # COMPUTE LD
  logger$log("Compute LD ...")
  ld <- gaston::LD(geno, c(from, to), measure = "r2")
  logger$log("Compute LD DONE")

  logger$log("Create LD plot ...")
  if (write) {
    imgFile <- tempfile(fileext = ".png", tmpdir = dir)
    png(imgFile, width = 2400, height = 1600)
    logger$log("Create create file:", imgFile)
  }
  gaston::LD.plot(ld, snp.positions = geno@snps$pos[from:to])
  if (write) {
    dev.off()
  }
  logger$log("Create LD plot DONE")

  logger$log("DONE, return output")
  if (write) {
  return(imgFile)
  } else {return(NULL)}
}
