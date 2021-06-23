#' create manhatan plot
#' @param gwas [data.frame] output of the gwas function
#' @param adj_method correction method: "holm", "hochberg", "hommel",
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
  logger$log("Adjust p-values ...")
  gwas$p_adj <- p.adjust(gwas$p,
                        method = adj_method,
                        n = nrow(gwas))
  maxSinifP <- max(gwas[gwas$p_adj < thresh_p, "p"], thresh_p/nrow(gwas)/100)
  minUnSinifP <- min(gwas[gwas$p_adj > thresh_p, "p"], 1)
  logger$log("Adjust p-values DONE")

  # threshold adjustment
  logger$log("Adjust p-values threshold...")
  thresh_pAdj <- uniroot(
    function(log10p){
      p <- 10^(log10p)
      p.adjust(p,
               method = adj_method,
               n = nrow(gwas)) - thresh_p
    },
    c(log10(maxSinifP), log10(minUnSinifP))
  )
  thresh_pAdj <- 10^(thresh_pAdj$root)
  logger$log("Adjust p-values threshold DONE")

  # filter according to "chr"
  if (!identical(chr, NA)) {
    gwas <- gwas[as.character(gwas$chr) %in% chr,]
  }

  # get significant SNP
  significantSNP <- gwas[gwas$p_adj <= thresh_p, "id"]

  p <- manhattanly(
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
#' @param write should the function write the plot in a png file. Default: `false`
#'
#' @details `from` should be lower than `to`, and the maximum ranger size is 50.
#' (In order to get a readable image). If write is `true`, the function will write the plot in a temporary png file, else it will plot it.
#' @return null if `write` is false, else the path of the png file.
#' @export
#'
#' @examples
LDplot <- function(geno, from, to, write = FALSE) {
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

  # COMPUTE LD
  logger$log("Compute LD ...")
  ld <- LD(geno, c(from, to), measure = "r2")
  logger$log("Compute LD DONE")

  logger$log("Create LD plot ...")
  imgFile <- tempfile(fileext = ".png")
  logger$log("Create create file:", imgFile)
  png(imgFile, width = 2400, height = 1600)
  LD.plot(ld, snp.positions = geno@snps$pos[from:to])
  dev.off()
  logger$log("Create LD plot DONE")

  logger$log("DONE, return output")
  return(imgFile)
}
