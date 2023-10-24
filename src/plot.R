#' create manhatan plot
#' @param gwas [data.frame] output of the gwas function
#' @param adj_method correction method: "holm", "hochberg",
#' "bonferroni", "BH", "BY", "fdr", "none" (see ?p.adjust for more details)
#' @param thresh_p p value significant threshold (default 0.05)
#' @param chr [char] name of the chromosome to show (show all if NA)
#' @param title [char] Title of the plot. Default is "Manhattan Plot"
#' @param interactive [bool] should the plot be interactive (the default) or not
#' @param filter_pAdj [numeric] threshold to remove points
#' with pAdj < filter_pAdj from the plot (default no filtering)
#' @param filter_nPoints [numeric] threshold to keep only the filter_nPoints
#' with the lowest p-values for the plot (default no filtering)
#' @param filter_quant [numeric] threshold to keep only the filter_quant*100 %
#' of the points with the lowest p-values for the plot (default no filtering)
#'
#' @details If several filtering rules are given, the filtering process apply
#' the filtering process sequentially (this lead to having the same result
#' that if only the strongest rules were given).
#' Moreover, the number of points kept for the plot will be display
#' in the plot title.
#' @return plotly graph if interactive is TRUE, or NULL if not.
manPlot <- function(gwas,
                    adj_method,
                    thresh_p = 0.05,
                    chr = NA,
                    title = "Manhattan Plot",
                    filter_pAdj = 1,
                    filter_nPoints = Inf,
                    filter_quant = 1,
                    interactive = TRUE) {

  logger <- Logger$new("r-manPlot()")


  # Check parameters ----
  logger$log("Check parameters...")
  if (!is.na(as.numeric(thresh_p))) {
    thresh_p <- as.numeric(thresh_p)
  } else {
    stop('`thresh_p` should be a numeric value')
  }
  logger$log("Check parameters DONE")


  # Check chromosome name ----
  logger$log("Check chromosome name ...")
  if (!is.na(chr)) {
    if (!chr %in% unique(gwas$chr)) {
      stop("`chr` should match the chromosome names of `gwas`")
    }
  }
  logger$log("Check chromosome name DONE")

  logger$log("Remove NAs ...")
  # remove NA lines form data (important for cases where all the data had
  # been removed at the gaws step)
  allNaLines <- apply(gwas, 1, function(x){all(is.na(x))})
  gwas <- gwas[!allNaLines,]
  logger$log("Remove NAs DONE")

  # P-Values adjustment ----
  logger$log("Adjust p-values ...")
  adj <- adjustPval(gwas$p, adj_method, thresh_p)
  gwas$p_adj <- adj$p_adj
  thresh_pAdj <- adj$thresh_adj
  logger$log("Adjust p-values DONE")


  # filter according to "chr" ----
  if (!is.na(chr)) {
    gwas <- gwas[as.character(gwas$chr) %in% chr,]
  }

  # plot functions require chr as numeric values:
  chrlabels <- unique(gwas$chr)
  chrlabels <- stringr::str_sort(chrlabels, numeric = TRUE)
  gwas$chr <- as.numeric(factor(gwas$chr,
                                levels = chrlabels))

  # manage duplicate in SNP's ID ----
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

  # get significant SNP ----
  if (highlightSinif) {
    logger$log("Extract significant SNP ...")
    significantSNP <- gwas[gwas$p_adj <= thresh_p, "id"]
    if (length(significantSNP) == 0) {
      significantSNP <- NULL
    }
    logger$log("Extract significant SNP DONE")
  } else significantSNP <- NULL

  # filter results ----
  nTotalSnp <- nrow(gwas)
  gwas <- filterGWAS(gwas = gwas,
                     filter_pAdj = filter_pAdj,
                     filter_nPoints = filter_nPoints,
                     filter_quant = filter_quant)

  # update plot title to give filtering feed back to the user
  remainPoints <- nrow(gwas)
  if (nTotalSnp != 0) {
    remainPercent <- round(remainPoints/nTotalSnp, digits = 3) * 100
    if (remainPoints != nTotalSnp && remainPercent != 0) {
      title <- paste(title,
                     paste0(remainPoints, ' points ~', remainPercent, '%)'),
                     sep = '\n(')
    } else if (remainPercent == 0) {
      title <- paste(title,
                     'filtering process removed all the points',
                     sep = '\n')
    }
  } else {
    warning('There is no points to display')
    title <- paste(title,
                   'no point to display',
                   sep = '\n')
  }

  # Draw plot ----
  logger$log("Draw plot ...")
  if (interactive) {
    if (remainPoints == 0) {
      # manhattanly::manhattanly can't handle empty data
      p <- plotly::plot_ly(type = "scatter",
                           mode = "lines+markers") %>%
        plotly::layout(title = title,
                       xaxis = list(title = 'Chromosome',
                                    zerolinecolor = '#ffff',
                                    gridcolor = 'ffff',
                                    showticklabels=FALSE),
                       yaxis = list(title = '-log10(p)',
                                    # zerolinecolor = '#ffff',
                                    # gridcolor = 'ffff',
                                    showticklabels=FALSE))
    } else {
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
        title = title)

      if (!is.na(chr)) {
        p <- plotly::layout(
          p,
          xaxis = list(title = paste('Chromosome', chr, 'position'))
        )
      }
    }
  } else {
    if (remainPoints == 0) {
      # qqman::manhattan can't handle empty data
      plot(1, type = "n",
           xlab = "Chromosome",
           ylab = expression(-log[10](italic(p))),
           xlim = c(0, length(chrlabels)),
           xaxt = 'n',
           yaxt = 'n',
           main = title)
      p <- NULL
    } else {
      if (length(chrlabels) == 1) {
        # qqman::manhattan raise error if there is only one chromosome
        # in the data and chrlabs is defined
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
  }
  logger$log("Draw plot DONE")
  logger$log("DONE, return output")
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
  logger <- Logger$new("r-LDplot()")

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





#' Heatmap of a relationship matrix
#'
#' @param relMat relationship matrix return by `pedRelMat`
#' @param interactive [bool] should the plot be interactive (the default) or not
#'
#' @return plotly graph if interactive is TRUE, or NULL if not.
relMatHeatmap <- function(relMat, interactive = TRUE) {
  logger <- Logger$new("r-manPlot()")

  # Check parameters ----
  logger$log('Check parameters ...')
  if (!is.matrix(relMat)) {
    stop('`relMat` should be a matrix.')
  }
  if (!isSymmetric(relMat)) {
    stop('`relMat` should be a symetric matrix.')
  }
  logger$log('Check parameters DONE')

  if (interactive) {
    # Create inreractive heatmap
    logger$log('Create interactive heatmap ...')
    newOrder <- hclust(dist(relMat))$order
    hm <- plotly::plot_ly(type = "heatmap",
                  x = colnames(relMat[newOrder, newOrder]),
                  y = rownames(relMat[newOrder, newOrder]),
                  z = relMat[newOrder, newOrder],
                  colors = hcl.colors(12, "YlOrRd", rev = TRUE)
    )
    hm <- plotly::layout(hm,
                         xaxis = list(
                           type = "category"
                         ),
                         yaxis = list(
                           type = "category",
                           autorange = "reversed")
    )
    logger$log('Create interactive heatmap DONE')
  } else {
    logger$log('Create static heatmap ...')
    heatmap(relMat,
            symm = TRUE)
    hm <- NULL
    logger$log('Create static heatmap DONE')
  }

  logger$log("DONE, return output")
  return(hm)
}





#' Draw an interactive Pedigree network
#'
#' @param ped List return by `readPedData` function
#'
#'
#' @return a `forceNetwork` object (`htmlwidget`)
pedNetwork <- function(ped) {
  logger <- Logger$new("r-pedNetwork()")

  # Check parameters ----
  logger$log('Check parameters ...')

  ### Check input ----
  logger$log('Check inputs ...')
  if (!is.list(ped)) {
    stop(
      "`ped` should be a list generated by `readPedData` function. ",
      "`ped` is \"",
      class(ped),
      "\""
    )
  }
  if (!identical(names(ped), c("data", "graph"))) {
    stop(
      "`ped` should be a list generated by `readPedData` function. ",
      "Therfore names of `ped` should be \"data\", \"graph\".",
      "It is: \"", paste(names(ped), collapse = '", "'),"\""
    )
  }
  if (!is.data.frame(ped$data)) {
    stop(
      "`ped` should be a list generated by `readPedData` function. ",
      "Therfore `ped$data` should be a \"data.frame\". It is: \"",
      class(ped$data),
      "\""
    )
  }
  if (!identical(colnames(ped$data), c('ind', 'parent1', 'parent2'))) {
    stop(
      "`ped` should be a list generated by `readPedData` function. ",
      'Therfore `colnames(ped$data)` should be:',
      ' "ind", "parent1", "parent2". It is: "',
      paste(colnames(ped$data), collapse = '", "'), '".'
    )
  }
  logger$log('Check inputs DONE')


  logger$log('Create network data ...')
  data <- data.frame(parent = c(ped$data[, 2], ped$data[, 3]),
                     ind = rep(ped$data[, 1], 2))
  data$ind <- as.factor(data$ind)
  data$parent <- factor(data$parent, levels = levels(data$ind))
  data <- na.omit(data)


  L <- data.frame(s = as.numeric(data$parent) - 1,
                  t = as.numeric(data$ind) - 1,
                  v = 1)
  s <- 's'
  t <- 't'
  v <- 'v'

  N <- data.frame(name = levels(data$ind),
                  size = 1,
                  grp = NA)
  nid <- 'name'
  nsize <- 'size'
  grp <- 'grp'
  logger$log('Create network data DONE')

  logger$log('Create network ...')

  supThisWarning({
    p <- networkD3::forceNetwork(L,N,s,t,v,nid,nsize,grp,
                                 zoom = TRUE,
                                 fontSize = 14,
                                 fontFamily = "sans-serif",
                                 opacity = 0.9,
                                 opacityNoHover = 0.8,
                                 charge = -50,
                                 arrows = TRUE)
  }, "It looks like Source/Target is not zero-indexed. This is required in JavaScript and so your plot may not render.")
  logger$log('Create network DONE')

  logger$log("DONE, return output.")
  return(p)
}






#' Draw a plotly graph of blups data
#'
#' X axis is the crosses, and Y axis the blups. The points are located at the
#' expected value and the error bar length is the standard deviation.
#'
#' @param blupDta list of data.frame of 4 columns: "ind1", "ind2", "blup_exp", "blup_var"
#' @param sorting method to sort the individuals (X axis) can be:
#'   - "asc": sort the BLUP expected value in ascending order (from left to right)
#'   - "dec": sort the BLUP expected value in decreasing order (from left to right)
#'   - any other value will sort the individuals in alphabetical order (from left to right)
#' @param y_axisName Name of the Y axis (default = `trait`)
#' @param errorBarInterval length of XX\% interval of interest represented by the error bars (default=0.95)
#' @param trait name of the trait to plot. This should be a name of the blupDta list. (optional if only one trait in `blupDta`, it will be set to the name of this trait)
#' @return plotly graph
plotBlup <- function(blupDta,
                     sorting = 'alpha',
                     y_axisName = NULL,
                     errorBarInterval = 0.95,
                     trait = NULL) {
  logger <- Logger$new("r-plotBlup()")

  ### Check input ----
  logger$log('Check inputs ...')

  if (!is.list(blupDta)) {
    stop('`blupDta` must be a list')
  }

  if (length(blupDta) == 1 && is.null(trait)) {
    trait <- names(blupDta)
  } else {
    if (!trait %in% names(blupDta)) {
      stop('`trait` must be a name of the list `blupDta`')
    }
  }
  blupDta <- blupDta[[trait]]

  expectedColumns <- c("ind1", "ind2", "blup_exp", "blup_var")
  if (!all(expectedColumns %in% colnames(blupDta))) {
    stop('`blupDta` must conntains the folowing columns: ',
         paste(expectedColumns, collapse = ', '))
  }
  if (nrow(blupDta) == 0) {
    stop('`blupDta` do not contain any row')
  }
  expectedSort <- c('alpha', 'asc', 'dec')
  if (!sorting %in% expectedSort) {
    logger$log('Sorting method not recognised, ',
            'x axis will be sorted alphabetically. Recognised methods are: ',
            paste(expectedSort, collapse = ', '))
  }
  logger$log('Check inputs DONE')


  # get cross' names ----
  blupDta$cross <- paste0(blupDta$ind1,
                          '_X_',
                          blupDta$ind2)

  # sort x axis values ----
  logger$log('sort x axis ...')
  if (sorting == 'asc') {
    blupDta$cross <- reorder(blupDta$cross, blupDta$blup_exp)
  } else if (sorting == 'dec') {
    blupDta$cross <- reorder(blupDta$cross, -blupDta$blup_exp)
  }
  logger$log('sort x axis DONE')



  # draw graph ----
  logger$log('draw plot ...')

  # Calculate the length of the error bar.
  # Calculation based on the value of `errorBarInterval`
  # (eg. 0.95 -> 95% of the data are included in this interval centered on the
  # expected value)
  quantileOfinterest <- (1 - errorBarInterval)/2
  mean <- blupDta$blup_exp
  sd <- sqrt(blupDta$blup_var)
  quantileMin <- qnorm(quantileOfinterest, mean, sd)
  quantileMax <- qnorm(1 - quantileOfinterest, mean, sd)
  errorBarHalfLenght <- (quantileMax - quantileMin)/2
  quantileMinName <- paste0("Quantile ", quantileOfinterest*100,'%')
  quantileMaxName <- paste0("Quantile ", (1-quantileOfinterest)*100,'%')

  # new values to show
  blupDta[, "Standard deviation"] <- sd
  blupDta[, quantileMinName] <- quantileMin
  blupDta[, quantileMaxName] <- quantileMax

  # update columns names
  colnames(blupDta)[colnames(blupDta) == "ind1"] <- "Parent 1"
  colnames(blupDta)[colnames(blupDta) == "ind2"] <- "Parent 2"
  colnames(blupDta)[colnames(blupDta) == "blup_exp"] <- "Expected value"
  colnames(blupDta)[colnames(blupDta) == "blup_var"] <- "Variance"
  colnames(blupDta)[colnames(blupDta) == "cross"] <- "Cross"

  tooltipValuesOrder <- c(
    "Cross",
    "Parent 1",
    "Parent 2",
    "Expected value",
    "Variance",
    "Standard deviation",
    quantileMinName,
    quantileMaxName)

  blupDta <- blupDta[, tooltipValuesOrder]


  # plot
  legend <- paste0('Progenies expected values\n',
                  '(Error bars represent a ',
                  errorBarInterval*100,
                  '% interval)')
  p <- plotly::plot_ly(
    type = 'scatter',
    mode = 'markers',
    name = legend,
    data = blupDta,
    x = ~ Cross,
    y = ~ `Expected value`,
    error_y = ~ list(
      name = "Error",
      array = errorBarHalfLenght,
      type = 'data',
      color = '#000000'),
    hoverinfo = 'text',
    text = apply(blupDta, 1, function(l) {
      paste(names(l), ":", l, collapse = "\n")
    })
  )

  if (is.null(y_axisName)) {
    y_axisName <- trait
  }

  p <- plotly::layout(p,
    showlegend=T,
    yaxis = list(title = y_axisName),
    xaxis = list(title = "Cross")
  )
  logger$log('draw plot DONE')

  # return output
  p

}
