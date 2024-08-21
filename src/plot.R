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

  logger$log("Remove NAs ...")
  # remove NA lines form data (important for cases where all the data had
  # been removed at the gwas step)
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
    bad_argument('relMat', must_be = "a matrix" , not = relMat, "type")
  }
  if (!isSymmetric(relMat)) {
    bad_argument('relMat', must_be = "a symetric matrix")
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
    bad_argument('ped', must_be = "a list generated by `readPedData` function", not = ped, "type")
  }
  expected_ped_names <- c("data", "graph")
  if (!identical(names(ped), expected_ped_names)) {
    bad_argument(
      'names(ped)',
      must_be = paste("one of ", paste(expected_ped_names, sep = ", ")),
      not = names(ped)
    )
  }
  if (!is.data.frame(ped$data)) {
    bad_argument('ped$data', must_be = "a data.frame", not = ped$data, "type")
  }
  expected_ped_colnames <- c('ind', 'parent1', 'parent2')
  if (!identical(colnames(ped$data), expected_ped_colnames)) {
    bad_argument(
      'colnames(ped$data)',
      must_be = expected_ped_colnames,
      not = colnames(ped$data)
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






#' Draw a plotly graph of blups data for 1 trait
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
plotBlup_1trait <- function(blupDta,
                            sorting = 'alpha',
                            y_axisName = NULL,
                            errorBarInterval = 0.95,
                            trait) {
  logger <- Logger$new("r-plotBlup_1trait()")

  blupDta <- do.call(rbind, lapply(blupDta, function(blupDta_cross){
    data.frame(
      ind1 = blupDta_cross$ind1,
      ind2 = blupDta_cross$ind2,
      blup_exp = blupDta_cross$blup_exp[[trait]],
      blup_var = blupDta_cross$blup_var[[trait]]
    )
  }))
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




#' Draw a plotly graph of blups data for 2 traits
#'
#' The points are located at the expected value and the ellipses
#' size represent the `confidenceLevel` prediction interval.
#'
#' @param blupDta list returned by calc_progenyBlupEstimation
#' @param x_trait name of the trait to show on the x axis
#' @param y_trait name of the trait to show on the y axis
#' @param confidenceLevel level of the prediction ellipses (default 0.95, ie 95%
#' ellypses)
#' @param x_suffix suffix to add to the x axis's name
#' @param y_suffix suffix to add to the y axis's name
#' @param ellipses_npoints number of points used to draw the ellipses (default 100)
#' @return plotly graph
plotBlup_2traits <- function(blupDta,
                             x_trait,
                             y_trait,
                             confidenceLevel = 0.95,
                             x_suffix = "",
                             y_suffix = "",
                             ellipses_npoints = 100) {

  logger <- Logger$new("r-plotBlup_2traits()")

  ### Draw plot ----
  logger$log('draw plot ...')
  p <- plotly::plot_ly(type = "scatter",
               mode = "markers",
               colors = "Set3")
  for (blupDta_cross in blupDta) {
    cov <- as.matrix(blupDta_cross$cov[c(x_trait, y_trait),
                                       c(x_trait, y_trait)])
    center <- c(blupDta_cross$blup_exp[[x_trait]],
                blupDta_cross$blup_exp[[y_trait]])
    cross <- paste0(blupDta_cross$ind1, " X ", blupDta_cross$ind2)

    ellipse_dta <- ellipse::ellipse(cov,
                                    centre = center,
                                    level = confidenceLevel,
                                    npoints = ellipses_npoints)

    ellipse_dta <- data.frame(
      x = unlist(ellipse_dta[, 1]),
      y = unlist(ellipse_dta[, 2])
    )
    colnames(ellipse_dta) <- c(x_trait, y_trait)

    center <- as.data.frame(t(center))
    colnames(center) <- paste0("expeced value ", c(x_trait, y_trait))
    center$`parent 1` <- blupDta_cross$ind1
    center$`parent 2` <- blupDta_cross$ind2
    center$cross <- cross
    center[, paste0("variance ", x_trait)] <- blupDta_cross$blup_var[x_trait]
    center[, paste0("variance ", y_trait)] <- blupDta_cross$blup_var[y_trait]
    center[, paste0("covariance ", x_trait, "/", y_trait)] <- blupDta_cross$cov[x_trait, y_trait]

    ellipse_dta$cross <- cross

    p <- plotly::add_markers(
      p,
      data = center,
      x = center[, paste0("expeced value ", x_trait)],
      y = center[, paste0("expeced value ", y_trait)],
      color = ~cross,
      hoverinfo = "text",
      text = apply(center, 1, function(l) {
        paste(names(l), ":", l, collapse = "\n")
      })
    )
    p <- plotly::add_trace(
      p,
      type = "scatter",
      mode = "lines",
      inherit = FALSE,
      data = ellipse_dta,
      x = ellipse_dta[, x_trait],
      y = ellipse_dta[, y_trait],
      color = ~cross,
      hoverinfo = "text",
      text = paste0(100 * confidenceLevel,
                    "% ellipse\nCross: ",
                    ellipse_dta$cross)
    )
  }
  p <- plotly::layout(
    p,
    xaxis = list(title = paste(x_trait, x_suffix)),
    yaxis = list(title = paste(y_trait, y_suffix))
  )

  logger$log('draw plot DONE')
  p
}




#' Draw a plotly graph of a GS model cross-validation evaluation
#'
#' This plots is composed of several subplots:
#' - Observed vs Predicted scatter plot for all the cross-validation folds
#' - Horizontal box plots of models metrics calculated during the
#' cross-validation
#'
#' @param evaluation_results list returned by `cross_validation_evaluation()`
#'
#' @return plotly graph
evaluation_plot <- function(evaluation_results) {

  pred_dta <- evaluation_results$predictions

  n_reps <- length(unique(pred_dta$repetition))
  n_folds <- length(unique(pred_dta$fold)) / n_reps

  scatter_plot <- plotly::plot_ly(type = "scatter", mode = "markers")

  min_x <- min(pred_dta$actual)
  max_x <- max(pred_dta$actual)
  scatter_plot <- plotly::add_lines(scatter_plot,
                                    inherit = FALSE,
                                    x = c(min_x, max_x),
                                    y = c(min_x, max_x),
                                    line = list(color = 'red', dash = 'dash'),
                                    name = "identity line"
  )

  scatter_plot <- plotly::add_markers(
    scatter_plot,
    data = pred_dta,
    x = ~actual,
    y = ~predicted,
    color = ~repetition,
    hoverinfo = "text",
    text = apply(pred_dta, 1, function(l) {
      paste(names(l), ":", l, collapse = "\n")
    })
  )

  scatter_plot <- plotly::layout(
    scatter_plot,
    yaxis = list(title = "Predicted values"),
    xaxis = list(title = "Observed values")
  )

  all_metrics <- evaluation_results$metrics

  metrics_long_names <- list(
    "rmse" = "Root Mean Square Error",
    "corel_pearson" = "Correlation (Pearson)",
    "corel_spearman" = "Correlation (Spearman)",
    "r2" = "R squared"
  )

  metrics <- colnames(all_metrics)[colnames(all_metrics) %in% names(metrics_long_names)]

  box_plots <- lapply(metrics, function(metric) {
    mean_val <- signif(mean(all_metrics[, metric], na.rm = TRUE), 3)
    box_plot <- plotly::plot_ly(
      type = "box",
      x = all_metrics[, metric],
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      hoverinfo = "text",
      name = metrics_long_names[[metric]],
      text = apply(all_metrics, 1, function(l) {
        paste(names(l), ":", l, collapse = "\n")
      })
    )
    box_plot <- plotly::layout(
      box_plot,
      yaxis = list(title = "", showticklabels = FALSE),
      xaxis = list(
        title = list(
          standoff = 0,
          text = paste0(metrics_long_names[[metric]], "\nMean value: ", mean_val)
        )
      )
    )
    box_plot %>% plotly::layout(title = paste0("Cross-Validation Results (",
                                               n_reps, " repetitions, ", n_folds, " folds)"),
                                plot_bgcolor = '#e5ecf6'
                                )
  })


  margin <- c(0.0, 0.0, 0.12, 0.0)
  heights = rep(1/length(box_plots), length(box_plots))
  # plotly::subplot is a bit "buggy" and because the margins are not applied on
  # the 1st and last plot, those one appears bigger than the others
  # we need to manually tweak the heights of the plots to get equal sizes
  heights[1] <- heights[1] - margin[3]
  heights <- heights + margin[3]/length(heights)

  box_plots <- plotly::subplot(box_plots,
                               nrows = length(box_plots),
                               titleY = TRUE,
                               titleX = TRUE,
                               heights = heights,
                               margin = margin)

  fig <- plotly::subplot(scatter_plot,
                         box_plots,
                         widths = c(0.7, 0.3),
                         titleY = TRUE,
                         titleX = TRUE)
  fig <- plotly::layout(fig,
                 title = paste0("Cross-Validation Results (",
                                n_reps, " repetitions, ", n_folds, " folds)"),
                 plot_bgcolor = '#e5ecf6')

  return(fig)
}
