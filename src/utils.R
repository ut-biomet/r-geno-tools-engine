# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Utilities function of this engine


#' Do not show a specific warning message
#'
#' @param expr expression to evaluate.
#' @param warnMessage warning message to catch
#'
#' @details This is based on the `base::suppressWarnings` function.
#'
supThisWarning <- function(expr, warnMessage) {
  withCallingHandlers(expr, warning = function(warn) {
    if (identical(warn$message, warnMessage)) {
      tryInvokeRestart("muffleWarning")
    }
  })
}
