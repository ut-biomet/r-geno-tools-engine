# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# function to check r-geno-tool-engine objects





checkRelMat <- function(relMat) {
  errMessage <- paste('Bad relationship information,',
                      'please be sure the relationship matrix file',
                      'have been created using `r-geno-tool-engine`.')

  if (!is.numeric(relMat)) {
    stop(errMessage, '`relMat` is not numeric.')
  }
  cols <- sort(colnames(relMat))
  rows <- sort(row.names(relMat))
  if (!identical(cols, rows)) {
    stop(errMessage, "`relMat`'s rows and columns names are different.")
  }
  # if (!isSymmetric(relMat)) {
  #   stop(errMessage, '`relMat` is not symmetric.')
  # }
  if (!all.equal.numeric(relMat, t(relMat))) {
    stop(errMessage, '`relMat` is not approximately symmetric.')
  }

  return(TRUE)
}
