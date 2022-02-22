# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# function to check r-geno-tool-engine objects





checkRelMat <- function(relMat) {
  if (!is.numeric(relMat)) {
    stop('Bad relationship information, please be sure the relationship matrix file have been created using `r-geno-tool-engine`.')
  }
  cols <- sort(colnames(relMat))
  rows <- sort(row.names(relMat))
  if (!identical(cols, rows)) {
    stop('Bad relationship information, please be sure the relationship matrix file have been created using `r-geno-tool-engine`.')
  }
  if (!isSymmetric(relMat)) {
    stop('Bad relationship information, please be sure the relationship matrix file have been created using `r-geno-tool-engine`.')
  }
  TRUE
}
