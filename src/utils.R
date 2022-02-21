# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Utilities function of this engine


#' R6 class use to log messages in this engine's function
#'
logger <- R6::R6Class(
  "logger",
  lock_objects = FALSE,
  public = list(
    # Public Fields ####
    #' @field context [char] context of the log
    context = NULL,

    # Public Methods ####
    #' @description Create a new logger object.
    #' @param context [char] context of the log, (eg. inside which function)
    #' @return A new `logger` object.
    #' @examples
    #' mylogger <- logger$new(context = NULL)
    initialize = function(context = NULL){
      self$context <- context
    },
    #' @description Description of log
    #' @param ... [char] arguments passed to "paste" for message
    #' @param time [bool] display time in log message
    #' @param context [bool] display context in log message
    log = function(...,
                   time = TRUE,
                   context = TRUE){
      if (time) {
        time <- as.character(Sys.time())
      } else {
        time <- "\t"
      }
      if (context) {
        context <- paste0(" - ", self$context, ": ")
      } else {
        context <- "  "
      }
      message <- paste(time, context, paste(...), "\n", sep = "")
      cat(message)

    }
  )
)







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




#' Write engine documentation
#'
#' @param srcDir path of R sources folder (default "./src")
#' @param docDir path of documentation folder (default "./doc")
#'
#' @return NULL
#' @details Write engine's functions documentation in a README.md file located in `docDir`.
writeDoc <- function(srcDir = "./src",
                     docDir = "./doc"){
  stopifnot(dir.exists(srcDir))
  stopifnot(dir.exists(docDir))
  tmpDir <- file.path(docDir, ".tmp")
  dir.create(tmpDir)

  outFile <- file.path(docDir, "README.md")
  file.remove(outFile)
  file.create(outFile)

  sourcefiles <- list.files(srcDir, pattern = ".R$", full.names = TRUE)
  for (sourcefile in sourcefiles) {
    source_env = roxygen2::env_file(sourcefile)
    source(sourcefile, local = source_env)
    rd_blocks = roxygen2::parse_file(sourcefile, env = source_env)
    help_topics = roxygen2::roclet_process(roxygen2::rd_roclet(),
                                           rd_blocks,
                                           env = source_env,
                                           # env = NULL,
                                           dirname(sourcefile))
    rd_code = lapply(help_topics, format)

    for (funName in names(rd_code)) {
      writeLines(rd_code[[funName]],
                 con = file.path(tmpDir, funName))
      Rd2md::Rd2markdown(file.path(tmpDir, funName), outFile, append = TRUE)
      file.remove(file.path(tmpDir, funName))
    }
  }
  file.remove(tmpDir)

  NULL
}
