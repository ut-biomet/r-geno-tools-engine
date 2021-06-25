#' R6 class use to log messages
#'
#' @import R6
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
    #' # create specie:
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
                   context = TRUE,
                   status = NULL){
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
