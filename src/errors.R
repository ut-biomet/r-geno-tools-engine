#' Title
#'
#' @param message error message
#' @param extra list of extra information
#' @param n_skip_caller (int, default 1) This is to catch the function where the
#' error happens. 0 will show this function, 1 will show the function calling this one and so on.
#'
engineError <- function(message, extra = list(), n_skip_caller = 1) {
  error_locations <- c()
  i <- n_skip_caller
  last_caller <- rlang::caller_call(i)
  while (!is.null(last_caller)) {
    error_locations <- c(error_locations, last_caller)
    i <- i + 1
    last_caller <- rlang::caller_call(i)
  }
  error_locations <- as.character(error_locations)
  last_location <- error_locations[1]
  error_locations <- rev(error_locations)
  error_locations <- paste(error_locations, collapse = " -> ")
  print(error_locations)

  message <- paste0("Error in ", last_location, ": ", message,"\nAdditional information:\n",
                    paste0("- ", unlist(lapply(names(extra), function(field) {
                      paste(field, ":", paste(as.character(extra[[field]]), collapse = ", "), collapse = "\n")
                    })), collapse = "\n")
  )

  err <- structure(
    list(
      message = message,
      extra = c(
        extra,
        list("error_location" = last_location,
             "trace" = error_locations)
      ),
      call = NULL
    ),
    class = c(extra$code, "engineError", "error", "condition")
  )
  stop(err)
}


#' Raise error
#'
#' To be used inside functions to check their arguments.
#' the error message will be:
#' "`arg` must be `must_be` not `not`"
#' (adapted from https://adv-r.hadley.nz/conditions.html#signalling)
#'
#' @param arg (character) tested argument name
#' @param must_be (character) expected value or type
#' @param not (any) provided value
#' @param errType "type" is the only recognised value. Use it to report an
#' error about the type of the argument
#' @param class use "engineError" (default) to raise an expected "engineError"
#' @param extra (list, default NULL) extra information to pass to engineError
#' @param n_skip_caller (int, default 2) see `engineError`
#'
#' @examples
#' x = "a"
#' bad_argument("x", must_be = 42, not = x)
#' # will stop with the error msg: "`x` must be 42 not `a`"
#' bad_argument("x", must_be = "numeric", not = x, "type")
#' # will stop with the error msg: "`x` must be numeric not `character`"
bad_argument = function(arg, must_be, not = NULL, errType = "value", class = "engineError", extra = NULL, n_skip_caller = 2) {
  must_be_msg <- must_be
  if (length(must_be) > 1) {
    must_be_msg <- paste0("`", paste(must_be, collapse = ", "), "`")
  }

  if (errType == "type") {
    not = typeof(not)
  }

  not_msg <- ""
  if (!is.null(not)) {
    if (length(not) >= 1) {
      not <- paste0("`", paste(not, collapse = ", "), "`")
    }
    not_msg = paste0(" not ", not)
  }

  msg <- paste0(arg, " must be ", must_be_msg, not_msg)
  extra <- c(list("expected" = must_be,
                  "provided" = not),
             extra)
  if (identical(class, "engineError")) {
    engineError(msg, extra = extra, n_skip_caller = n_skip_caller)
  }
  stop(msg)
}

#' list of error codes as function so that unexpected error code raise and error
errorCode <- function(code) {
  codes <- list(
    "BAD_ARG_ADJ_METHOD",
    "BAD_ARG_CHR",
    "BAD_ARG_CONFIDENCELEVEL",
    "BAD_ARG_FILTER_NPOINTS",
    "BAD_ARG_FILTER_PADJ",
    "BAD_ARG_FILTER_QUANT",
    "BAD_ARG_FIXED",
    "BAD_ARG_FROM_TO",
    "BAD_ARG_INTERACTIVE",
    "BAD_ARG_METHOD",
    "BAD_ARG_NCROSS",
    "BAD_ARG_RESPONSE",
    "BAD_ARG_SUFFIX",
    "BAD_ARG_TAU_OMEGA",
    "BAD_ARG_TEST",
    "BAD_ARG_THRESH_CALLRATE",
    "BAD_ARG_THRESH_MAF",
    "BAD_ARG_THRESH_P",
    "BAD_ARG_TRAIT",

    "BAD_CROSSTABLE_COLUMNS",
    "BAD_CROSSTABLE_EMPTY",
    "BAD_CROSSTABLE_EMPTY",
    "BAD_CROSSTABLE_MISSING_VALUES",
    "BAD_CROSSTABLE_UNAVAILABLE_INDS",

    "BAD_GENOTYPE_UNPHASED",

    "BAD_GWAS_FILE",
    "BAD_PROGENY_BLUP_FILE",

    "BAD_MARKER_EFFECTS_1ST_COLUMN",
    "BAD_MARKER_EFFECTS_DUPLICATED_ID",
    "BAD_MARKER_EFFECTS_FORMAT_EMPTY",
    "BAD_MARKER_EFFECTS_MISSING_VALUES",

    "BAD_OUTPUT_FILE_FORMAT",

    "BAD_PEDIGREE_FORMAT_EMPTY",
    "BAD_PEDIGREE_FORMAT_INCONSISTENT_GENEALOGY",
    "BAD_PEDIGREE_FORMAT_INCONSISTENT_ID",
    "BAD_PEDIGREE_FORMAT_MISSING_IND_ID",
    "BAD_PEDIGREE_FORMAT_N_COLUMNS",

    "BAD_PHENOTYPE_DUPLICATES",

    "BAD_SNPCOORD_DUPLICATED_ID",
    "BAD_SNPCOORD_MISSING_LINK_MAP_POS_COLUMN",
    "BAD_SNPCOORD_MISSING_SNPID_COLUMNS",
    "BAD_SNPCOORD_MISSING_VALUE",
    "BAD_SNPCOORD_MISSING_SNP",
    "BAD_SNPCOORD_INCONSISTENT_SNP_ORDER",

    "INCONSISTENT_RELATIONSHIP_MATRICES",

    "FILE_EXIST",
    "FILE_NOT_FOUND",
    "BAD_FILE_FORMAT",
    "MULTIPLE_INPUTFILE",
    "MULTIPLE_OUTFILE",
    "NO_OUTPUT_DIR"

  )
  names(codes) <- codes

  if (is.null(codes[[code]])) {
    err <- structure(
      list(
        message = paste0("Error code `", code, "` does not exist."),
      ),
      class = c("missingErrorCode", "error", "condition")
    )
    stop(err)
  }
  return(codes[[code]])
}
