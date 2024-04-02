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
    bad_argument('relMat',
                 must_be = "numeric matrix",
                 not = relMat,
                 "type",
                 class = NULL)
  }
  cols <- sort(colnames(relMat))
  rows <- sort(row.names(relMat))
  if (!identical(cols, rows)) {
    bad_argument("colnames(relMat) and row.names(relMat)",
                 must_be = "identical",
                 class = NULL)
  }
  if (!all.equal.numeric(relMat, t(relMat))) {
    bad_argument("relMat",
                 must_be = "approximately symmetric",
                 class = NULL)
  }
  invisible(TRUE)
}


#' Check compatibility between 2 snps coordinates data set
#' and keep only genotypes SNPs
#'
#' @param user_SNPcoord SNPs coordinates data coming from the user (`.csv` file)
#' @param vcf_SNPcoord SNPs coordinates data coming from the `.vcf` file
#'
#' @description If the `user_SNPcoord` data miss some SNPs defined in the `.vcf`
#' file, an error is raised. If the `user_SNPcoord` data have additional SNPs,
#' those SNPs will be removed.
#' If the data between those two data-set are inconsistent, an error is raised.
#'
#' @return The filtered `user_SNPcoord` data.frame
checkAndFilterSNPcoord <- function(user_SNPcoord, vcf_SNPcoord) {
  # SNP ids
  user_SNPid <- user_SNPcoord$SNPid
  vcf_SNPid <- vcf_SNPcoord$SNPid
  missingSNP <- which(!vcf_SNPid %in% user_SNPid)
  if (length(missingSNP) != 0) {
    msg <- paste('The SNPs coordinate file miss', length(missingSNP),
                 'SNPs when compared with the provided vcf file.')
    engineError(msg,
                extra = list(
                  "code" = errorCode("BAD_SNPCOORD_MISSING_SNP"),
                  "n_missing_SNP" = length(missingSNP),
                  "n_provided_SNP" = length(user_SNPid),
                  "n_expected_SNP" = length(vcf_SNPid),
                  "missing_SNP" = vcf_SNPid[missingSNP],
                  "reference" = "vcf file"
                ))
  }

  additionalSNP <- which(!user_SNPid %in% vcf_SNPid)
  if (length(additionalSNP) != 0) {
    msg <- paste('The SNPs coordinate file have', length(additionalSNP),
                 'SNPs not defined in the `.vcf` file,',
                 'those SNP will not be considered:',
                 paste(user_SNPid[additionalSNP], collapse = ', '))
    warning(msg)
    user_SNPcoord <- user_SNPcoord[-additionalSNP,]
  }

  inconsistent_order_chr <- unlist(lapply(unique(user_SNPcoord$chr), function(chr) {
    snpcoord <- user_SNPcoord[user_SNPcoord$chr == chr,]
    snpcoord <- snpcoord[order(snpcoord$physPos),]
    if (is.unsorted(snpcoord$linkMapPos, strictly = TRUE)){
      return(chr)
    }
    return(NULL)
  }))

  if (length(inconsistent_order_chr) != 0) {
    engineError("SNP's position order should be similar when sorted by physical position and by linkage map position",
                extra = list(
                  "code" = errorCode("BAD_SNPCOORD_INCONSISTENT_SNP_ORDER"),
                  "affected_chr" = inconsistent_order_chr
                ))
  }


  # add vcf_file snp physical position to SNPcoord data frame
  user_SNPcoord[order(user_SNPcoord$SNPid), 'physPos'] <- vcf_SNPcoord[order(vcf_SNPcoord$SNPid), 'physPos']

  return(user_SNPcoord)
}



#' Check individuals in the crossing table are in the haplotype data
#'
#' @param crossTable the crossing table
#' @param haplo the haplotype data given by the function `readPhasedGeno`
#'
#' @return NULL, raise error if missing individuals are detected.
checkIndNamesConsistency <- function(crossTable, haplo) {
  crossTableInds <- unique(c(crossTable$ind1, crossTable$ind2))
  haploInds <- gsub('_[12]$', '', colnames(haplo)[1:(ncol(haplo)/2)])

  missIndsId <- which(!crossTableInds %in% haploInds)
  if (length(missIndsId) != 0) {
    msg <- paste(
      length(missIndsId),
      'individuals are defined in the crossing table',
      'but not in the genoype file:',
      paste(crossTableInds[missIndsId], collapse = ', ')
    )
    engineError(msg, extra = list(
      code = errorCode("BAD_CROSSTABLE_UNAVAILABLE_INDS"),
      n_unavailable_ind = length(missIndsId),
      unavailable_inds = crossTableInds[missIndsId]
    ))
  }
  invisible(NULL)
}


check_test <- function(test){
  code <- "BAD_ARG_TEST"

  if (!is.character(test)) {
    bad_argument("test",
                 must_be = 'a character string of length 1',
                 not = test,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(test) != 1) {
    bad_argument("length(test)",
                 must_be = '1',
                 not = length(test),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  accepted_tests <- c("score", "wald", "lrt")
  if (!test %in% accepted_tests) {
    bad_argument(
      "test",
      must_be = paste("one of ", paste0(accepted_tests, collapse = ", ")),
      not = test,
      n_skip_caller = 3,
      extra = list(code = errorCode(code)))
  }
  invisible(TRUE)
}

check_fixed <- function(fixed, test) {
  code <- "BAD_ARG_FIXED"

  if (test %in% c("wald", "lrt")) {
    if (length(fixed) != 1) {
      bad_argument("length(fixed)",
                   must_be = '1 when `test` is "wald" or "lrt"',
                   not = length(fixed),
                   n_skip_caller = 3,
                   extra = list(code = errorCode(code),
                                test = test))
    }
    if (!is.numeric(fixed)) {
      bad_argument("fixed",
                   must_be = 'a positive integer when `test` is "wald" or "lrt"',
                   not = fixed,
                   "type",
                   n_skip_caller = 3,
                   extra = list(code = errorCode(code),
                                test = test))
    }
    if (fixed < 0 || fixed %% 1 != 0) {
      bad_argument("fixed",
                   must_be = 'a positive integer when `test` is "wald" or "lrt"',
                   not = fixed,
                   n_skip_caller = 3,
                   extra = list(code = errorCode(code),
                                test = test))
    }
  } else if (!(is.null(fixed) | fixed == 0)) {
    warning('`fixed` seems to have been set but it will not be taken in account when `test` is "',
            test, '"')
  }
  invisible(TRUE)
}

check_response <- function(response) {
  code <- "BAD_ARG_RESPONSE"

  if (!is.character(response)) {
    bad_argument("response",
                 must_be = 'a character string of length 1',
                 not = response, "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(response) != 1) {
    bad_argument("length(response)",
                 must_be = '1',
                 not = length(response),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  accepted_responses <- c("quantitative", "binary")
  if (!response %in% accepted_responses) {
    bad_argument(
      "response",
      must_be = paste("one of ", paste0(accepted_responses, collapse = ", ")),
      not = response,
      n_skip_caller = 3,
      extra = list(code = errorCode(code)))
  }
  invisible(TRUE)
}

check_thresh_maf <- function(thresh_maf, geno = NULL) {
  code <- "BAD_ARG_THRESH_MAF"

  if (!is.numeric(thresh_maf)) {
    bad_argument("thresh_maf",
                 must_be = "a numeric between 0 and 0.5",
                 not = thresh_maf,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(thresh_maf) != 1) {
    bad_argument("length(thresh_maf)",
                 must_be = "`1`",
                 not = length(thresh_maf),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (thresh_maf < 0 || thresh_maf > 0.5) {
    bad_argument("thresh_maf",
                 must_be = "a numeric between 0 and 0.5",
                 not = thresh_maf,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (!is.null(geno)) {
    if (!any(geno@snps$maf > thresh_maf)) {
      engineError("The provided `thresh_maf` is too high and filter out all SNPs.",
                  n_skip_caller = 3,
                  extra = list(code = errorCode("BAD_ARG_THRESH_MAF_FILTER_OUT_ALL_MARKERS"),
                               provided = thresh_maf,
                               max_geno_maf = max(geno@snps$maf)))
    }
  }
  invisible(TRUE)
}

check_thresh_callrate <- function(thresh_callrate, geno = NULL) {
  code <- "BAD_ARG_THRESH_CALLRATE"

  if (!is.numeric(thresh_callrate)) {
    bad_argument("thresh_callrate",
                 must_be = "a numeric between 0 and 1",
                 not = thresh_callrate,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(thresh_callrate) != 1) {
    bad_argument("length(thresh_callrate)",
                 must_be = "`1`",
                 not = length(thresh_callrate),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (thresh_callrate < 0 || thresh_callrate > 1) {
    bad_argument("thresh_callrate",
                 must_be = "a numeric between 0 and 1",
                 not = thresh_callrate,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  if (!is.null(geno)) {
    if (!any(geno@snps$callrate > thresh_callrate)) {
      engineError("The provided `thresh_callrate` is too high and filter out all SNPs.",
                  n_skip_caller = 3,
                  extra = list(code = errorCode("BAD_ARG_THRESH_CALLRATE_FILTER_OUT_ALL_MARKERS"),
                               provided = thresh_callrate,
                               max_geno_collrate = max(geno@snps$callrate)))
    }
  }
  invisible(TRUE)
}

check_file_extention <- function(file, expected_exts) {
  code <- "BAD_FILE_FORMAT"

  f <- file
  full_ext <- tools::file_ext(f)
  while (full_ext != "" && tools::file_ext(tools::file_path_sans_ext(f)) != "") {
    f <- tools::file_path_sans_ext(f)
    full_ext <- paste(tools::file_ext(f), full_ext, sep = ".")
  }

  expected_exts <- gsub("^\\.", "", expected_exts)

  if (!full_ext %in% expected_exts) {
    bad_argument('file extention',
                 must_be = paste0("one of [", paste0("`", expected_exts, "`", collapse = ", "), "]"),
                 not = full_ext,
                 extra = list(file = file,
                              code = errorCode(code)),
                 n_skip_caller = 4)
  }
  invisible(TRUE)
}


check_outFile <- function(outFile, expected_exts = NULL, allow_existing = TRUE, accept_null = FALSE) {

  if (accept_null && is.null(outFile)) {
    return(invisible(TRUE))
  }

  if (length(outFile) != 1) {
    bad_argument("length(outFile)",
                 must_be = "1" ,
                 not = length(outFile),
                 n_skip_caller = 3,
                 extra = list(code = errorCode("MULTIPLE_OUTFILE")))
  }
  if (file.exists(outFile)) {
    if (allow_existing) {
      warning('output "file" already exists. This file will be overwritten.')
    } else {
      bad_argument("outFile",
                   must_be = "innexistant",
                   n_skip_caller = 3,
                   extra = list(code = errorCode("FILE_EXIST"),
                                file = outFile))
    }
  }

  if (!is.null(expected_exts)) {
    check_file_extention(outFile, expected_exts)
  }

  invisible(TRUE)
}



check_inputFile <- function(file, expected_exts = NULL) {
  if (length(file) != 1) {
    bad_argument("length(file)",
                 must_be = "1" ,
                 not = length(file),
                 n_skip_caller = 3,
                 extra = list(code = errorCode("MULTIPLE_INPUTFILE")))
  }
  if (!file.exists(file)) {
    engineError(
      "file not found",
      extra = list(code = errorCode("FILE_NOT_FOUND"), file = file),
      n_skip_caller = 3
    )
  }
  if (!is.null(expected_exts)) {
    check_file_extention(file, expected_exts)
  }
  invisible(TRUE)
}


check_trait <- function(trait, available_traits) {
  code <- "BAD_ARG_TRAIT"
  if (!is.character(trait)) {
    bad_argument("trait",
                 must_be = 'character string of length 1',
                 not = trait,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(trait) != 1) {
    bad_argument("length(trait)",
                 must_be = '1',
                 not = length(trait),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (!trait %in% available_traits) {
    bad_argument("trait",
                 must_be = paste("one of ", paste(available_traits, collapse = ", ")),
                 not = trait,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  invisible(TRUE)
}

check_adj_method <- function(adj_method) {
  code <- "BAD_ARG_ADJ_METHOD"
  if (!is.character(adj_method)) {
    bad_argument("adj_method",
                 must_be = 'a character string of length 1',
                 not = adj_method,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(adj_method) != 1) {
    bad_argument("length(adj_method)",
                 must_be = '1',
                 not = length(adj_method),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  accetped_adj_methods <- c("holm", "hochberg","bonferroni", "BH", "BY", "fdr", "none")
  if (!adj_method %in% accetped_adj_methods) {
    bad_argument('adj_method',
      must_be = paste("one of ", paste0(accetped_adj_methods, collapse = ", ")),
      not = adj_method,
      n_skip_caller = 3,
      extra = list(code = errorCode(code)))
  }
  invisible(TRUE)
}


check_thresh_p <- function(thresh_p) {
  code <- "BAD_ARG_THRESH_P"
  if (!is.numeric(thresh_p)) {
    bad_argument('thresh_p',
                 must_be = "a numeric",
                 not = thresh_p, "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(thresh_p) != 1) {
    bad_argument("length(thresh_p)",
                 must_be = '1',
                 not = length(thresh_p),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (thresh_p <= 0 || thresh_p >= 1) {
    bad_argument("thresh_p",
                 must_be = 'between 0 and 1',
                 not = thresh_p,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  invisible(TRUE)
}

check_interactive <- function(interactive) {
  code <- "BAD_ARG_INTERACTIVE"
  if (!is.logical(interactive)) {
    bad_argument('interactive',
                 must_be = "a boolean",
                 not = interactive, "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(interactive) != 1) {
    bad_argument("length(interactive)",
                 must_be = '1',
                 not = length(interactive),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  invisible(TRUE)
}


check_filter_pAdj <- function(filter_pAdj) {
  code <- "BAD_ARG_FILTER_PADJ"
  if (!is.numeric(filter_pAdj)) {
    bad_argument("filter_pAdj",
                 must_be = "a numeric",
                 not = filter_pAdj,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(filter_pAdj) != 1) {
    bad_argument("length(filter_pAdj)",
                 must_be = '1',
                 not = length(filter_pAdj),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (filter_pAdj < 0 || filter_pAdj > 1) {
    bad_argument("filter_pAdj",
                 must_be = 'between [0, 1]',
                 not = filter_pAdj,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  invisible(TRUE)
}


check_filter_nPoints <- function(filter_nPoints) {
  code <- "BAD_ARG_FILTER_NPOINTS"
  if (!is.numeric(filter_nPoints)) {
    bad_argument("filter_nPoints",
                 must_be = "integer",
                 not = filter_nPoints,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(filter_nPoints) != 1) {
    bad_argument("length(filter_nPoints)",
                 must_be = '1',
                 not = length(filter_nPoints),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (filter_nPoints != Inf) {
    if (filter_nPoints <= 0 || filter_nPoints %% 1 != 0 ) {
      bad_argument("filter_nPoints",
                   must_be = 'a positive integer',
                   not = filter_nPoints,
                   n_skip_caller = 3,
                   extra = list(code = errorCode(code)))
    }
  }
  invisible(TRUE)
}

check_filter_quant <-  function(check_filter_quant) {
  code <- "BAD_ARG_FILTER_QUANT"
  if (!is.numeric(check_filter_quant)) {
    bad_argument("check_filter_quant",
                 must_be = "a numeric",
                 not = check_filter_quant,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(check_filter_quant) != 1) {
    bad_argument("length(check_filter_quant)",
                 must_be = '1',
                 not = length(check_filter_quant),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (check_filter_quant < 0 || check_filter_quant > 1) {
    bad_argument("check_filter_quant",
                 must_be = 'between [0, 1]',
                 not = check_filter_quant,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  invisible(TRUE)
}





check_chr <- function(chr, chr_list) {
  code <- "BAD_ARG_CHR"
  if (!is.na(chr)) {
    if (length(chr) != 1) {
      bad_argument("length(chr)",
                   must_be = '1',
                   not = length(chr),
                   n_skip_caller = 3,
                   extra = list(code = errorCode(code)))
    }
    if (!chr %in% unique(chr_list)) {
      bad_argument("chr",
                   must_be = paste("one of ", paste0(unique(chr_list), collapse = ", ")),
                   not = chr,
                   n_skip_caller = 3,
                   extra = list(code = errorCode(code)))
    }
  }
  invisible(TRUE)
}



check_from_to <- function(from, to, max_diff) {
  code <- "BAD_ARG_FROM_TO"
  # from
  if (!is.numeric(from)) {
    bad_argument("from",
                 must_be = 'a positive integer',
                 not = from,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(from) != 1) {
    bad_argument("length(from)",
                 must_be = '1',
                 not = length(from),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (from < 0 || from %% 1 != 0) {
    bad_argument("from",
                 must_be = 'a positive integer',
                 not = from,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  # to
  if (!is.numeric(to)) {
    bad_argument("to",
                 must_be = 'a positive integer',
                 not = to,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(to) != 1) {
    bad_argument("length(to)",
                 must_be = '1',
                 not = length(to),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (to < 0 || to %% 1 != 0) {
    bad_argument("to",
                 must_be = 'a positive integer',
                 not = from,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  # from to
  if (from >= to) {
    bad_argument("from",
                 must_be = paste('lower or equal than', to, '(value of `to`)'),
                 not = from,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  if (to - from > max_diff) {
    bad_argument('to - from',
                 must_be = paste("lower or equal than", max_diff) ,
                 not = to - from,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  invisible(TRUE)
}


check_method <- function(method, tau, omega) {
  code <- "BAD_ARG_METHOD"
  if (length(method) != 1) {
    bad_argument("length(method)", must_be = "1", not = length(method),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  expected_methods <- c('Martini', 'Legarra')
  if (!method %in% expected_methods) {
    bad_argument("method",
                 must_be = paste("one of ", paste(expected_methods, collapse = ", ")),
                 not = method,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  if (method != "Martini") {
    if (!is.null(tau) || !is.null(omega)) {
      warnMsg <- paste('`tau` and `omega` are only used for "Martini" method.',
                       'Those parameters will be ignored.')
      warning(warnMsg)
    }
  }
  invisible(TRUE)
}


check_tau_omega <- function(tau, omega) {
  code <- "BAD_ARG_TAU_OMEGA"
  if (tau < 0) {
    bad_argument("tau",
                 must_be = ">=0",
                 not = tau,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (omega > 1) {
    bad_argument("omega",
                 must_be = "<=1",
                 not = omega,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (identical(c(tau, omega), c(0, 1))) {
    bad_argument("c(tau, omega)",
                 must_be = "different that c(0,1)",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  invisible(TRUE)
}


check_nCross <- function(nCross) {
  code <- "BAD_ARG_NCROSS"
  if (!is.numeric(nCross)) {
    bad_argument('thresh_p',
                 must_be = "a numeric",
                 not = nCross, "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(nCross) != 1) {
    bad_argument("length(nCross)",
                 must_be = '1',
                 not = length(nCross),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (nCross <= 0 || nCross %% 1 != 0) {
    bad_argument("nCross",
                 must_be = 'a strictly positive integer',
                 not = nCross,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }

  invisible(TRUE)
}



check_sorting <- function(sorting) {
  expectedSort <- c('alpha', 'asc', 'dec')
  if (!sorting %in% expectedSort) {
    warning('Sorting method not recognised, ',
            'x axis will be sorted alphabetically. Recognised methods are: ',
            paste(expectedSort, collapse = ', '))
  }
  invisible(TRUE)
}


check_confidenceLevel <- function(confidenceLevel) {
  code <- "BAD_ARG_CONFIDENCELEVEL"
  if (!is.numeric(confidenceLevel)) {
    bad_argument('confidenceLevel',
                 must_be = "a numeric",
                 not = confidenceLevel, "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(confidenceLevel) != 1) {
    bad_argument("length(confidenceLevel)",
                 must_be = '1',
                 not = length(confidenceLevel),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (confidenceLevel <= 0 || confidenceLevel >= 1) {
    bad_argument("confidenceLevel",
                 must_be = 'between 0 and 1',
                 not = confidenceLevel,
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  invisible(TRUE)
}

check_suffix <- function(suffix) {
  code <- "BAD_ARG_SUFFIX"
  if (!is.character(suffix)) {
    bad_argument("suffix",
                 must_be = 'character string of length 1',
                 not = suffix,
                 "type",
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  if (length(suffix) != 1) {
    bad_argument("length(suffix)",
                 must_be = '1',
                 not = length(suffix),
                 n_skip_caller = 3,
                 extra = list(code = errorCode(code)))
  }
  invisible(TRUE)
}
