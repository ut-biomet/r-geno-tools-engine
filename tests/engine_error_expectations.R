# Author: Julien Diot juliendiot@ut-biomet.org
# 2024 The University of Tokyo
#
# Description:
# engineError expectation

expect_engineError <- function(code) {
  error <- expect_error(
    {
      code
    },
    class = "engineError"
  )

  if (is.null(error$extra$code)) {
    return(invisible(error))
  }

  expected_fields <- list(
    "BAD_ARG_ADJ_METHOD" = c("expected", "provided"),
    "BAD_ARG_CHR" = c("expected", "provided"),
    "BAD_ARG_CONFIDENCELEVEL" = c("expected", "provided"),
    "BAD_ARG_FILTER_NPOINTS" = c("expected", "provided"),
    "BAD_ARG_FILTER_PADJ" = c("expected", "provided"),
    "BAD_ARG_FILTER_QUANT" = c("expected", "provided"),
    "BAD_ARG_FIXED" = c("expected", "provided", "test"),
    "BAD_ARG_FROM_TO" = c("expected", "provided"),
    "BAD_ARG_INTERACTIVE" = c("expected", "provided"),
    "BAD_ARG_METHOD" = c("expected", "provided"),
    "BAD_ARG_NCROSS" = c("expected", "provided"),
    "BAD_ARG_RESPONSE" = c("expected", "provided"),
    "BAD_ARG_SUFFIX" = c("expected", "provided"),
    "BAD_ARG_TAU_OMEGA" = c("expected", "provided"),
    "BAD_ARG_TEST" = c("expected", "provided"),
    "BAD_ARG_THRESH_CALLRATE" = c("expected", "provided"),
    "BAD_ARG_THRESH_CALLRATE_FILTER_OUT_ALL_MARKERS" = c("provided", "max_geno_callrate"),
    "BAD_ARG_THRESH_MAF" = c("expected", "provided"),
    "BAD_ARG_THRESH_MAF_FILTER_OUT_ALL_MARKERS" = c("provided", "max_geno_maf"),
    "BAD_ARG_THRESH_P" = c("expected", "provided"),
    "BAD_ARG_TRAIT" = c("expected", "provided"),
    "BAD_CROSSTABLE_COLUMNS" = c("expected", "provided"),
    "BAD_CROSSTABLE_EMPTY" = c("expected", "provided"),
    "BAD_CROSSTABLE_MISSING_VALUES" = c("expected", "provided", "n_missing"),
    "BAD_CROSSTABLE_UNAVAILABLE_INDS" = c("n_unavailable_ind", "unavailable_inds"),
    "BAD_GENOTYPE_UNPHASED" = c(),
    "BAD_GENOTYPE_DUPLICATED_SNP_IDS" = c("dup_snp_chrom", "dup_snp_pos", "dup_snp_ids"),
    "BAD_GENO_ALL_MONOMORPHIC_SNP" = c(),
    "BAD_GENO_PHENO_NO_COMMON_INDS" = c(),
    "BAD_GWAS_FILE" = c(),
    "BAD_PROGENY_BLUP_FILE" = c("expected", "provided"),
    "BAD_MARKER_EFFECTS_1ST_COLUMN" = c("expected", "provided"),
    "BAD_MARKER_EFFECTS_DUPLICATED_ID" = c("n_duplicated_ids", "duplicated_ids"),
    "BAD_MARKER_EFFECTS_FORMAT_EMPTY" = c(),
    "BAD_MARKER_EFFECTS_MISSING_VALUES" = c(),
    "BAD_MARKER_EFFECTS_JSON_KEYS" = c("expected", "provided"),
    "BAD_OUTPUT_FILE_FORMAT" = c("expected", "provided"),
    "BAD_PEDIGREE_FORMAT_EMPTY" = c(),
    "BAD_PEDIGREE_FORMAT_INCONSISTENT_GENEALOGY" = c("n_inconsistent", "affected_individuals"),
    "BAD_PEDIGREE_FORMAT_INCONSISTENT_ID" = c("n_inconsistent", "affected_lines"),
    "BAD_PEDIGREE_FORMAT_MISSING_IND_ID" = c("n_missing"),
    "BAD_PEDIGREE_FORMAT_N_COLUMNS" = c("expected", "provided"),
    "BAD_PHENOTYPE_DUPLICATES" = c("n_duplicates", "duplicates"),
    "BAD_PHENO_DATA_CLASS" = c("columnWithWrongTypes", "detectedClass"),
    "BAD_SNPCOORD_DUPLICATED_ID" = c("n_duplicated_ids", "duplicated_ids"),
    "BAD_SNPCOORD_MISSING_LINK_MAP_POS_COLUMN" = c("expected", "provided"),
    "BAD_SNPCOORD_MISSING_SNPID_COLUMNS" = c(),
    "BAD_SNPCOORD_MISSING_VALUE" = c("columns", "expected", "provided"),
    "BAD_SNPCOORD_MISSING_SNP" = c("n_missing_SNP", "n_provided_SNP", "n_expected_SNP", "missing_SNP", "reference"),
    "BAD_SNPCOORD_INCONSISTENT_SNP_ORDER" = c("affected_chr"),
    "BAD_SNPCOORD_SNP_LINKMAP_POSITION_NOT_NUMERIC" = c("n_snp", "n_not_numeric_linkMapPos", "not_numeric_linkMapPos_ids"),
    "BAD_SNPCOORD_PHYSICAL_POS_VCF_MISSMATCH" = c("n_physPos_missmatch", "n_snp", "physPos_missmatch_ids"),
    "INCONSISTENT_RELATIONSHIP_MATRICES" = c(),
    "DOMINANCE_MODEL_NOT_APPLICABLE" = c("n_homozygous_ind", "n_homozygous_snp", "n_ind", "n_snp", "homozygous_threshold"),
    "DOMINANCE_MODEL_FOR_BLUP_ESTIMATION" = c(),
    "FILE_EXIST" = c("file"),
    "FILE_NOT_FOUND" = c("file"),
    "INPUT_FILE_NOT_PROVIDED" = c("input_file"),
    "BAD_FILE_FORMAT" = c("expected", "provided", "file"),
    "MULTIPLE_INPUTFILE" = c("expected", "provided"),
    "MULTIPLE_OUTFILE" = c("expected", "provided"),
    "NO_OUTPUT_DIR" = c("dir")
  )

  if (!error$extra$code %in% names(expected_fields)) {
    stop("extra info expectations for code `", error$extra$code, "` unknown.")
  }

  if (is.null(UNTESTED_ERROR_CODES)) {
    UNTESTED_ERROR_CODES <<- names(expected_fields)
  } else {
    UNTESTED_ERROR_CODES <<- UNTESTED_ERROR_CODES[UNTESTED_ERROR_CODES != error$extra$code]
  }

  for (expected_field in expected_fields[[error$extra$code]]) {
    expect(
      expected_field %in% names(error$extra),
      sprintf("engineError with code `%s` miss extra info: `%s`.", error$extra$code, expected_field)
    )
  }

  engine_error_default_fields <- c("code", "error_location", "trace")
  for (provided_field in names(error$extra)) {
    if (!provided_field %in% c(
      expected_fields[[error$extra$code]],
      engine_error_default_fields
    )) {
      warning(sprintf(
        "engineError with code `%s` have unexpected field: `%s`.",
        error$extra$code,
        provided_field
      ))
    }
  }

  invisible(error)
}
