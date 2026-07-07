# build_lookup_aim2.R
# Generate CSV lookup tables for Aim 2 (site-specific metastasis prediction).
# Uses random-forest (tree_model) only. Inputs: Aim 1 covariates + mets_* (Yes/No)
# for the 3 sites OTHER than the target site, + mate_other (Yes/No).
# Usage: Rscript build_lookup_aim2.R [<cancer>_<site> | all]
#        Models loaded from aim2_models/<file>.RData

suppressPackageStartupMessages({
  library(caret)
  library(ranger)
})

args   <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[1] else "all"

model_dir <- "aim2_models"
out_dir   <- "lookup_tables_aim2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Age bucket midpoints (10 buckets covering 0-99)
AGE_BUCKETS <- seq(5, 95, by = 10)

# Numeric input buckets (re-used from Aim 1 build_lookup.R)
PSA_BUCKETS   <- c(2, 7, 15, 35, 75, 150)
CORES_BUCKETS <- c(0.1, 0.3, 0.5, 0.7, 0.9)
CEA_BUCKETS   <- c(1, 3.75, 7.5, 15, 35, 100)

CANCERS <- c("breast","prostate","colon","rectum","urine","esophagu",
             "melanoma","liver","kidney","ovary","retroper","testis","lnsc","lsc",
             "uterine","thyroid","pancreas","stomach","cervix")
SITES   <- c("bone","brain","liver","lung")

# File naming: most cancers are "<cancer>_<site>.RData"; prostate is "prostate_decision_<site>.RData".
model_file <- function(cancer, site) {
  if (cancer == "prostate") {
    return(file.path(model_dir, sprintf("prostate_decision_%s.RData", site)))
  }
  file.path(model_dir, sprintf("%s_%s.RData", cancer, site))
}

# Per-cancer column order + numeric/categorical kind. "cat" pulls levels from model$xlevels.
# Mirrors Aim 1 SPECS but with mets_* (Yes/No) added per model below.
AIM1_SPECS <- list(
  breast   = list(AGE="age", TNM_N_cat="cat", TNM_T_cat="cat",
                  breast_er_cat="cat", breast_pr_cat="cat", breast_her2_cat="cat", breast_br_cat="cat"),
  prostate = list(AGE="age", TNM_N_cat="cat", TNM_T_cat="cat",
                  prostate_psa="psa", prostate_cores="cores", prostate_gs_cat="cat"),
  colon    = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat", colon_cea="cea"),
  rectum   = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat", rectum_cea_cat="cat"),
  urine    = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  esophagu = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  melanoma = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  liver    = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat",
                  liver_afp_cat="cat", liver_fs_cat="cat"),
  kidney   = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat",
                  rcc_sf_cat="cat", rcc_fng_cat="cat"),
  ovary    = list(AGE="age", TNM_N_cat="cat", TNM_T_cat="cat", ovary_ca125_cat="cat"),
  retroper = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat", retroper_sg_cat="cat"),
  testis   = list(AGE="age", TNM_N_cat="cat", TNM_T_cat="cat",
                  testis_afp_cat="cat", testis_hcg_cat="cat", testis_ldh_cat="cat"),
  lnsc     = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  lsc      = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  # New 5 cancers. Note: Aim 2 thyroid does NOT use thyroid_cat (Aim 1 only).
  uterine  = list(AGE="age", TNM_N_cat="cat", TNM_T_cat="cat"),
  thyroid  = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  pancreas = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  stomach  = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  cervix   = list(AGE="age", TNM_N_cat="cat", TNM_T_cat="cat")
)

NUMERIC_VALUES <- list(
  age   = AGE_BUCKETS,
  psa   = PSA_BUCKETS,
  cores = CORES_BUCKETS,
  cea   = CEA_BUCKETS
)

YES_NO <- function() factor(c("No","Yes"), levels = c("No","Yes","Unknown"))

# Build expand.grid for one (cancer, site) pair. Includes mets_* for the 3 OTHER sites
# plus mate_other, all Yes/No only.
build_grid <- function(model, cancer, site) {
  spec <- AIM1_SPECS[[cancer]]
  cols <- list()
  for (var in names(spec)) {
    kind <- spec[[var]]
    if (kind == "cat") {
      lv <- model$xlevels[[var]]
      if (is.null(lv)) stop(sprintf("[%s_%s] model has no xlevels for %s", cancer, site, var))
      cols[[var]] <- factor(lv, levels = lv)
    } else {
      cols[[var]] <- NUMERIC_VALUES[[kind]]
    }
  }
  # Add mets_* inputs (Yes/No only — skip "Unknown" per design)
  other_sites <- setdiff(SITES, site)
  for (s in other_sites) {
    v <- paste0("mets_", s)
    if (is.null(model$xlevels[[v]])) {
      stop(sprintf("[%s_%s] model has no xlevels for %s", cancer, site, v))
    }
    # Preserve the model's full factor levels (so prediction sees a recognized factor)
    cols[[v]] <- factor(c("No","Yes"), levels = model$xlevels[[v]])
  }
  # mate_other (note: model uses the typo'd name)
  if (!is.null(model$xlevels[["mate_other"]])) {
    cols[["mate_other"]] <- factor(c("No","Yes"), levels = model$xlevels[["mate_other"]])
  }
  do.call(expand.grid, cols)
}

run_one <- function(cancer, site) {
  tag <- sprintf("%s_%s", cancer, site)
  f   <- model_file(cancer, site)
  if (!file.exists(f)) {
    cat(sprintf("[%s] SKIP — model file %s not found\n", tag, f))
    return(invisible(NULL))
  }
  cat(sprintf("[%s] loading %s\n", tag, f))
  e <- new.env()
  load(f, envir = e)
  model <- e$tree_model
  if (is.null(model)) stop(sprintf("[%s] no tree_model in %s", tag, f))
  rm(e); gc(verbose = FALSE)

  grid <- build_grid(model, cancer, site)
  cat(sprintf("[%s] %d combinations — predicting...\n", tag, nrow(grid)))

  t0   <- Sys.time()
  pred <- predict(model, grid, type = "prob")
  dt   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[%s] predict() took %.1f sec\n", tag, dt))

  out <- cbind(grid,
               p_mets = round(pred$Yes, 4))

  out_file <- file.path(out_dir, sprintf("%s_%s_lookup.csv.gz", cancer, site))
  gz <- gzfile(out_file, "w")
  write.csv(out, gz, row.names = FALSE)
  close(gz)
  sz_mb <- file.info(out_file)$size / 1024 / 1024
  cat(sprintf("[%s] wrote %s (%.2f MB, %d rows)\n", tag, out_file, sz_mb, nrow(out)))

  rm(model, grid, pred, out); gc(verbose = FALSE)
  invisible(NULL)
}

# Targets
if (target == "all") {
  for (c in CANCERS) for (s in SITES) run_one(c, s)
} else {
  parts <- strsplit(target, "_")[[1]]
  if (length(parts) != 2) stop("Argument must be '<cancer>_<site>' or 'all'")
  run_one(parts[1], parts[2])
}
