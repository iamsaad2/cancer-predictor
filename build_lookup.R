# build_lookup.R
# Enumerate input combinations per cancer using each model's actual training
# levels, predict, write CSV lookup tables.
# Usage: Rscript build_lookup.R [cancer_type | all]   (defaults to "all")

suppressPackageStartupMessages({
  library(caret)
  library(ranger)
})

args   <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[1] else "all"

model_dir <- "models"
out_dir   <- "lookup_tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

NUMERIC_VALUES <- list(
  age   = 0:99,
  psa   = c(2, 7, 15, 35, 75, 150),     # PSA bucket midpoints
  cores = c(0.1, 0.3, 0.5, 0.7, 0.9),   # core_ratio bucket midpoints
  cea   = c(1, 3.75, 7.5, 15, 35, 100)  # CEA bucket midpoints
)

# Per-cancer column order + type. "cat" = pulled from model$xlevels.
SPECS <- list(
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
  LNSC     = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  LSC      = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  uterine  = list(AGE="age", TNM_N_cat="cat", TNM_T_cat="cat"),
  thyroid  = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat", thyroid_cat="cat"),
  pancreas = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  stomach  = list(AGE="age", SEX="cat", TNM_N_cat="cat", TNM_T_cat="cat"),
  cervix   = list(AGE="age", TNM_N_cat="cat", TNM_T_cat="cat")
)

build_grid <- function(model, cancer) {
  spec <- SPECS[[cancer]]
  cols <- list()
  for (var in names(spec)) {
    kind <- spec[[var]]
    if (kind == "cat") {
      lv <- model$xlevels[[var]]
      if (is.null(lv)) stop(sprintf("[%s] model has no xlevels for %s", cancer, var))
      cols[[var]] <- factor(lv, levels = lv)
    } else {
      vals <- NUMERIC_VALUES[[kind]]
      if (is.null(vals)) stop(sprintf("[%s] unknown numeric kind '%s' for %s", cancer, kind, var))
      cols[[var]] <- vals
    }
  }
  do.call(expand.grid, cols)
}

run_one <- function(cancer) {
  model_file <- file.path(model_dir, paste0(cancer, "_tree.rds"))
  if (!file.exists(model_file)) {
    cat(sprintf("[%s] SKIP — model file %s not found\n", cancer, model_file))
    return(invisible(NULL))
  }
  cat(sprintf("[%s] loading model...\n", cancer))
  model <- readRDS(model_file)

  grid <- build_grid(model, cancer)
  cat(sprintf("[%s] %d combinations — predicting...\n", cancer, nrow(grid)))

  t0   <- Sys.time()
  pred <- predict(model, grid, type = "prob")
  dt   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[%s] predict() took %.1f sec\n", cancer, dt))

  out <- cbind(grid,
               no_metastasis = round(pred$No,  4),
               metastasis    = round(pred$Yes, 4),
               risk_level    = ifelse(pred$Yes > 0.5, "HIGH", "LOW"))

  out_file <- file.path(out_dir, paste0(cancer, "_lookup.csv.gz"))
  gz <- gzfile(out_file, "w")
  write.csv(out, gz, row.names = FALSE)
  close(gz)
  sz_mb <- file.info(out_file)$size / 1024 / 1024
  cat(sprintf("[%s] wrote %s (%.2f MB, %d rows)\n", cancer, out_file, sz_mb, nrow(out)))

  rm(model, grid, pred, out); gc(verbose = FALSE)
  invisible(NULL)
}

ALL_CANCERS <- c(
  "urine","esophagu","melanoma","ovary","retroper","testis","LNSC","LSC",
  "rectum","colon","liver","kidney","breast","prostate",
  "uterine","thyroid","pancreas","stomach","cervix"
)

if (target == "all") {
  for (c in ALL_CANCERS) run_one(c)
} else {
  run_one(target)
}
