# predict_once_meta.R
# Standalone script: load model, predict, write result, exit (frees RAM)
# Aim 1 only: any metastasis prediction
library(jsonlite)
library(caret)
library(ranger)

# ── Args ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript predict_once_meta.R <cancer_type> <input.json> <output.json>")

cancer_type <- args[1]
input_file  <- args[2]
output_file <- args[3]

model_dir <- Sys.getenv("MODEL_DIR", unset = "/app/models")

# ── Helpers ───────────────────────────────────────────────────────────────
convert_sex <- function(s) {
  if (s == "male") return("1")
  if (s == "female") return("2")
  return(s)
}

# ── Patient data builder ─────────────────────────────────────────────────
create_patient_data <- function(cancer_type, d) {
  
  if (cancer_type == "breast") {
    return(data.frame(
      AGE            = as.numeric(d$age),
      TNM_N_cat      = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat      = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4")),
      breast_er_cat  = factor(d$er_status, levels = c("No","Yes")),
      breast_pr_cat  = factor(d$pr_status, levels = c("No","Yes")),
      breast_her2_cat= factor(d$her2_status, levels = c("No","Yes")),
      breast_br_cat  = factor(d$grade, levels = c("Grade 1","Grade 2","Grade 3","Grade 4"))
    ))
  }
  
  if (cancer_type == "prostate") {
    return(data.frame(
      AGE            = as.numeric(d$age),
      TNM_N_cat      = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat      = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4")),
      prostate_psa   = as.numeric(d$psa),
      prostate_cores = as.numeric(d$core_ratio),
      prostate_gs_cat= factor(d$gleason, levels = c("Grade group 1","Grade group 2","Grade group 3","Grade group 4","Grade group 5"))
    ))
  }
  
  if (cancer_type == "colon") {
    return(data.frame(
      AGE       = as.numeric(d$age),
      SEX       = factor(convert_sex(d$sex), levels = c("1","2")),
      TNM_N_cat = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4")),
      colon_cea = as.numeric(d$cea)
    ))
  }
  
  if (cancer_type == "rectum") {
    return(data.frame(
      AGE            = as.numeric(d$age),
      SEX            = factor(convert_sex(d$sex), levels = c("1","2")),
      TNM_N_cat      = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat      = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4")),
      rectum_cea_cat = factor(d$cea_status, levels = c("Negative/normal; within normal limits","Positive/elevated"))
    ))
  }
  
  if (cancer_type == "urine") {
    return(data.frame(
      AGE       = as.numeric(d$age),
      SEX       = factor(convert_sex(d$sex), levels = c("1","2")),
      TNM_N_cat = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4"))
    ))
  }
  
  if (cancer_type == "esophagu") {
    return(data.frame(
      AGE       = as.numeric(d$age),
      SEX       = factor(convert_sex(d$sex), levels = c("1","2")),
      TNM_N_cat = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4"))
    ))
  }
  
  if (cancer_type == "melanoma") {
    return(data.frame(
      AGE       = as.numeric(d$age),
      SEX       = factor(convert_sex(d$sex), levels = c("1","2")),
      TNM_N_cat = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4"))
    ))
  }
  
  if (cancer_type == "liver") {
    return(data.frame(
      AGE           = as.numeric(d$age),
      SEX           = factor(convert_sex(d$sex), levels = c("1","2")),
      TNM_N_cat     = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat     = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4")),
      liver_afp_cat = factor(d$afp_status, levels = c("Negative/normal","Positive/elevated")),
      liver_fs_cat  = factor(d$fibrosis_score, levels = c("Fibrosis score 0-4","Fibrosis score 5-6"))
    ))
  }
  
  if (cancer_type == "kidney") {
    return(data.frame(
      AGE        = as.numeric(d$age),
      SEX        = factor(convert_sex(d$sex), levels = c("1","2")),
      TNM_N_cat  = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat  = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4")),
      rcc_sf_cat = factor(d$surgical_factors, levels = c("Yes","No")),
      rcc_fng_cat= factor(d$fuhrman_grade, levels = c("1","2","3","4"))
    ))
  }
  
  if (cancer_type == "ovary") {
    return(data.frame(
      AGE             = as.numeric(d$age),
      TNM_N_cat       = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat       = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4")),
      ovary_ca125_cat = factor(d$ca125_status, levels = c("Negative/normal; within normal limits","Positive/elevated"))
    ))
  }
  
  if (cancer_type == "retroper") {
    return(data.frame(
      AGE            = as.numeric(d$age),
      SEX            = factor(convert_sex(d$sex), levels = c("1","2")),
      TNM_N_cat      = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat      = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4")),
      retroper_sg_cat= factor(d$surgical_grade, levels = c("1","2","3"))
    ))
  }
  
  if (cancer_type == "testis") {
    return(data.frame(
      AGE            = as.numeric(d$age),
      TNM_N_cat      = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat      = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4")),
      testis_afp_cat = factor(d$afp_status, levels = c("Within normal limits","Elevated")),
      testis_hcg_cat = factor(d$hcg_status, levels = c("Within normal limits","Elevated")),
      testis_ldh_cat = factor(d$ldh_status, levels = c("Within normal limits","Elevated"))
    ))
  }
  
  if (cancer_type %in% c("LNSC","LSC")) {
    return(data.frame(
      AGE       = as.numeric(d$age),
      SEX       = factor(convert_sex(d$sex), levels = c("1","2")),
      TNM_N_cat = factor(d$tnm_n, levels = c("N0","N1","N2","N3")),
      TNM_T_cat = factor(d$tnm_t, levels = c("T0","T1","T2","T3","T4"))
    ))
  }
  
  stop(paste("Unsupported cancer type:", cancer_type))
}

# ── Main ──────────────────────────────────────────────────────────────────

inp <- fromJSON(input_file)

model_file <- file.path(model_dir, paste0(cancer_type, ".RData"))
if (!file.exists(model_file)) stop(paste("Model not found:", model_file))

# Build patient data once (tiny, stays in memory)
patient_data <- create_patient_data(cancer_type, inp)

# ── Step 1: Load logistic model, predict, then free it ──
model_env <- new.env()
load(model_file, envir = model_env)
logistic_model <- model_env$logistic_model
rm(model_env); gc()

logistic_pred <- predict(logistic_model, patient_data, type = "prob")
rm(logistic_model); gc()

# ── Step 2: Load tree model, predict, then free it ──
model_env2 <- new.env()
load(model_file, envir = model_env2)
tree_model <- model_env2$tree_model
rm(model_env2); gc()

tree_pred <- predict(tree_model, patient_data, type = "prob")
rm(tree_model); gc()

# Build result
result <- list(
  cancer_type = cancer_type,
  patient_data = inp,
  predictions = list(
    logistic_regression = list(
      no_metastasis = round(logistic_pred$No, 4),
      metastasis    = round(logistic_pred$Yes, 4),
      risk_level    = ifelse(logistic_pred$Yes > 0.5, "HIGH", "LOW")
    ),
    random_forest = list(
      no_metastasis = round(tree_pred$No, 4),
      metastasis    = round(tree_pred$Yes, 4),
      risk_level    = ifelse(tree_pred$Yes > 0.5, "HIGH", "LOW")
    )
  ),
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)

write(toJSON(result, auto_unbox = TRUE), output_file)