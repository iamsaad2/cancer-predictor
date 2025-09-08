# Cancer Metastasis Prediction API using Plumber
# Lazy loading approach to handle memory constraints

library(plumber)
library(caret)
library(dplyr)
library(pROC)
library(ranger)
library(jsonlite)

#* @filter cors
cors <- function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "GET,POST,PUT,DELETE,OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type")
  
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$status <- 200
    return(list())
  }
  
  plumber::forward()
}

# Global variables - all models are now in the same directory
MODELS_DIR <- "./models" # Current directory where cancer_api.R is located

# Available cancer types
CANCER_TYPES <- c(
  "breast", "prostate", "colon", "rectum", "urine", "esophagu", 
  "melanoma", "liver", "kidney", "ovary", "retroper", "testis",
  "LNSC", "LSC"  # Lung cancer types
)

# Function to get model directory - simplified since all models are in same folder
get_model_dir <- function(cancer_type) {
  return(MODELS_DIR)
}

# Function to load specific cancer model
load_cancer_model <- function(cancer_type) {
  model_dir <- get_model_dir(cancer_type)
  model_file <- file.path(model_dir, paste0(cancer_type, ".RData"))
  
  if (!file.exists(model_file)) {
    stop(paste("Model file not found:", model_file))
  }
  
  # Create new environment to load models
  model_env <- new.env()
  load(model_file, envir = model_env)
  
  return(list(
    logistic = model_env$logistic_model,
    tree = model_env$tree_model
  ))
}

# Helper function to convert sex values
convert_sex_to_numeric <- function(sex_value) {
  if (sex_value == "male") return("1")
  if (sex_value == "female") return("2")
  return(sex_value)  # fallback for existing numeric values
}

# Complete create_patient_data function
create_patient_data <- function(cancer_type, input_data) {
  
  if (cancer_type == "breast") {
    return(data.frame(
      AGE = as.numeric(input_data$age),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4")),
      breast_er_cat = factor(input_data$er_status, levels = c("No", "Yes")),
      breast_pr_cat = factor(input_data$pr_status, levels = c("No", "Yes")),
      breast_her2_cat = factor(input_data$her2_status, levels = c("No", "Yes")),
      breast_br_cat = factor(input_data$grade, levels = c("Grade 1", "Grade 2", "Grade 3", "Grade 4"))
    ))
  }
  
  if (cancer_type == "prostate") {
    return(data.frame(
      AGE = as.numeric(input_data$age),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4")),
      prostate_psa = as.numeric(input_data$psa),
      prostate_cores = as.numeric(input_data$core_ratio),
      prostate_gs_cat = factor(input_data$gleason, levels = c("Grade group 1", "Grade group 2", "Grade group 3", "Grade group 4", "Grade group 5"))
    ))
  }
  
  if (cancer_type == "colon") {
    sex_numeric <- convert_sex_to_numeric(input_data$sex)
    return(data.frame(
      AGE = as.numeric(input_data$age),
      SEX = factor(sex_numeric, levels = c("1", "2")),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4")),
      colon_cea = as.numeric(input_data$cea)
    ))
  }
  
  if (cancer_type == "rectum") {
    sex_numeric <- convert_sex_to_numeric(input_data$sex)
    return(data.frame(
      AGE = as.numeric(input_data$age),
      SEX = factor(sex_numeric, levels = c("1", "2")),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4")),
      rectum_cea_cat = factor(input_data$cea_status, levels = c("Negative/normal; within normal limits", "Positive/elevated"))
    ))
  }
  
  if (cancer_type == "urine") {
    sex_numeric <- convert_sex_to_numeric(input_data$sex)
    return(data.frame(
      AGE = as.numeric(input_data$age),
      SEX = factor(sex_numeric, levels = c("1", "2")),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4"))
    ))
  }
  
  if (cancer_type == "esophagu") {
    sex_numeric <- convert_sex_to_numeric(input_data$sex)
    return(data.frame(
      AGE = as.numeric(input_data$age),
      SEX = factor(sex_numeric, levels = c("1", "2")),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4"))
    ))
  }
  
  if (cancer_type == "melanoma") {
    sex_numeric <- convert_sex_to_numeric(input_data$sex)
    return(data.frame(
      AGE = as.numeric(input_data$age),
      SEX = factor(sex_numeric, levels = c("1", "2")),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4"))
    ))
  }
  
  if (cancer_type == "liver") {
    sex_numeric <- convert_sex_to_numeric(input_data$sex)
    return(data.frame(
      AGE = as.numeric(input_data$age),
      SEX = factor(sex_numeric, levels = c("1", "2")),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4")),
      liver_afp_cat = factor(input_data$afp_status, levels = c("Negative/normal", "Positive/elevated")),
      liver_fs_cat = factor(input_data$fibrosis_score, levels = c("Fibrosis score 0-4", "Fibrosis score 5-6"))
    ))
  }
  
  if (cancer_type == "kidney") {
    sex_numeric <- convert_sex_to_numeric(input_data$sex)
    return(data.frame(
      AGE = as.numeric(input_data$age),
      SEX = factor(sex_numeric, levels = c("1", "2")),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4")),
      rcc_sf_cat = factor(input_data$surgical_factors, levels = c("Yes", "No")),
      rcc_fng_cat = factor(input_data$fuhrman_grade, levels = c("1", "2", "3", "4"))
    ))
  }
  
  if (cancer_type == "ovary") {
    return(data.frame(
      AGE = as.numeric(input_data$age),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4")),
      ovary_ca125_cat = factor(input_data$ca125_status, levels = c("Negative/normal; within normal limits", "Positive/elevated"))
    ))
  }
  
  if (cancer_type == "retroper") {
    sex_numeric <- convert_sex_to_numeric(input_data$sex)
    return(data.frame(
      AGE = as.numeric(input_data$age),
      SEX = factor(sex_numeric, levels = c("1", "2")),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4")),
      retroper_sg_cat = factor(input_data$surgical_grade, levels = c("1", "2", "3"))
    ))
  }
  
  if (cancer_type == "testis") {
    return(data.frame(
      AGE = as.numeric(input_data$age),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4")),
      testis_afp_cat = factor(input_data$afp_status, levels = c("Within normal limits", "Elevated")),
      testis_hcg_cat = factor(input_data$hcg_status, levels = c("Within normal limits", "Elevated")),
      testis_ldh_cat = factor(input_data$ldh_status, levels = c("Within normal limits", "Elevated"))
    ))
  }
  
  if (cancer_type %in% c("LNSC", "LSC")) {
    sex_numeric <- convert_sex_to_numeric(input_data$sex)
    return(data.frame(
      AGE = as.numeric(input_data$age),
      SEX = factor(sex_numeric, levels = c("1", "2")),
      TNM_N_cat = factor(input_data$tnm_n, levels = c("N0", "N1", "N2", "N3")),
      TNM_T_cat = factor(input_data$tnm_t, levels = c("T0", "T1", "T2", "T3", "T4"))
    ))
  }
  
  # Default fallback
  stop(paste("Cancer type not supported:", cancer_type))
}

#* @apiTitle Cancer Metastasis Prediction API
#* @apiDescription API for predicting cancer metastasis risk

#* Get available cancer types
#* @get /cancer-types
function() {
  return(list(cancer_types = CANCER_TYPES))
}

#* Get required inputs for a specific cancer type
#* @param cancer_type Cancer type
#* @get /inputs/<cancer_type>
function(cancer_type) {
  
  if (!cancer_type %in% CANCER_TYPES) {
    return(list(error = "Invalid cancer type"))
  }
  
  base_inputs <- list(
    age = "Patient age (number)",
    tnm_n = "Lymph node status (N0, N1, N2, N3)",
    tnm_t = "Tumor status (T0, T1, T2, T3, T4)"
  )
  
  if (cancer_type == "breast") {
    return(c(base_inputs, list(
      er_status = "ER status (No, Yes)",
      pr_status = "PR status (No, Yes)", 
      her2_status = "HER2 status (No, Yes)",
      grade = "Grade (Grade 1, Grade 2, Grade 3, Grade 4)"
    )))
  }
  
  if (cancer_type == "prostate") {
    return(c(base_inputs, list(
      psa = "PSA level (number)",
      core_ratio = "Core ratio (number 0-1)",
      gleason = "Gleason score (Grade group 1, Grade group 2, Grade group 3, Grade group 4, Grade group 5)"
    )))
  }
  
  if (cancer_type == "colon") {
    return(c(base_inputs, list(
      sex = "Sex (male, female)",
      cea = "CEA level (number)"
    )))
  }
  
  if (cancer_type == "rectum") {
    return(c(base_inputs, list(
      sex = "Sex (male, female)",
      cea_status = "CEA status (Negative/normal; within normal limits, Positive/elevated)"
    )))
  }
  
  if (cancer_type == "urine") {
    return(c(base_inputs, list(
      sex = "Sex (male, female)"
    )))
  }
  
  if (cancer_type == "esophagu") {
    return(c(base_inputs, list(
      sex = "Sex (male, female)"
    )))
  }
  
  if (cancer_type == "melanoma") {
    return(c(base_inputs, list(
      sex = "Sex (male, female)"
    )))
  }
  
  if (cancer_type == "liver") {
    return(c(base_inputs, list(
      sex = "Sex (male, female)",
      afp_status = "AFP status (Negative/normal, Positive/elevated)",
      fibrosis_score = "Fibrosis score (Fibrosis score 0-4, Fibrosis score 5-6)"
    )))
  }
  
  if (cancer_type == "kidney") {
    return(c(base_inputs, list(
      sex = "Sex (male, female)",
      surgical_factors = "Surgical factors (Yes, No)",
      fuhrman_grade = "Fuhrman grade (1, 2, 3, 4)"
    )))
  }
  
  if (cancer_type == "ovary") {
    return(c(base_inputs, list(
      ca125_status = "CA-125 status (Negative/normal; within normal limits, Positive/elevated)"
    )))
  }
  
  if (cancer_type == "retroper") {
    return(c(base_inputs, list(
      sex = "Sex (male, female)",
      surgical_grade = "Surgical grade (1, 2, 3)"
    )))
  }
  
  if (cancer_type == "testis") {
    return(c(base_inputs, list(
      afp_status = "AFP status (Within normal limits, Elevated)",
      hcg_status = "HCG status (Within normal limits, Elevated)",
      ldh_status = "LDH status (Within normal limits, Elevated)"
    )))
  }
  
  if (cancer_type %in% c("LNSC", "LSC")) {
    return(c(base_inputs, list(
      sex = "Sex (male, female)"
    )))
  }
  
  # Default fallback
  return(c(base_inputs, list(
    sex = "Sex (male, female)"
  )))
}

#* Predict metastasis risk
#* @param cancer_type Cancer type
#* @post /predict/<cancer_type>
function(req, cancer_type) {
  
  # Validate cancer type
  if (!cancer_type %in% CANCER_TYPES) {
    return(list(error = "Invalid cancer type"))
  }
  
  # Parse request body
  input_data <- jsonlite::fromJSON(req$postBody)
  
  tryCatch({
    # Load models for this cancer type
    models <- load_cancer_model(cancer_type)
    
    # Create patient data
    patient_data <- create_patient_data(cancer_type, input_data)
    
    # Make predictions
    logistic_pred <- predict(models$logistic, patient_data, type = "prob")
    tree_pred <- predict(models$tree, patient_data, type = "prob")
    
    # Format response
    result <- list(
      cancer_type = cancer_type,
      patient_data = input_data,
      predictions = list(
        logistic_regression = list(
          no_metastasis = round(logistic_pred$No, 4),
          metastasis = round(logistic_pred$Yes, 4),
          risk_level = ifelse(logistic_pred$Yes > 0.5, "HIGH", "LOW")
        ),
        random_forest = list(
          no_metastasis = round(tree_pred$No, 4),
          metastasis = round(tree_pred$Yes, 4),
          risk_level = ifelse(tree_pred$Yes > 0.5, "HIGH", "LOW")
        )
      ),
      timestamp = Sys.time()
    )
    
    # Clean up memory
    rm(models)
    gc()
    
    return(result)
    
  }, error = function(e) {
    return(list(error = paste("Prediction failed:", e$message)))
  })
}

#* Health check endpoint
#* @get /health
function() {
  return(list(
    status = "healthy",
    timestamp = Sys.time(),
    available_cancer_types = length(CANCER_TYPES)
  ))
}

#* Get example request for a cancer type
#* @param cancer_type Cancer type
#* @get /example/<cancer_type>
function(cancer_type) {
  
  if (cancer_type == "breast") {
    return(list(
      example_request = list(
        age = 55,
        tnm_n = "N1",
        tnm_t = "T2", 
        er_status = "Yes",
        pr_status = "Yes",
        her2_status = "No",
        grade = "Grade 2"
      )
    ))
  }
  
  if (cancer_type == "kidney") {
    return(list(
      example_request = list(
        age = 56,
        sex = "male",
        tnm_n = "N0",
        tnm_t = "T1",
        surgical_factors = "Yes",
        fuhrman_grade = "2"
      )
    ))
  }
  
  if (cancer_type == "LNSC") {
    return(list(
      example_request = list(
        age = 65,
        sex = "male",
        tnm_n = "N1",
        tnm_t = "T2"
      )
    ))
  }
  
  # Default example
  return(list(
    example_request = list(
      age = 60,
      sex = "male",
      tnm_n = "N1",
      tnm_t = "T2"
    )
  ))
}