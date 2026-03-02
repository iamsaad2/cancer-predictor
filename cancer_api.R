# cancer_api.R  – Lightweight API: spawns separate Rscript per prediction to keep RAM low
library(plumber)
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

# Available cancer types
CANCER_TYPES <- c(
  "breast", "prostate", "colon", "rectum", "urine", "esophagu",
  "melanoma", "liver", "kidney", "ovary", "retroper", "testis",
  "LNSC", "LSC"
)

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
    age   = "Patient age (number)",
    tnm_n = "Lymph node status (N0, N1, N2, N3)",
    tnm_t = "Tumor status (T0, T1, T2, T3, T4)"
  )
  
  extra <- switch(cancer_type,
                  breast = list(
                    er_status   = "ER status (No, Yes)",
                    pr_status   = "PR status (No, Yes)",
                    her2_status = "HER2 status (No, Yes)",
                    grade       = "Grade (Grade 1, Grade 2, Grade 3, Grade 4)"
                  ),
                  prostate = list(
                    psa        = "PSA level (number)",
                    core_ratio = "Core ratio (number 0-1)",
                    gleason    = "Gleason score (Grade group 1, Grade group 2, Grade group 3, Grade group 4, Grade group 5)"
                  ),
                  colon = list(
                    sex = "Sex (male, female)",
                    cea = "CEA level (number)"
                  ),
                  rectum = list(
                    sex        = "Sex (male, female)",
                    cea_status = "CEA status (Negative/normal; within normal limits, Positive/elevated)"
                  ),
                  urine = , esophagu = , melanoma = list(
                    sex = "Sex (male, female)"
                  ),
                  liver = list(
                    sex            = "Sex (male, female)",
                    afp_status     = "AFP status (Negative/normal, Positive/elevated)",
                    fibrosis_score = "Fibrosis score (Fibrosis score 0-4, Fibrosis score 5-6)"
                  ),
                  kidney = list(
                    sex              = "Sex (male, female)",
                    surgical_factors = "Surgical factors (Yes, No)",
                    fuhrman_grade    = "Fuhrman grade (1, 2, 3, 4)"
                  ),
                  ovary = list(
                    ca125_status = "CA-125 status (Negative/normal; within normal limits, Positive/elevated)"
                  ),
                  retroper = list(
                    sex            = "Sex (male, female)",
                    surgical_grade = "Surgical grade (1, 2, 3)"
                  ),
                  testis = list(
                    afp_status = "AFP status (Within normal limits, Elevated)",
                    hcg_status = "HCG status (Within normal limits, Elevated)",
                    ldh_status = "LDH status (Within normal limits, Elevated)"
                  ),
                  LNSC = , LSC = list(
                    sex = "Sex (male, female)"
                  ),
                  list(sex = "Sex (male, female)")
  )
  
  return(c(base_inputs, extra))
}

#* Predict metastasis risk (spawns separate process to keep RAM low)
#* @param cancer_type Cancer type
#* @post /predict/<cancer_type>
function(req, res, cancer_type) {
  
  if (!cancer_type %in% CANCER_TYPES) {
    res$status <- 400
    return(list(error = "Invalid cancer type"))
  }
  
  inp <- tryCatch(jsonlite::fromJSON(req$postBody), error = function(e) NULL)
  if (is.null(inp)) {
    res$status <- 400
    return(list(error = "Invalid JSON body"))
  }
  
  # Write input to temp file
  input_file  <- tempfile(fileext = ".json")
  output_file <- tempfile(fileext = ".json")
  jsonlite::write_json(inp, input_file, auto_unbox = TRUE)
  
  # Spawn separate Rscript process — model loads and dies with the process
  status <- system2(
    "Rscript",
    args = c("/app/predict_once_meta.R", cancer_type, input_file, output_file),
    wait    = TRUE,
    timeout = 300
  )
  
  if (status != 0 || !file.exists(output_file)) {
    unlink(c(input_file, output_file))
    res$status <- 500
    return(list(error = "Prediction process failed. Check cancer type and input values."))
  }
  
  # Read result and clean up
  result <- jsonlite::fromJSON(output_file)
  unlink(c(input_file, output_file))
  return(result)
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
  
  example <- switch(cancer_type,
                    breast = list(age=55, tnm_n="N1", tnm_t="T2", er_status="Yes", pr_status="Yes", her2_status="No", grade="Grade 2"),
                    kidney = list(age=56, sex="male", tnm_n="N0", tnm_t="T1", surgical_factors="Yes", fuhrman_grade="2"),
                    LNSC   = list(age=65, sex="male", tnm_n="N1", tnm_t="T2"),
                    list(age=60, sex="male", tnm_n="N1", tnm_t="T2")
  )
  
  return(list(example_request = example))
}