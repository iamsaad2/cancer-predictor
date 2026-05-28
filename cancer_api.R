# cancer_api.R — CSV-lookup backed API
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

CANCER_TYPES <- c(
  "breast","prostate","colon","rectum","urine","esophagu",
  "melanoma","liver","kidney","ovary","retroper","testis","LNSC","LSC"
)

# Clinical bucket midpoints for numeric inputs (must match build_lookup.R)
PSA_BUCKETS   <- c(2, 7, 15, 35, 75, 150)
CORES_BUCKETS <- c(0.1, 0.3, 0.5, 0.7, 0.9)
CEA_BUCKETS   <- c(1, 3.75, 7.5, 15, 35, 100)

# Aim 2 — age is bucketed (must match build_lookup_aim2.R)
AGE_BUCKETS <- seq(5, 95, by = 10)
AIM2_SITES  <- c("bone", "brain", "liver", "lung")

# Per-cancer mapping: API field → (CSV column, value type)
# Types: "exact_str", "sex", "age", "psa", "cores", "cea"
INPUT_MAPS <- list(
  breast = list(
    age=list(col="AGE",type="age"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    er_status=list(col="breast_er_cat",type="exact_str"),
    pr_status=list(col="breast_pr_cat",type="exact_str"),
    her2_status=list(col="breast_her2_cat",type="exact_str"),
    grade=list(col="breast_br_cat",type="exact_str")
  ),
  prostate = list(
    age=list(col="AGE",type="age"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    psa=list(col="prostate_psa",type="psa"),
    core_ratio=list(col="prostate_cores",type="cores"),
    gleason=list(col="prostate_gs_cat",type="exact_str")
  ),
  colon = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    cea=list(col="colon_cea",type="cea")
  ),
  rectum = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    cea_status=list(col="rectum_cea_cat",type="exact_str")
  ),
  urine = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str")
  ),
  esophagu = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str")
  ),
  melanoma = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str")
  ),
  liver = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    afp_status=list(col="liver_afp_cat",type="exact_str"),
    fibrosis_score=list(col="liver_fs_cat",type="exact_str")
  ),
  kidney = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    surgical_factors=list(col="rcc_sf_cat",type="exact_str"),
    fuhrman_grade=list(col="rcc_fng_cat",type="exact_str")
  ),
  ovary = list(
    age=list(col="AGE",type="age"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    ca125_status=list(col="ovary_ca125_cat",type="exact_str")
  ),
  retroper = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    surgical_grade=list(col="retroper_sg_cat",type="exact_str")
  ),
  testis = list(
    age=list(col="AGE",type="age"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    afp_status=list(col="testis_afp_cat",type="exact_str"),
    hcg_status=list(col="testis_hcg_cat",type="exact_str"),
    ldh_status=list(col="testis_ldh_cat",type="exact_str")
  ),
  LNSC = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str")
  ),
  LSC = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str")
  )
)

nearest_bucket <- function(value, buckets) {
  buckets[which.min(abs(buckets - value))]
}

transform_value <- function(value, type) {
  if (type == "exact_str") return(as.character(value))
  if (type == "sex") {
    v <- as.character(value)
    if (v == "male")   return("1")
    if (v == "female") return("2")
    return(v)  # already "1"/"2"
  }
  v <- suppressWarnings(as.numeric(value))
  if (is.na(v)) stop(sprintf("Expected number, got: '%s'", value))
  if (type == "age")        return(min(max(as.integer(round(v)), 0), 99))
  if (type == "age_bucket") return(nearest_bucket(min(max(as.integer(round(v)), 0), 99), AGE_BUCKETS))
  if (type == "psa")        return(nearest_bucket(v, PSA_BUCKETS))
  if (type == "cores")      return(nearest_bucket(v, CORES_BUCKETS))
  if (type == "cea")        return(nearest_bucket(v, CEA_BUCKETS))
  stop(sprintf("Unknown input type: %s", type))
}

# Load CSVs once at startup. Keep numeric columns numeric so == compares cleanly.
LOOKUPS <- new.env()

prepare_df <- function(df, cancer) {
  map <- INPUT_MAPS[[cancer]]
  for (field in names(map)) {
    m <- map[[field]]
    if (m$type %in% c("age","psa","cores","cea")) {
      df[[m$col]] <- as.numeric(df[[m$col]])
    } else {
      # factor stores low-cardinality strings as 4-byte int codes — saves RAM
      df[[m$col]] <- as.factor(df[[m$col]])
    }
  }
  df$no_metastasis <- as.numeric(df$no_metastasis)
  df$metastasis    <- as.numeric(df$metastasis)
  df$risk_level    <- as.factor(df$risk_level)
  df
}

lookup_dir <- Sys.getenv("LOOKUP_DIR", unset = "/app/lookup_tables")
cat(sprintf("Loading lookup tables from %s\n", lookup_dir))
for (cancer in CANCER_TYPES) {
  f <- file.path(lookup_dir, paste0(cancer, "_lookup.csv"))
  if (!file.exists(f)) {
    cat(sprintf("  WARN: missing %s\n", f))
    next
  }
  df <- read.csv(f, colClasses = "character", stringsAsFactors = FALSE)
  LOOKUPS[[cancer]] <- prepare_df(df, cancer)
  cat(sprintf("  loaded %s (%d rows)\n", cancer, nrow(LOOKUPS[[cancer]])))
}

lookup_predict <- function(cancer, inp) {
  df  <- LOOKUPS[[cancer]]
  if (is.null(df)) stop(sprintf("No lookup table loaded for cancer: %s", cancer))
  map <- INPUT_MAPS[[cancer]]

  filtered <- df
  for (field in names(map)) {
    raw <- inp[[field]]
    if (is.null(raw) || (is.character(raw) && nchar(raw) == 0)) {
      stop(sprintf("Missing required field: %s", field))
    }
    m   <- map[[field]]
    val <- transform_value(raw, m$type)
    filtered <- filtered[filtered[[m$col]] == val, , drop = FALSE]
    if (nrow(filtered) == 0) {
      allowed <- sort(unique(df[[m$col]]))
      stop(sprintf("Invalid %s='%s' (mapped to '%s'). Expected one of: %s",
                   field, raw, val, paste(allowed, collapse = ", ")))
    }
  }

  if (nrow(filtered) != 1) {
    stop(sprintf("Internal: expected 1 matching row, got %d", nrow(filtered)))
  }
  filtered[1, ]
}

# ── Aim 2 — site-specific metastasis prediction ──────────────────────────
# Reuses Aim 1's per-cancer covariate maps but bucketing AGE, and adds mets_* fields per site.

AIM2_INPUT_MAPS <- lapply(INPUT_MAPS, function(m) {
  if (!is.null(m$age)) m$age$type <- "age_bucket"
  m
})

aim2_key <- function(cancer) tolower(cancer)  # CSVs use lowercase (lnsc, lsc)

AIM2_LOOKUPS <- new.env()

prepare_aim2_df <- function(df, cancer) {
  map <- AIM2_INPUT_MAPS[[cancer]]
  for (field in names(map)) {
    m <- map[[field]]
    if (m$type %in% c("age_bucket","psa","cores","cea")) {
      df[[m$col]] <- as.numeric(df[[m$col]])
    } else {
      df[[m$col]] <- as.factor(df[[m$col]])
    }
  }
  for (col in c("mets_bone","mets_brain","mets_liver","mets_lung","mate_other")) {
    if (col %in% names(df)) df[[col]] <- as.factor(df[[col]])
  }
  df$p_mets <- as.numeric(df$p_mets)
  df
}

aim2_lookup_dir <- Sys.getenv("AIM2_LOOKUP_DIR", unset = "/app/lookup_tables_aim2")
cat(sprintf("Loading Aim 2 lookup tables from %s\n", aim2_lookup_dir))
for (cancer in CANCER_TYPES) {
  for (site in AIM2_SITES) {
    key <- paste0(aim2_key(cancer), "_", site)
    f   <- file.path(aim2_lookup_dir, paste0(key, "_lookup.csv.gz"))
    if (!file.exists(f)) {
      cat(sprintf("  WARN: missing %s\n", f))
      next
    }
    df <- read.csv(gzfile(f), colClasses = "character", stringsAsFactors = FALSE)
    AIM2_LOOKUPS[[key]] <- prepare_aim2_df(df, cancer)
  }
}
cat(sprintf("  loaded %d Aim 2 tables\n", length(ls(AIM2_LOOKUPS))))

# Predict probability of mets to a single target site.
predict_sites_one <- function(cancer, site, inp) {
  key <- paste0(aim2_key(cancer), "_", site)
  df  <- AIM2_LOOKUPS[[key]]
  if (is.null(df)) stop(sprintf("No Aim 2 lookup table for %s/%s", cancer, site))
  map <- AIM2_INPUT_MAPS[[cancer]]

  filtered <- df
  # Aim 1 covariates (age bucketed)
  for (field in names(map)) {
    raw <- inp[[field]]
    if (is.null(raw) || (is.character(raw) && nchar(raw) == 0)) {
      stop(sprintf("Missing required field: %s", field))
    }
    m   <- map[[field]]
    val <- transform_value(raw, m$type)
    filtered <- filtered[filtered[[m$col]] == val, , drop = FALSE]
    if (nrow(filtered) == 0) {
      stop(sprintf("Invalid %s='%s' (mapped to '%s') for %s/%s", field, raw, val, cancer, site))
    }
  }

  # mets_<other 3 sites>
  for (s in setdiff(AIM2_SITES, site)) {
    col <- paste0("mets_", s)
    raw <- inp[[col]]
    if (is.null(raw)) stop(sprintf("Missing required field: %s", col))
    val <- as.character(raw)
    if (!val %in% c("Yes","No")) stop(sprintf("%s must be 'Yes' or 'No', got '%s'", col, raw))
    filtered <- filtered[filtered[[col]] == val, , drop = FALSE]
  }

  # mets_other (CSV column is "mate_other" — typo carried from model)
  raw <- inp$mets_other
  if (is.null(raw)) stop("Missing required field: mets_other")
  val <- as.character(raw)
  if (!val %in% c("Yes","No")) stop(sprintf("mets_other must be 'Yes' or 'No', got '%s'", raw))
  filtered <- filtered[filtered$mate_other == val, , drop = FALSE]

  if (nrow(filtered) != 1) {
    stop(sprintf("Internal: expected 1 row for %s/%s, got %d", cancer, site, nrow(filtered)))
  }
  filtered$p_mets[1]
}

#* @apiTitle Cancer Metastasis Prediction API (lookup-table backed)

#* @get /cancer-types
#* @serializer unboxedJSON
function() list(cancer_types = CANCER_TYPES)

#* @get /health
#* @serializer unboxedJSON
function() list(
  status = "healthy",
  timestamp = Sys.time(),
  available_cancer_types = length(CANCER_TYPES),
  loaded_aim1_tables = length(ls(LOOKUPS)),
  loaded_aim2_tables = length(ls(AIM2_LOOKUPS))
)

#* @get /inputs/<cancer_type>
#* @serializer unboxedJSON
function(cancer_type) {
  if (!cancer_type %in% CANCER_TYPES) return(list(error = "Invalid cancer type"))
  base_inputs <- list(
    age   = "Patient age (number)",
    tnm_n = "Lymph node status (N0, N1, N2, N3)",
    tnm_t = "Tumor status (T0, T1, T2, T3, T4)"
  )
  extra <- switch(cancer_type,
    breast = list(er_status="ER status (No, Yes)", pr_status="PR status (No, Yes)",
                  her2_status="HER2 status (No, Yes)",
                  grade="Grade (Grade 1, Grade 2, Grade 3, Grade 4)"),
    prostate = list(psa="PSA level (number)", core_ratio="Core ratio (number 0-1)",
                    gleason="Gleason score (Grade group 1..5)"),
    colon = list(sex="Sex (male, female)", cea="CEA level (number)"),
    rectum = list(sex="Sex (male, female)",
                  cea_status="CEA status (Negative/normal; within normal limits, Positive/elevated)"),
    urine = , esophagu = , melanoma = list(sex="Sex (male, female)"),
    liver = list(sex="Sex (male, female)",
                 afp_status="AFP status (Negative/normal, Positive/elevated)",
                 fibrosis_score="Fibrosis score (Fibrosis score 0-4, Fibrosis score 5-6)"),
    kidney = list(sex="Sex (male, female)",
                  surgical_factors="Surgical factors (Yes, No)",
                  fuhrman_grade="Fuhrman grade (1, 2, 3, 4)"),
    ovary = list(ca125_status="CA-125 status (Negative/normal; within normal limits, Positive/elevated)"),
    retroper = list(sex="Sex (male, female)", surgical_grade="Surgical grade (1, 2, 3)"),
    testis = list(afp_status="AFP status (Within normal limits, Elevated)",
                  hcg_status="HCG status (Within normal limits, Elevated)",
                  ldh_status="LDH status (Within normal limits, Elevated)"),
    LNSC = , LSC = list(sex="Sex (male, female)"),
    list(sex="Sex (male, female)")
  )
  c(base_inputs, extra)
}

#* @post /predict/<cancer_type>
#* @serializer unboxedJSON
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
  row <- tryCatch(lookup_predict(cancer_type, inp), error = function(e) e)
  if (inherits(row, "error")) {
    res$status <- 400
    return(list(error = conditionMessage(row)))
  }
  list(
    cancer_type = cancer_type,
    patient_data = inp,
    predictions = list(
      random_forest = list(
        no_metastasis = round(row$no_metastasis, 4),
        metastasis    = round(row$metastasis, 4),
        risk_level    = row$risk_level
      )
    ),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
}

#* @post /predict-sites/<cancer_type>
#* @serializer unboxedJSON
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

  predictions <- list()
  for (site in AIM2_SITES) {
    p <- tryCatch(predict_sites_one(cancer_type, site, inp), error = function(e) e)
    if (inherits(p, "error")) {
      res$status <- 400
      return(list(error = conditionMessage(p)))
    }
    predictions[[site]] <- list(
      p_mets     = round(p, 4),
      risk_level = ifelse(p > 0.5, "HIGH", "LOW")
    )
  }

  list(
    cancer_type  = cancer_type,
    patient_data = inp,
    predictions  = predictions,
    timestamp    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
}

#* @get /example/<cancer_type>
#* @serializer unboxedJSON
function(cancer_type) {
  example <- switch(cancer_type,
    breast   = list(age=55, tnm_n="N1", tnm_t="T2", er_status="Yes", pr_status="Yes", her2_status="No", grade="Grade 2"),
    prostate = list(age=65, tnm_n="N0", tnm_t="T2", psa=7.5, core_ratio=0.4, gleason="Grade group 2"),
    colon    = list(age=60, sex="male", tnm_n="N1", tnm_t="T3", cea=4.5),
    kidney   = list(age=56, sex="male", tnm_n="N0", tnm_t="T1", surgical_factors="Yes", fuhrman_grade="2"),
    LNSC     = list(age=65, sex="male", tnm_n="N1", tnm_t="T2"),
    list(age=60, sex="male", tnm_n="N1", tnm_t="T2")
  )
  list(example_request = example)
}
