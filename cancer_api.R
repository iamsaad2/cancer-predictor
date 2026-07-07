# cancer_api.R — SQLite-lookup backed API
# Predictions are precomputed lookup grids stored in an indexed SQLite database
# (lookups.db, built by build_sqlite.R). Each request is a single-row keyed
# lookup, so the process holds almost no data in RAM regardless of table count.
library(plumber)
library(jsonlite)
library(DBI)
library(RSQLite)

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
  "melanoma","liver","kidney","ovary","retroper","testis","LNSC","LSC",
  "uterine","thyroid","pancreas","stomach","cervix"
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
  ),
  uterine = list(
    age=list(col="AGE",type="age"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str")
  ),
  thyroid = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str"),
    focality=list(col="thyroid_cat",type="exact_str")   # Aim 1 only; dropped in Aim 2
  ),
  pancreas = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str")
  ),
  stomach = list(
    age=list(col="AGE",type="age"),
    sex=list(col="SEX",type="sex"),
    tnm_n=list(col="TNM_N_cat",type="exact_str"),
    tnm_t=list(col="TNM_T_cat",type="exact_str")
  ),
  cervix = list(
    age=list(col="AGE",type="age"),
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

# ── SQLite lookup store ──────────────────────────────────────────────────
# One read-only connection; every prediction is a single keyed row lookup.
# Nothing is held in RAM beyond SQLite's small page cache.
db_path <- Sys.getenv("LOOKUP_DB", unset = "/app/lookups.db")
cat(sprintf("Opening lookup DB %s\n", db_path))
DB        <- dbConnect(RSQLite::SQLite(), dbname = db_path, flags = RSQLite::SQLITE_RO)
DB_TABLES <- dbListTables(DB)
cat(sprintf("  %d lookup tables available\n", length(DB_TABLES)))

# Fetch the row(s) matching the given column = value pairs (all TEXT).
db_lookup <- function(tbl, cols, vals, select = "*") {
  where <- paste(sprintf('"%s" = ?', cols), collapse = " AND ")
  sql   <- sprintf('SELECT %s FROM "%s" WHERE %s', select, tbl, where)
  dbGetQuery(DB, sql, params = as.list(as.character(vals)))
}

# On a failed lookup, return the allowed values for a column if `val` isn't one
# of them (used to build a helpful error), else NULL.
db_offending <- function(tbl, col, val) {
  allowed <- dbGetQuery(DB, sprintf('SELECT DISTINCT "%s" AS v FROM "%s"', col, tbl))$v
  if (as.character(val) %in% allowed) NULL else sort(allowed)
}

lookup_predict <- function(cancer, inp) {
  tbl <- paste0("aim1_", tolower(cancer))
  if (!tbl %in% DB_TABLES) stop(sprintf("No lookup table for cancer: %s", cancer))
  map <- INPUT_MAPS[[cancer]]

  cols <- character(0); vals <- character(0); fields <- character(0)
  for (field in names(map)) {
    raw <- inp[[field]]
    if (is.null(raw) || (is.character(raw) && nchar(raw) == 0)) {
      stop(sprintf("Missing required field: %s", field))
    }
    m      <- map[[field]]
    cols   <- c(cols, m$col)
    vals   <- c(vals, as.character(transform_value(raw, m$type)))
    fields <- c(fields, field)
  }

  row <- db_lookup(tbl, cols, vals, select = "no_metastasis, metastasis, risk_level")
  if (nrow(row) == 0) {
    for (i in seq_along(cols)) {
      allowed <- db_offending(tbl, cols[i], vals[i])
      if (!is.null(allowed)) {
        stop(sprintf("Invalid %s='%s' (mapped to '%s'). Expected one of: %s",
                     fields[i], inp[[fields[i]]], vals[i], paste(allowed, collapse = ", ")))
      }
    }
    stop("No matching row for the given inputs")
  }
  if (nrow(row) != 1) stop(sprintf("Internal: expected 1 matching row, got %d", nrow(row)))
  list(
    no_metastasis = as.numeric(row$no_metastasis[1]),
    metastasis    = as.numeric(row$metastasis[1]),
    risk_level    = row$risk_level[1]
  )
}

# ── Aim 2 — site-specific metastasis prediction ──────────────────────────
# Reuses Aim 1's per-cancer covariate maps but bucketing AGE, and adds mets_* fields per site.

AIM2_INPUT_MAPS <- lapply(INPUT_MAPS, function(m) {
  if (!is.null(m$age)) m$age$type <- "age_bucket"
  m
})
# Aim 2 thyroid model does not use tumor focality (that covariate is Aim 1 only).
AIM2_INPUT_MAPS$thyroid$focality <- NULL

# Predict probability of mets to a single target site (queries the SQLite store).
predict_sites_one <- function(cancer, site, inp) {
  tbl <- paste0("aim2_", tolower(cancer), "_", site)
  if (!tbl %in% DB_TABLES) stop(sprintf("No Aim 2 lookup table for %s/%s", cancer, site))
  map <- AIM2_INPUT_MAPS[[cancer]]

  cols <- character(0); vals <- character(0)
  # Aim 1 covariates (age bucketed)
  for (field in names(map)) {
    raw <- inp[[field]]
    if (is.null(raw) || (is.character(raw) && nchar(raw) == 0)) {
      stop(sprintf("Missing required field: %s", field))
    }
    m    <- map[[field]]
    cols <- c(cols, m$col)
    vals <- c(vals, as.character(transform_value(raw, m$type)))
  }

  # mets_<other 3 sites>
  for (s in setdiff(AIM2_SITES, site)) {
    col <- paste0("mets_", s)
    raw <- inp[[col]]
    if (is.null(raw)) stop(sprintf("Missing required field: %s", col))
    val <- as.character(raw)
    if (!val %in% c("Yes","No")) stop(sprintf("%s must be 'Yes' or 'No', got '%s'", col, raw))
    cols <- c(cols, col); vals <- c(vals, val)
  }

  # mets_other (DB column is "mate_other" — typo carried from model)
  raw <- inp$mets_other
  if (is.null(raw)) stop("Missing required field: mets_other")
  val <- as.character(raw)
  if (!val %in% c("Yes","No")) stop(sprintf("mets_other must be 'Yes' or 'No', got '%s'", raw))
  cols <- c(cols, "mate_other"); vals <- c(vals, val)

  row <- db_lookup(tbl, cols, vals, select = "p_mets")
  if (nrow(row) != 1) {
    stop(sprintf("Internal: expected 1 row for %s/%s, got %d", cancer, site, nrow(row)))
  }
  as.numeric(row$p_mets[1])
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
  loaded_aim1_tables = sum(startsWith(DB_TABLES, "aim1_")),
  loaded_aim2_tables = sum(startsWith(DB_TABLES, "aim2_"))
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
    thyroid = list(sex="Sex (male, female)",
                   focality="Tumor focality (Solitary, Multifocal, Other/Unknown)"),
    uterine = , cervix = list(),
    pancreas = , stomach = , LNSC = , LSC = list(sex="Sex (male, female)"),
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
    thyroid  = list(age=50, sex="female", tnm_n="N0", tnm_t="T2", focality="Solitary"),
    uterine  = , cervix = list(age=55, tnm_n="N1", tnm_t="T2"),
    LNSC     = list(age=65, sex="male", tnm_n="N1", tnm_t="T2"),
    list(age=60, sex="male", tnm_n="N1", tnm_t="T2")
  )
  list(example_request = example)
}
