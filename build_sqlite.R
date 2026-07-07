# build_sqlite.R
# Ingest all lookup tables (Aim 1 + Aim 2, gzipped CSVs) into a single SQLite
# database the API queries per-request. This replaces loading every table into
# RAM at startup — runtime RAM stays small and constant regardless of how many
# cancers/rows exist.
#
# Each table is a WITHOUT ROWID table whose PRIMARY KEY is the full set of input
# columns. The rows are clustered on that key, so a full-key lookup is a single
# B-tree seek and there is no separate (duplicate) index — keeping the file lean.
#
# All columns are stored as TEXT and matched as strings. The numeric bucket
# columns (age, psa, cores, cea) already have canonical short text forms in the
# CSVs (e.g. "7", "0.1", "3.75"), which match as.character() on the API side.
#
# Usage: Rscript build_sqlite.R [out_db]   (defaults to lookups.db)

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
})

args   <- commandArgs(trailingOnly = TRUE)
out_db <- if (length(args) >= 1) args[1] else "lookups.db"

aim1_dir <- "lookup_tables"
aim2_dir <- "lookup_tables_aim2"

# Output (non-input) columns — everything else is part of the lookup key.
AIM1_OUT <- c("no_metastasis", "metastasis", "risk_level")
AIM2_OUT <- c("p_mets")

if (file.exists(out_db)) file.remove(out_db)
con <- dbConnect(RSQLite::SQLite(), out_db)
on.exit(dbDisconnect(con))
dbExecute(con, "PRAGMA journal_mode = OFF")
dbExecute(con, "PRAGMA synchronous = OFF")

ingest <- function(gz_file, table, out_cols) {
  df <- read.csv(gzfile(gz_file), colClasses = "character",
                 stringsAsFactors = FALSE, check.names = FALSE)
  in_cols  <- setdiff(names(df), out_cols)
  col_defs <- paste(sprintf('"%s" TEXT', names(df)), collapse = ", ")
  pk       <- paste(sprintf('"%s"', in_cols), collapse = ", ")
  dbExecute(con, sprintf('DROP TABLE IF EXISTS "%s"', table))
  dbExecute(con, sprintf('CREATE TABLE "%s" (%s, PRIMARY KEY (%s)) WITHOUT ROWID',
                         table, col_defs, pk))
  dbAppendTable(con, table, df)
  nrow(df)
}

# ── Aim 1 ────────────────────────────────────────────────────────────────
aim1_files <- list.files(aim1_dir, pattern = "_lookup\\.csv\\.gz$", full.names = TRUE)
cat(sprintf("Aim 1: %d tables\n", length(aim1_files)))
for (f in aim1_files) {
  stem  <- sub("_lookup\\.csv\\.gz$", "", basename(f))
  table <- paste0("aim1_", tolower(stem))
  n <- ingest(f, table, AIM1_OUT)
  cat(sprintf("  %s (%d rows)\n", table, n))
}

# ── Aim 2 ────────────────────────────────────────────────────────────────
aim2_files <- list.files(aim2_dir, pattern = "_lookup\\.csv\\.gz$", full.names = TRUE)
cat(sprintf("Aim 2: %d tables\n", length(aim2_files)))
for (f in aim2_files) {
  stem  <- sub("_lookup\\.csv\\.gz$", "", basename(f))
  table <- paste0("aim2_", tolower(stem))
  n <- ingest(f, table, AIM2_OUT)
  cat(sprintf("  %s (%d rows)\n", table, n))
}

dbExecute(con, "PRAGMA journal_mode = DELETE")
dbExecute(con, "VACUUM")
dbExecute(con, "ANALYZE")

sz_mb <- file.info(out_db)$size / 1024 / 1024
cat(sprintf("Wrote %s (%.1f MB, %d tables)\n",
            out_db, sz_mb, length(aim1_files) + length(aim2_files)))
