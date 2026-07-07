############################
# 1)  Builder  – clone repo
############################
FROM alpine/git AS lfsstage

ARG REPO_URL
ARG REF=main
ARG GH_TOKEN

# Lookup CSVs are not LFS-tracked, so we don't need git-lfs.
# Skip smudge so LFS pointer files for legacy model artifacts stay as stubs.
ENV GIT_LFS_SKIP_SMUDGE=1

RUN if [ -z "$GH_TOKEN" ]; then \
      git clone --depth 1 --branch "$REF" "$REPO_URL" /repo; \
    else \
      git clone --depth 1 --branch "$REF" "https://${GH_TOKEN}@${REPO_URL#https://}" /repo; \
    fi && \
    ls -lh /repo/lookup_tables

############################
# 2)  Runtime image
############################
FROM rocker/r-ver:4.3.3

RUN install2.r --error plumber jsonlite DBI RSQLite

WORKDIR /app

COPY cancer_api.R build_sqlite.R /app/
COPY --from=lfsstage /repo/lookup_tables /app/lookup_tables
COPY --from=lfsstage /repo/lookup_tables_aim2 /app/lookup_tables_aim2

# Bake the compressed lookup CSVs into a single indexed SQLite DB, then drop the
# CSVs. At runtime the API queries the DB per-request, so RAM stays small and
# constant regardless of how many lookup tables exist.
RUN Rscript build_sqlite.R /app/lookups.db && \
    rm -rf /app/lookup_tables /app/lookup_tables_aim2

ENV LOOKUP_DB=/app/lookups.db

EXPOSE 8000
CMD ["R", "-e", "pr <- plumber::plumb('cancer_api.R'); pr$run(host='0.0.0.0', port=8000)"]
