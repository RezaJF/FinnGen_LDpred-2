# ldpred2 / bigsnpr PGS scoring image.
#
# Pure SCORING engine: applies precomputed per-SNP weights to imputed bgen
# dosages via bigsnpr (snp_readBGEN + snp_match + big_prodVec). No LDpred2
# inference, no LD reference is installed or used.
#
# Build + push (see docker/cloudbuild.yaml for the Cloud Build recipe):
#   docker build -f docker/Dockerfile \
#     -t REGION-docker.pkg.dev/PROJECT/REGISTRY/ldpred2:0.1 .
#   docker push REGION-docker.pkg.dev/PROJECT/REGISTRY/ldpred2:0.1

# rocker/r-ver pins R itself; bigsnpr pulls bigstatsr as a dependency.
FROM rocker/r-ver:4.3.2

LABEL maintainer="FinnGen long_gwas" \
      description="bigsnpr score-only application of precomputed LDAK MegaPRS PGS weights to FinnGen R13 bgen" \
      version="0.1"

ENV DEBIAN_FRONTEND=noninteractive

# System libraries: bigsnpr/bigstatsr build (C++), data.table (OpenMP),
# RSQLite (reads the .bgi index), zlib for bgen, plus curl/ssl for R pkg fetch.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        libz-dev \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libsqlite3-dev \
        libomp-dev \
        ca-certificates \
        procps \
    && rm -rf /var/lib/apt/lists/*

# R packages. Versions pinned where reasonable via a fixed CRAN snapshot
# (Posit Package Manager dated snapshot => reproducible builds).
# NOTE: bigsnpr::snp_readBGEN reads the .bgi SQLite index through dplyr's
# SQLite backend, which requires dbplyr (both are bigsnpr *Suggests*, so
# installing bigsnpr alone does NOT pull them). Without them snp_readBGEN
# aborts with "Please install package 'dbplyr'". Install explicitly.
RUN R -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2024-01-15')); \
          install.packages(c('optparse','data.table','bigreadr','RSQLite','Rcpp','RcppEigen','dplyr','dbplyr','Matrix'), Ncpus = 4); \
          install.packages('bigsnpr', Ncpus = 4); \
          stopifnot(requireNamespace('bigsnpr', quietly = TRUE), \
                    requireNamespace('bigstatsr', quietly = TRUE), \
                    requireNamespace('bigsparser', quietly = TRUE), \
                    requireNamespace('data.table', quietly = TRUE), \
                    requireNamespace('RSQLite', quietly = TRUE), \
                    requireNamespace('dplyr', quietly = TRUE), \
                    requireNamespace('dbplyr', quietly = TRUE), \
                    requireNamespace('Matrix', quietly = TRUE), \
                    requireNamespace('optparse', quietly = TRUE))"

# Bake all pipeline scripts into the image (scoring + inference + LD reference).
COPY scripts/ /scripts/
COPY tests/ /tests/
RUN chmod a+x /scripts/*.R

WORKDIR /scripts
