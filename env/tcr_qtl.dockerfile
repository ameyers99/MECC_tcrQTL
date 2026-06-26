# ---- Base image ----
FROM rocker/tidyverse:4.2.2

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# ---- System dependencies ----
RUN apt-get update && \
apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    curl \
    git \
    unzip \
    cmake \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    libbz2-dev \
    liblapack-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    pkg-config \
    libblas-dev \
    libglpk-dev \
    libgit2-dev \
    glpk-utils \
    gfortran \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ---- Install PLINK 1.9 ----
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip && \
    unzip plink_linux_x86_64_20250819.zip && \
    mv plink /usr/local/bin && chmod +x /usr/local/bin/plink && \
    rm plink_linux_x86_64_20250819.zip
   
# ---- Install GCTA 1.94.1 ----
RUN wget https://github.com/jianyangqt/gcta/releases/download/v1.94.1/gcta-1.94.1-linux-x86_64-static && \
    mv gcta-1.94.1-linux-x86_64-static /usr/local/bin/gcta64 && chmod +x /usr/local/bin/gcta64

# ---- Install GEMMA 0.98.5 ----
RUN wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz && \
    gunzip gemma-0.98.5-linux-static-AMD64.gz && \
    cp gemma-0.98.5-linux-static-AMD64 /usr/local/bin/gemma && \
    chmod +x /usr/local/bin/gemma && \
    rm gemma-0.98.5-linux-static-AMD64

# ---- Install R packages ----
RUN R -e "options(error = function() quit(status = 1)); \
          install.packages('remotes', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('remotes', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('R.utils', '2.12.2', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('R.utils', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('optparse', '1.7.5', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('optparse', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('rbibutils', '2.4.1', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('rbibutils', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('numDeriv', '2016.8-1.1', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('numDeriv', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('ucminf', '1.2.3', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('ucminf', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('RcppEigen', '0.3.4.0.2', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('RcppEigen', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('reformulas', '0.4.4', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('reformulas', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('nloptr', '2.2.1', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('nloptr', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('minqa', '1.2.8', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('minqa', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('Rdpack', '2.6.6', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('Rdpack', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('ordinal', '2025.12-29', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('ordinal', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('lme4', '2.0-1', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('lme4', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('jomo', '2.7-6', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('jomo', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('pan', '1.9', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('pan', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('iterators', '1.0.14', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('iterators', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('shape', '1.4.6.1', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('shape', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('foreach', '1.5.2', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('foreach', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('operator.tools', '1.6.3.1', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('operator.tools', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('mitml', '0.4-5', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('mitml', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('glmnet', '5.0', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('glmnet', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('formula.tools', '1.7.1', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('formula.tools', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('mice', '3.19.0', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('mice', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('logistf', '1.24.1', repos='https://cloud.r-project.org', dependencies = TRUE, upgrade = FALSE); \
          stopifnot(requireNamespace('logistf', quietly=TRUE))"

RUN R -e "options(error = function() quit(status = 1)); \
          remotes::install_version('data.table', '1.14.8', repos='https://cloud.r-project.org'); \
          stopifnot(requireNamespace('data.table', quietly=TRUE))"

# github packages
RUN R -e "remotes::install_github('bhattacharya-a-bt/isotwas@a1ddc1c873803c8866f111c03c24b052e7c40d71'); \
          remotes::install_github('gabraham/plink2R/plink2R@d74be015e8f54d662b96c6c2a52a614746f9030d')"

# ---- Default command ----
CMD ["/bin/bash"]
