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
    liblzma-dev \
    libcurl4-openssl-dev \
    pkg-config \
    libgit2-dev \
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
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org'); \
          remotes::install_version('R.utils', '2.12.2', repos='https://cloud.r-project.org'); \
          remotes::install_version('optparse', '1.7.5', repos='https://cloud.r-project.org'); \
          remotes::install_version('glmnet', '4.1-6', repos='https://cloud.r-project.org'); \
          remotes::install_version('data.table', '1.14.8', repos='https://cloud.r-project.org'); \
          remotes::install_version('logistf', '1.24.1', repos='https://cloud.r-project.org')"

# github packages
RUN R -e "remotes::install_github('bhattacharya-a-bt/isotwas@a1ddc1c873803c8866f111c03c24b052e7c40d71'); \
          remotes::install_github('gabraham/plink2R/plink2R@d74be015e8f54d662b96c6c2a52a614746f9030d')"

# ---- Default command ----
CMD ["/bin/bash"]