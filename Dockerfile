FROM rocker/tidyverse:4.4.0
LABEL maintainer="jrokita@childrensnational.org"
WORKDIR /rocker-build/

ARG github_pat=$GITHUB_PAT

ENV GITHUB_PAT=$github_pat

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

### Install apt-getable packages to start
#########################################

# Installing all apt required packages at once
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    jq \
    libgmp3-dev \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev \
    libpoppler-cpp-dev \
    libglpk-dev \
    libncurses5 \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libgdbm-dev \
    libdb5.3-dev \
    libbz2-dev \
    libexpat1-dev \
    liblzma-dev \
    libffi-dev \
    libuuid1 \
    wget \
    xorg \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and install Python 3.11
RUN cd /usr/src && \
    wget https://www.python.org/ftp/python/3.11.0/Python-3.11.0.tgz && \
    tar xzf Python-3.11.0.tgz && \
    cd Python-3.11.0 && \
    ./configure --enable-optimizations && \
    make altinstall && \
    rm -rf /usr/src/Python-3.11.0.tgz

# Setup the default python commands to use Python 3.11
RUN ln -s /usr/local/bin/python3.11 /usr/local/bin/python3 && \
    ln -s /usr/local/bin/python3.11 /usr/local/bin/python
RUN python3 -m pip install --upgrade pip

# Set working directory
WORKDIR /home/rstudio


# Install python packages
##########################

# Install python3 tools and ALL dependencies
RUN pip3 install \
    "appdirs==1.4.4" \
    "attrs==23.1.0" \
    "certifi==2023.5.7" \
    "chardet==5.1.0" \
    "ConfigArgParse==1.7" \
    "CrossMap==0.6.5" \
    "Cython==0.29.15" \
    "defusedxml==0.7.1" \
    "docutils==0.20" \
    "gitdb==4.0.10" \
    "GitPython==3.1.31" \
    "idna==3.4" \
    "importlib-metadata==6.6.0" \
    "ipykernel==6.23.0" \
    "ipython==8.13.2" \
    "ipython-genutils==0.2.0" \
    "jsonschema==4.17.3" \
    "jupyter-client==8.2.0" \
    "jupyter-core==5.3.0" \
    "MarkupSafe==2.1.2" \
    "matplotlib==3.7.1" \
    "nbconvert==7.4.0" \
    "nbformat==5.8.0" \
    "notebook==6.5.4" \
    "numpy==1.24.3" \
    "openpyxl==3.1.2" \
    "packaging==23.1" \
    "palettable==3.3.3" \
    "pandas==2.0.1" \
    "pandocfilters==1.5.0" \
    "parso==0.8.3" \
    "patsy==0.5.3" \
    "pexpect==4.8.0" \
    "pickleshare==0.7.5" \
    "plotnine==0.12.1" \
    "prompt-toolkit==3.0.38" \
    "psutil==5.9.5" \
    "ptyprocess==0.7.0" \
    "PuLP==2.8.0" \
    "Pygments==2.15.1" \
    "pyparsing==3.0.9" \
    "python-dateutil==2.8.2" \
    "pytz==2023.3" \
    "PyYAML==6.0" \
    "pyzmq==25.0.2" \
    "ratelimiter==1.2.0.post0" \
    "requests==2.30.0" \
    "rpy2==3.5.0" \
    "scikit-learn==1.2.2" \
    "scipy==1.10.1" \
    "seaborn==0.12.2" \
    "setuptools==46.3.0" \
    "six==1.16.0" \
    "snakemake==8.11.6" \
    "statsmodels==0.14.0" \
    "tornado==6.3.1" \
    "traitlets==5.9.0" \
    "urllib3==2.0.2" \
    "utils==1.0.1" \
    "webencodings==0.5.1" \
    "widgetsnbextension==4.0.7" \
    "wheel==0.34.2" \
    "wrapt==1.15.0" \
    "zipp==3.15.0" \
    && rm -rf /root/.cache/pip/wheels

# Install java
RUN apt-get update && apt-get install -y openjdk-11-jdk

# Required for running matplotlib in Python in an interactive session
RUN apt-get -y --no-install-recommends install \
    python3-tk

# Standalone tools and libraries
################################

# Required for mapping segments to genes
# Add bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && rm -f bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin && \
    cd .. && rm -rf bedtools2

# Add bedops per the BEDOPS documentation
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.37/bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    tar -jxvf bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    rm -f bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    mv bin/* /usr/local/bin

# HTSlib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -jxvf htslib-1.9.tar.bz2 && rm -f htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf htslib-1.9

# GenomeTools
RUN wget http://genometools.org/pub/genometools-1.6.2.tar.gz  && \
    tar -zxvf genometools-1.6.2.tar.gz  && rm -f genometools-1.6.2.tar.gz && \
    cd genometools-1.6.2 && \
    make cairo=no && \
    make prefix=/usr/local cairo=no install && \
    cd .. && rm -rf genometools-1.6.2

#### R packages
###############

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.19', ask = FALSE)"

# Install R packages
RUN R -e 'BiocManager::install(c( \
  "annotatr", \
  "AnnotationDbi", \
  "arrow", \
  "bedr", \
  "BSgenome.Hsapiens.UCSC.hg19", \
  "BSgenome.Hsapiens.UCSC.hg38", \
  "caret", \
  "class", \
  "cluster", \
  "ComplexHeatmap", \
  "corrplot", \
  "d3r", \
  "data.table", \
  "DESeq2", \
  "dplyr", \
  "DT", \
  "e1071", \ 
  "EDASeq",  \
  "edgeR", \
  "EnsDb.Hsapiens.v86", \
  "ensembldb", \
  "flextable", \
  "foreign", \
  "GenomicRanges", \
  "GenVisR", \
  "GGally", \
  "ggbio", \
  "ggfortify", \
  "ggplot2", \
  "ggpubr", \
  "ggsignif", \
  "glmnet", \
  "glmnetUtils", \
  "gplots", \
  "gridGraphics", \
  "GSVA", \
  "ids", \
  "irlba", \
  "lattice", \
  "maftools", \
  "MASS", \
  "Matrix", \
  "msigdbr", \
  "multipanelfigure", \
  "mygene", \
  "openxlsx", \
  "optparse", \
  "org.Hs.eg.db", \
  "pheatmap", \
  "preprocessCore", \
  "qdapRegex", \
  "R.utils", \
  "RColorBrewer", \
  "rJava", \
  "rlist", \
  "rpart", \
  "rprojroot", \
  "rtracklayer", \
  "Rtsne", \
  "RUVSeq", \
  "spatial", \
  "survival", \
  "survminer", \
  "survMisc", \
  "sva", \
  "tidyr", \
  "TxDb.Hsapiens.UCSC.hg38.knownGene", \
  "umap" , \
  "UpSetR", \
  "uwot", \
  "VennDiagram", \
  "viridis", \
  "vroom" \
  ))'


# package required for immune deconvolution
RUN R -e "remotes::install_github('omnideconv/immunedeconv', ref = 'a7e4ee9993aa94f268e862263eaf226a251514f9', dependencies = TRUE)"

RUN R -e "remotes::install_github('const-ae/ggupset', ref = '7a33263cc5fafdd72a5bfcbebe5185fafe050c73', dependencies = TRUE)"

# Need this package to make plots colorblind friendly
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"

# package required for shatterseek
RUN R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), remotes::install_github('parklab/ShatterSeek', ref = '83ab3effaf9589cc391ecc2ac45a6eaf578b5046', dependencies = TRUE))"

# Need this specific version of circlize so it has hg38
RUN R -e "remotes::install_github('jokergoo/circlize', ref = 'b7d86409d7f893e881980b705ba1dbc758df847d', dependencies = TRUE)"

# signature.tools.lib needed for mutational-signatures 
RUN R -e "remotes::install_github('Nik-Zainal-Group/signature.tools.lib', ref = '59a3a3236f16f0c1383d0ab125fec8a251d7f42d', dependencies = TRUE)"

# Molecular subtyping MB
RUN R -e "remotes::install_github('d3b-center/medullo-classifier-package', ref = 'e3d12f64e2e4e00f5ea884f3353eb8c4b612abe8', dependencies = TRUE, upgrade = FALSE)"  
    
# More recent version of sva required for molecular subtyping MB
RUN R -e "remotes::install_github('jtleek/sva-devel@123be9b2b9fd7c7cd495fab7d7d901767964ce9e', dependencies = FALSE, upgrade = FALSE)"

# Packages required for de novo mutational signatures
#RUN install2.r --error --deps TRUE \
#    lsa

# Packages for fusion annotation using annoFuse package
RUN R -e "remotes::install_github('d3b-center/annoFuseData',ref = '321bc4f6db6e9a21358f0d09297142f6029ac7aa', dependencies = TRUE)"
RUN R -e "remotes::install_github('d3b-center/annoFuse',ref = '55b4b862429fe886790d087b2f1c654689c691c4', dependencies = TRUE)"

# Latest deconstructSigs release for mut sigs analyses
RUN R -e "remotes::install_github('raerose01/deconstructSigs', ref = '41a705c5d80848121347d448cf9e2c58ff6b81ac', dependencies = TRUE)"

# MATLAB Compiler Runtime is required for GISTIC, MutSigCV
# Install steps are adapted from usuresearch/matlab-runtime
# https://hub.docker.com/r/usuresearch/matlab-runtime/dockerfile

# This is the version of MCR required to run the precompiled version of GISTIC
RUN mkdir /mcr-install-v83 && \
    mkdir /opt/mcr && \
    cd /mcr-install-v83 && \
    wget https://www.mathworks.com/supportfiles/downloads/R2014a/deployment_files/R2014a/installers/glnxa64/MCR_R2014a_glnxa64_installer.zip && \
    unzip -q MCR_R2014a_glnxa64_installer.zip && \
    rm -f MCR_R2014a_glnxa64_installer.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install-v83

WORKDIR /home/rstudio/
# GISTIC installation
RUN mkdir -p gistic_install && \
    cd gistic_install && \
    wget -q https://anaconda.org/HCC/gistic2/2.0.23/download/linux-64/gistic2-2.0.23-0.tar.bz2 && \
    tar jxf gistic2-2.0.23-0.tar.bz2 && \
    rm -f gistic2-2.0.23-0.tar.bz2 && \
    rm -rf MCR_Installer && \
    chown -R rstudio:rstudio /home/rstudio/gistic_install && \
    chmod 755 /home/rstudio/gistic_install

WORKDIR /home/rstudio/
# AWS CLI installation
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    sudo ./aws/install && \
    rm -rf aws*

# Install Desal latest release (v2.1.1)- converter for JSON, TOML, YAML, XML and CSV data formats
RUN sudo wget -qO /usr/local/bin/dasel "https://github.com/TomWright/dasel/releases/download/v2.1.1/dasel_linux_amd64" && \
    sudo chmod a+x /usr/local/bin/dasel

# Reset the frontend variable for interactive
ENV DEBIAN_FRONTEND=

WORKDIR /rocker-build/
