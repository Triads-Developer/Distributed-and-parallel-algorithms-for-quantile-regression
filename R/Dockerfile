FROM ubuntu:20.04

# install packages using apt-get

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update --fix-missing && \
apt-get install -y build-essential \
libcurl4-gnutls-dev libxml2-dev \
libssl-dev software-properties-common r-cran-rcpparmadillo && \
apt-get clean

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/'
RUN apt update
RUN apt install -y r-base


# extend image to include tidyverse and devtools R packages

RUN R -e "install.packages(c('devtools','tidyverse', 'doMC','doParallel', 'bigmemory', 'ff'), dependencies=TRUE, repos='https://cran.rstudio.com')"
