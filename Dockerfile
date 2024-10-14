FROM rocker/rstudio:4.3.2

LABEL info='Melanoma brain metastases docker environment'
	
RUN mkdir -p home/rstudio/data home/rstudio/bin  
RUN mkdir -p home/rstudio/data/stdata home/rstudio/data/stutility_coordinates home/rstudio/data/scrnaseq/
RUN mkdir -p dstu/setup/R

RUN apt update && apt install -y  \
	cmake \
	libpq-dev \
	libudunits2-dev \
	libgdal-dev  \
	libgeos-dev  \
	libproj-dev   \
	libglpk-dev \
	cmake \
	build-essential  \
	libxml2-dev  \
	libssl-dev   \
	libmagick++-dev  \
	libstdc++6 \
	libgsl-dev  \
	jags \
	htop \
	pip \
	libcurl4-gnutls-dev \
	libxml2-dev \
	libgit2-dev \
	libharfbuzz-dev\
	libfribidi-dev \
	libhdf5-dev


RUN pip install Jupyter
RUN apt update
RUN apt install -y python3-venv python3-pip python3-virtualenv
RUN R -e 'install.packages("remotes", repo = "https://cloud.r-project.org")' \
	-e 'remotes::install_github("rstudio/reticulate")' \
  	-e 'reticulate::virtualenv_create()'

EXPOSE 1227
ADD setup/install_packages.R dstu/setup/R
RUN Rscript dstu/setup/R/install_packages.R

ADD data/stutility_alignment/ home/rstudio/data/stutility_alignment/
ADD data/st_data/  home/rstudio/data/stdata/
ADD data/scrnaseq.h5ad home/rstudio/data/scrnaseq/
COPY data/gencode_v19_gene_pos.txt home/rstudio/data/gencode_v19_gene_pos.txt
COPY bin/* home/rstudio/bin/
RUN chmod -R 777 /home/rstudio/data
RUN mkdir /home/rstudio/results/
RUN chown -R rstudio /home/rstudio/results/
