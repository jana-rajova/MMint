# Install and load required packages
install_and_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

# List of packages to install and load
packages <- c(
  "Matrix", "Seurat", "ggplot2", "gridExtra", "RColorBrewer", 
  "ggplotify", "reticulate", "remotes", "magrittr", "dplyr", 
  "spdep", "kableExtra", "harmony", "tibble", "tidyr", 
  "factoextra", "viridis", "stringr", "pheatmap", "stats", 
  "patchwork", "ape", "fossil", "msigdbr"
)

# Install and load each package
sapply(packages, install_and_load)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GSVA", "decoupleR", "infercnv", "OmnipathR"), update = FALSE, ask = FALSE)

# Install packages from GitHub
remotes::install_github("saezlab/OmnipathR")
remotes::install_github("saezlab/decoupleR")
remotes::install_github("carmonalab/STACAS")
remotes::install_github("data2intelligence/SpaCET")
remotes::install_github("immunogenomics/presto")
remotes::install_github("theMILOlab/SPATA2")

# Install anndata (Python package)
reticulate::install_miniconda(update = TRUE, force = TRUE)
reticulate::conda_install('r-reticulate', 'anndata')

# Install leidenalg (Python package)
reticulate::conda_install('r-reticulate', 'leidenalg')

# Set Seurat object assay version
options(Seurat.object.assay.version = "v5")
print("All packages installed and loaded successfully!")
