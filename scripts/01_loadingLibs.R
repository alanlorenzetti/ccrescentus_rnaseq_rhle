# alorenzetti 202106

# description ####
# this script will load all
# the libs required to perform
# the differential expression analysis
# and load essential variables

# loading libs ####
# biocmanager is required for loading and installing packages
if(!require("BiocManager")){install.packages("BiocManager"); library("BiocManager")}

# pacman is a nice package manager; make it easier to load and install packages
if(!require("pacman")){install.packages("pacman"); library("pacman")}

# list of required packages
requiredPacks = c("tidyverse",
                  "ggpubr",
                  "ggthemes",
                  "pheatmap",
                  "RColorBrewer",
                  "genefilter",
                  "DT",
                  "ggplot2",
                  "plotly",
                  "Rqc",
                  "rtracklayer",
                  "QuasR",
                  "GenomicAlignments",
                  "Rsamtools",
                  "GenomicFeatures",
                  "BiocParallel",
                  "DESeq2",
                  "gage",
                  "magrittr",
                  "UniProt.ws",
                  "ggthemes",
                  "readtext")

# loading required cran packages
p_load(char=requiredPacks)

# learning tableau colors
tab10 = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`

# registering multiple processors
register(MulticoreParam(workers = threads))

# setting ggtheme
theme_set(theme_bw())
