# alorenzetti 202106

# description ####
# this script is the scaffold
# running all the scripts
# required to perform
# the differential expression analysis

# loading required variables ####
# setting up number of threads
threads=6

# should we run quality control?
reportFlag="no"

# DE analysis thresholds
# DESeq2 adjusted pval
padjthreshold = 0.01

# DESeq2 log2FoldChange
log2fcthreshold = 1

# hypergeometric p threshold
qthr = 0.01

# creating data directory
if(!dir.exists("data")){dir.create("data")}

# creating results directory
if(!dir.exists("results")){dir.create("results")}

# creating plots directory
if(!dir.exists("plots")){dir.create("plots")}

# sourcing ####
# loading libs
source("scripts/01_loadingLibs.R")
if(reportFlag == "yes"){source("scripts/02_qualityControl.R")}
source("scripts/03_functionalCategorization.R")
source("scripts/04_deAnalysis.R")
source("scripts/05_results.R")
source("scripts/06_enrichmentAnalysisAndFigs.R")