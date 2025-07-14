# trial of functions

# load dependencies
library(tidyverse)
library(future)
library(future.apply)

# load functions for package
source("./R/makecwm.R") # taken from the ecole package (deprecated on CRAN)
source("./R/specifiR.R")


# LOAD DATA ####
otu <- readRDS("./sandbox/soils_otu_low_24.rds")
groups <- readRDS("./sandbox/soils_env_low_24.rds") %>% pluck("species")


# TRY FUNCTION ####
out <- specifiR(comm = otu,
                groups = groups,
                seed = 123,
                n.perm = 99,
                pval.cutoff = 0.05,
                max.ratio = 0)

out$community_specificity_index
out$taxon_specificity_index
out$isa_results
