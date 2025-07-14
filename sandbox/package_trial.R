# trial of functions

# load dependencies
library(tidyverse)


if (!requireNamespace("specifiR", quietly = TRUE)){
  devtools::install_github("gzahn/specifiR")
}

library(specifiR)

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
