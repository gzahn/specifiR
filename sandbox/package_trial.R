# trial of functions

# load dependencies (for this script)
library(tidyverse)
library(phyloseq)
source("./R/specifiR_physeq.R")
# install specifiR if not already installed
if (!requireNamespace("specifiR", quietly = TRUE)){
  devtools::install_github("gzahn/specifiR")
}

# load specifiR
library(specifiR)

# LOAD DATA ####
otu <- readRDS("./sandbox/soils_otu_low_24.rds")
groups <- readRDS("./sandbox/soils_env_low_24.rds") %>% pluck("species")
ps <- readRDS("../Biocrust_MIP/data/physeq_18S_clean.RDS")


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


# TRY PHYLOSEQ VERSION ####
out2 <- specifiR_physeq(physeq = ps,groups = "depth",n.perm = 99)

