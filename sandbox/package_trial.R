# trial of functions

# load dependencies (for this script)
library(tidyverse)
library(phyloseq)

# install specifiR if not already installed
if (!requireNamespace("specifiR", quietly = TRUE)){
  devtools::install_github("gzahn/specifiR")
}

# load specifiR
library(specifiR)

# LOAD DATA ####
otu <- readRDS("./data/soils_otu_low_24.rds")
groups <- readRDS("./data/soils_env_low_24.rds") %>% pluck("species")
ps <- readRDS("./data/example_physeq.RDS")

# TRY FUNCTION ####
out <- specifiR(comm = otu,
                groups = groups,
                seed = 123,
                n.perm = 99,
                pval.cutoff = 0.05,
                max.ratio = 0,
                ovp.plot = TRUE)

out$community_specificity_index
out$taxon_specificity_index
out$isa_results
out$process_summary

# TRY PHYLOSEQ VERSION ####
out2 <- specifiR_physeq(physeq = ps,groups = "invasion",n.perm = 99,ovp.plot=TRUE)
out2$process_summary

