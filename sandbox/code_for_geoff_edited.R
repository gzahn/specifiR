##load packages
library(tidyverse)
# load makecwm() function
source("./R/makecwm.R")

# using data files from Abbey to start with.
# once this code makes sense, re-run everything starting with multipatt()
# ... to get the ISA table

##read in data tables

# load data
otu_low <- read_csv("./sandbox/otu_low.csv")
names(otu_low)[1] <- "sample"
isa_low <- read_csv("./sandbox/isa_low.csv")
names(isa_low)[1] <- "otu"

# get group vector
groups <- otu_low$sample %>% substr(1,1)
# get sample vector
sample_id <- otu_low$sample
# add row names
row.names(otu_low) <- otu_low$sample
row.names(isa_low) <- isa_low$otu


# remove sample column from otu table for now
otu_low$sample <- NULL
# convert to presence-absence 'matrix'
otu_low[otu_low > 0] <- 1
otu_low <- as.matrix(otu_low)
row.names(otu_low) <- sample_id

# get occurrence values for each otu
isa_low$occurrences <- colSums(otu_low)

##remove rare taxa
# we REALLY need to automate this or give reasonable default or
# somehow allow users to select the value

isa_low_sub <- isa_low %>%
  dplyr::filter(occurrences > 3)
# add "index" column
isa_low_sub <- isa_low_sub %>%
  dplyr::mutate(index = 1-p)

# subset otu table to remaining "non-rare" taxa
otu_isa_low <- otu_low[,which(colnames(otu_low) %in% isa_low_sub$otu)]

##perform cwm analysis
# CWM should be using relative abundance values, but I think we're giving it presence-absence data

isa_low_list <- isa_low_sub
isa_low_cwm <- makecwm(otu_isa_low, isa_low_sub$index)
