# This is a trial walkthrough starting from the otu_low data set

# We need to consider:
# 1. What data object will people be bringing to this package?
#    - An OTU table and a previously done multipatt object (they do that themselves, can set custom stuff)
#    - A phyloseq object (make it easy for lots of users)
#    - OTU table and a grouping vector (package agnostic)
#    - All 3 of the above (we would have 3 functions for import, then second function for rest of steps)
# 2. How can we automate the selection of a "rare taxon" cutoff?
#    - You chose >3, but that's based on eyeballing a graph
#    - What would a good default be? Can we base it on data properties?
#    - Or is it just the taxa that don't have valid p-values?
# 3. Doesn't CWM use relative abundance data?
#    - We're giving it presence-absence data though...
#    - Users could submit raw count data, we could do transformations for them
#    - Like, we would generate both PA and RA versions of their matrix


# SETUP ####
library(tidyverse) # we can re-write to remove this dependency eventually
library(indicspecies) # we can probably just include the multipatt() function and remove this dependency

# load makecwm() function
source("./R/makecwm.R") # taken from the ecole package (deprecated on CRAN)



# LOAD DATA ####

# load data
otu <- read_csv("./sandbox/otu_low.csv")
# fix missing first colname
names(otu)[1] <- "sample"
# pull out sample ids
sample_id <- otu[["sample"]]
# get grouping variable for indicator species analysis
# here, I'm just assuming this info comes from sample names
# we will need to have a way for users to indicate a grouping variable
groups <- substr(sample_id,1,1)
# convert to matrix and add row names
otu[["sample"]] <- NULL # remove character column
otu <- as(otu,'matrix')
row.names(otu) <- sample_id


# RUN MULTIPATT() ####
indval <- multipatt(otu, groups,
                     control = how(nperm=99)) # speed up w/ low permutations




# PREP FOR CWM ####

# extract p values
isa <- indval[["sign"]]
isa[["sample_id"]] <- row.names(isa)
# convert to presence-absence
otu[otu>0] <- 1

# get occurrence values
isa[["occurrence"]] <- colSums(otu)




# REMOVE RARE TAXA ####

# need to rethink how to make this non-arbitrary
# how to determine this value from data
# how to let users decide in function parameters
rare.cutoff <- 3

# find non-rare taxa
nonrare_otus <- isa[["sample_id"]][isa[["occurrence"]] > rare.cutoff]

# subset isa
isa_nonrare <- isa[isa[["sample_id"]] %in% nonrare_otus,]

# add index column
# multipatt() returns df with "index" column already
isa_nonrare[["taxon_index"]] <- (1 - isa_nonrare[["p.value"]])

# subset otu table to match
otu_nonrare <- otu[,colnames(otu) %in% nonrare_otus]

# convert NA values to 0
# need to think about this...
# this goes back to how to automate "rare" taxa removal...
# should taxa that can't get a valid pvalue just be removed???
isa_nonrare[["taxon_index"]][is.na(isa_nonrare[["taxon_index"]])] <- 0

# PERFORM CWM ANALYSIS ####
cwm <- makecwm(otu_nonrare, isa_nonrare[["taxon_index"]])

# clean it up to make a more usable object
output <- data.frame(sample_id = row.names(cwm),
           cwm = cwm[["V1"]])
output
