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
library(future)
library(future.apply)
plan(multisession)  # cross-platform parallelism

# load makecwm() function
source("./R/makecwm.R") # taken from the ecole package (deprecated on CRAN)
source("./R/first_consecutive.R")
source("./R/calc_ivmax.R")
'%ni%' <- Negate("%in%")

# user options ...
set.seed(666)
n_perm <- 999
pval.cutoff <- 0.05
max.ratio <- 0

# if users have physeq object, use that...
# if users have otu table and grouping vector, use that!


# LOAD DATA ####
otu <- readRDS("./sandbox/soils_otu_low_24.rds")
groups <- readRDS("./sandbox/soils_env_low_24.rds") %>% pluck("species")



# RE-ENGINEER PC-ORD ####

## Step 1. Calculate relative abundance of each taxon in each group ####
otu_sums <- colSums(otu)
group_sums <- data.frame(otu_name = colnames(otu))
otu$group <- groups
i <- "tshe"

for(i in unique(groups)){
  otu_subset <- otu %>% dplyr::filter(group %in% i)
  otu_subset$group <- NULL
  otu_sums <- colSums(otu_subset)
  group_sums[[i]] <- otu_sums
}

df_rel_abund <- group_sums %>%
  rowwise() %>%
  mutate(across(where(is.numeric), ~ .x / sum(c_across(where(is.numeric))))) %>%
  ungroup()

## Step 2. calculate fidelity ####
# proportion of samples within a given group that each taxon is found in

# divide up by group
# count number of samples in each group
# count number of samples within each group that a given taxon is present
# divide for proportion

# make presence-absence otu table
otu_pa <-
  otu %>%
  mutate(across(where(is.numeric), ~ if_else(.x > 0, 1, 0)))

group_fidelity <- data.frame(otu_name = colnames(otu_pa %>% dplyr::select(-group)))

i <- "tshe"

for(i in unique(groups)){
  otu_subset <- otu_pa %>% dplyr::filter(group %in% i)
  otu_subset$group <- NULL
  otu_sums <- colSums(otu_subset)
  group_fidelity[[i]] <- otu_sums / nrow(otu_subset)
}

## Step 3. multiply group_sums by group_fidelity ####
df_rel_abund
group_fidelity

# make matrices and add row names
group_sum_mat <- as.matrix(df_rel_abund %>% dplyr::select(-otu_name))
row.names(group_sum_mat) <- group_sums$otu_name
group_fidelity_mat <- as.matrix(group_fidelity %>% dplyr::select(-otu_name))
row.names(group_fidelity_mat) <- group_fidelity$otu_name

# multiply matrices
indicator_values <- group_sum_mat * group_fidelity_mat
indicator_values <- indicator_values * 100

# step 4. identify the highest indicator value (iv.max)
# of the N columns, get the biggest one
# for package, return both of these matrices (optional)


iv.max <-
  indicator_values %>%
  apply(1, max)
# make data frame with all indicator values and iv.max
indicator_values <- bind_cols(indicator_values,iv.max=iv.max)
# add back otu names
indicator_values$otu_name <- names(iv.max)

## Step 4. monte carlo permutations ####

# Precompute matrices
otu_numeric <- otu %>% dplyr::select(where(is.numeric))
otu_matrix <- as.matrix(otu_numeric)
otu_pa <- (otu_matrix > 0) * 1
sample_groups <- groups
group_levels <- unique(sample_groups)



# Observed iv.max
iv.max_obs <- calc_ivmax(sample_groups)

# Permuted iv.max values
iv.max_perm <- future_replicate(n_perm, {
  calc_ivmax(sample(sample_groups))
})

# Empirical p-values
iv.pval <- rowMeans(iv.max_perm >= iv.max_obs)

# Output result
indicator_results <- tibble(
  otu_name = colnames(otu_matrix),
  iv.max = iv.max_obs,
  p.value = iv.pval
)


# PREP FOR CWM ####

## convert to presence-absence ####
otu_pa <- otu %>% dplyr::select(-group)
otu_pa[otu_pa>0] <- 1

## get occurrence values ####
indicator_results[["occurrence"]] <- colSums(otu_pa)

# REMOVE RARE TAXA ####

occurrence_groups <- factor(indicator_results[["occurrence"]])
indicator_results[["occurrence_groups"]] <- occurrence_groups

# make ratio dataframe that connects Nsites groups to ratio of sig/non-sig pvalues
ratio_df <-
  indicator_results %>%
  mutate(significant = p.value <= 0.05) %>%
  group_by(occurrence_groups) %>%
  summarize(ratio = sum(significant) / sum(!significant))

# find which taxa to remove (first N taxa at or below max.ratio)
to_remove <- which(ratio_df$ratio <= max.ratio) %>% first_consecutive()
# subset indicator results
isa_subset <-
  indicator_results %>%
  dplyr::filter(occurrence %ni% to_remove)

# subset otu table to match
otu_subset <-
  otu[,colnames(otu) %in% isa_subset$otu_name]


# calculate taxon index (so bigger indicates more indicative)
isa_subset$taxon_index <- 1 - isa_subset$p.value

# PERFORM CWM ANALYSIS ####
cwm <- makecwm(otu_subset, isa_subset[["taxon_index"]])

# clean it up to make a more usable object
output <- data.frame(sample_id = row.names(cwm),
           cwm = cwm[["V1"]])

# add taxon index to results
indicator_results$taxon_index <- 1 - indicator_results$p.value

# create output object (list)
out <- list(community_specificity_index = output,
            taxon_specificity_index = indicator_results,
            isa_results = indicator_values)

# return "out"

