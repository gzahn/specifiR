##load packages
library(tidyverse)
library(ecole)

##read in data tables
otu_low <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soils_otu_low_24.rds")
isa_low <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\isa_soils_low_24.csv")

##LOW##
isa_low <- (isa_low %>% 
  mutate(across(
    everything(),
    ~ map_chr(.x, ~ gsub("\"", "", .x))
  )))

otu_low[otu_low > 0] = 1
occur <- as.data.frame(colSums(otu_low))

isa_low <- column_to_rownames(isa_low, var = "out.id")
isa_low$p <- as.numeric(isa_low$p)
isa_low$value.IV <- as.numeric(isa_low$value.IV)
isa_low_merge <- merge(isa_low, occur, by = 0)
isa_low_merge$occurances <- isa_low_merge$`colSums(otu_low)`
isa_low_merge <- subset(isa_low_merge, select = -c(7))

p_low <- ggplot(isa_low_merge, aes(x = occurances, y = p)) + geom_point(aes(color = value.IV)) + geom_vline(xintercept = 4) + geom_hline(yintercept = 0.05) + scale_x_continuous(breaks = seq(0, 60, by=5))
p_low

##remove rare taxa
isa_low_sub <- subset(isa_low_merge, occurances > 3)
isa_low_sub$index <- 1 - isa_low_sub$p
saveRDS(isa_low_sub, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\isa_low_sub_soils.rds" )

##perform cwm analysis
isa_low_list <- isa_low_sub$Row.names
isa_low_sub <- rownames_to_column(isa_low_sub, "numbers")
isa_low_sub <- column_to_rownames(isa_low_sub, "Row.names")
otu_isa_low <- otu_low[, (colnames(otu_low) %in% isa_low_list)]
isa_low_cwm <- makecwm(otu_isa_low, isa_low_sub$index)
