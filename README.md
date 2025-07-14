# specifiR

<img src="https://github.com/gzahn/specifiR/blob/main/assets/sticker.png" alt="hex_sticker" width="200"/>

(Replace with better hex sticker...)

Implements a "Community Weighted Mean Indicator" Analysis for quantifying group specificity of entire communities.
This method generates a specificity index associated with each sample community, with more specialist communities having higher values and more generalist communities having lower values.
This can be used to compare the specificity of communities associated with groups.
For example, you can compare the specificity of microbial communities associated with plant hosts.
Abbey add further, better, description of why and how?

## Installation:

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("gzahn/specifiR")
```

## Requirements

Depends on:

- tidyverse >= 2.0.0
- future >= 1.33.2
- future.apply >= 1.11.2

## Citation

(Link to publication)



## Example usage
