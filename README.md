# MicroMCA

## 1. Install packages
```
devtools::install_github('hzaurzli/PhagesLDA')
library(PhagesLDA)
```
**We recommend using R 4.0+ for operation**

## 2. MCA dimensionality reduction
***Step 1: Loading data***
```
data(data)
load(data)
```

Species abundance matrix or genes abundance matrix
***For species abundance matrix:***
| |  sample_1  |  sample_2  |  sample_3  |  sample_4 | sample_5 | sample_6 | ... | sample_m |
|  ----  | ----  |  ----  | ----  |  ----  | ----  | ---- | ---- | ---- |
| otu_1  | 8 | 9 | 11 | 4 | 2 | 9 | ... | 5 |
| otu_2  | 5 | 9 | 9 | 2 | 3 | 3 | ... | 6 |
| otu_3  | 12 | 6 | 17 | 16 | 16 | 17 | ... | 1 |
| otu_4  | 17 | 18 | 16 | 15 | 1 | 2 | ... | 6 |
| otu_5  | 15 | 16 | 14 | 2 | 6 | 5 | ... | 2 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| otu_n  | 7 | 8 | 5 | 2 | 6 | 17 | ... | 15 |

Rows of example data represent different otu and each column represents a different samples, count file

***For genes abundance matrix:***
| |  sample_1  |  sample_2  |  sample_3  |  sample_4 | sample_5 | sample_6 | ... | sample_m |
|  ----  | ----  |  ----  | ----  |  ----  | ----  | ---- | ---- | ---- |
| gene_1  | 8 | 9 | 11 | 4 | 2 | 9 | ... | 5 |
| gene_2  | 5 | 9 | 9 | 2 | 3 | 3 | ... | 6 |
| gene_3  | 12 | 6 | 17 | 16 | 16 | 17 | ... | 1 |
| gene_4  | 17 | 18 | 16 | 15 | 1 | 2 | ... | 6 |
| gene_5  | 15 | 16 | 14 | 2 | 6 | 5 | ... | 2 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| gene_n  | 7 | 8 | 5 | 2 | 6 | 17 | ... | 15 |

Rows of example data represent different genes and each column represents a different samples, count file
