# pagaR
a R version PAGA

origin code:
* https://github.com/theislab/paga

pagaR package can be easily installed from Github using devtools:
```
devtools::install_github("CaoWei-UM/pagaR")
```

How to use:
```
library(Seurat)
library(pagaR)
load(pbmc) #pbmc is a seurat object
CalculatePAGA(pbmc,group.by='celltype')
PlotPAGA(pbmc,group.by='celltype')
```
