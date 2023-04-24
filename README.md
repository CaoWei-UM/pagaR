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
library(pagaR)
load(pbmc)
CalculatePAGA(pbmc,group.by='celltype')
PlotPAGA(pbmc,group.by='celltype')
```
