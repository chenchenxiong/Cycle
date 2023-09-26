# LCSlncR
**Inferring cell type-specific lncRNA regulation**



## Schematic illustration

A schematic illustration of **LCSlncR** is shown in the folowing.







# Installation
```{r echo=FALSE, results='hide', message=FALSE}
install.packages("devtools")
library(devtools)
install_github("chenchenxiong/LCSlncR")
```
## :zap: Quick example to use LCSlncR
For inferring cell type-specific lncRNA-mRNA regulatory networks, users should prepare matched lncRNA and mRNA snRNA-seq expression data. Users can use the following scripts to infer cell type-specific lncRNA regulation. 

```{r echo=FALSE, results='hide', message=FALSE}
## Load LCSlncR package
library(LCSlncR)

## Load prepared datasets in LCSlncR package, it includes matched lncRNA and mRNA snRNA-seq expression data of three cell types (Microglia, ASTFB, Neumat).
data(ASD_exp_3cell_types)

## Inferring cell type-specific lncRNA-mRNA regulatory networks 
data(ASD_exp_3cell_types)

microglia_LCSlncR_network = LCSlncR_net(ASD_Microglia_ncR_data, ASD_Microglia_mR_data, boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
microglia_LCSlncR_networks = do.call(rbind, microglia_LCSlncR_network)

ASTFB_LCSlncR_network = LCSlncR_net(ASD_ASTFB_ncR_data,ASD_ASTFB_mR_data, boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
ASTFB_LCSlncR_networks = do.call(rbind, ASTFB_LCSlncR_networks)

Neumat_LCSlncR_network = LCSlncR_net(ASD_Neumat_ncR_data,ASD_Neumat_mR_data, boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
Neumat_LCSlncR_networks = do.call(rbind, Neumat_LCSlncR_network)

```

## License
GPL-3
