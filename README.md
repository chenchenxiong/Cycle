# Cycle

**Contructing cell type-specific lncRNA-mRNA regulatory networks in ASD with Cycle**

# Introduction

The rapid development of single-nucleus RNA-sequencing (snRNA-seq) technology has revolutionized the depth and resolution, which could better capture and analyze cell-specific gene expression, and deeper understand lncRNA regulation mechanism of ASD.

To explore the functions and roles of lncRNA regulation specific to different cell types in the pathogenesis of ASD, we construct cell type-specific lncRNA regulatory networks (Cycle) at the level of single cell nucleus.


## Schematic illustration

A schematic illustration of **LCSlncR** is shown in the folowing.







# Installation

```{r echo=FALSE, results='hide', message=FALSE}
install.packages("devtools")
library(devtools)
install_github("chenchenxiong/Cycle")
```

## Quick example to use LCSlncR

For inferring cell type-specific lncRNA-mRNA regulatory networks, users should prepare matched lncRNA and mRNA snRNA-seq expression data. Users can use the following scripts to infer cell type-specific lncRNA regulation. 

```{r echo=FALSE, results='hide', message=FALSE}
## Load LCSlncR package
library(LCSlncR)

## Load the matched lncRNA and mRNA snRNA-seq expression data of three cell types (Microglia, ASTFB, Neumat) from https://drive.google.com/drive/folders/1_xjA_S77POIh28w49trPnuFUCTWagBzG
load('ASD_exp_3cell_types.rda')

microglia_Cycle_network <- Cycle_network(ASD_Microglia_ncR_data, ASD_Microglia_mR_data, boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
microglia_Cycle_networks = do.call(rbind, microglia_Cycle_network)

ASTFB_Cycle_network <- Cycle_network(ASD_ASTFB_ncR_data,ASD_ASTFB_mR_data, boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
ASTFB_Cycle_networks = do.call(rbind, ASTFB_Cycle_networks)

Neumat_Cycle_network <- Cycle_network(ASD_Neumat_ncR_data,ASD_Neumat_mR_data, boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
Neumat_Cycle_networks = do.call(rbind, Neumat_Cycle_network)

```

## License

GPL-3
