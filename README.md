# Cycle
**Contructing cell type-specific lncRNA regulatory network in autism with Cycle**

# Introduction

As the pathogenesis of ASD involves multiple cell types and mul-tiple biological processes regulated by lncRNAs, thus it is crucial to study cell type-specific lncRNA regulation in ASD.

To explore the dynamic lncRNA regulation across different ASD cell types, we develop a novel method, Cycle, to construct cell type-specific lncRNA regulatory networks in ASD. 


## Flowchart illustration

A flowchart illustration of **Cycle** is shown in the following.

<p align="center">
  <img src="https://github.com/chenchenxiong/Cycle/blob/main/Cycle_flowchart.png" alt="Cycle flowchart illustration" border="0.1">
</p>

Firstly, we perform data integration and data imputation, and extract the matched lncRNA and mRNA expression data by using gene annotation information from HGNC. Secondly, we identify 17 cell type-specific lncRNA regulatory networks for 17 ASD cell types. Finally, we perform downstream analysis consisting of topological analysis, uniqueness analysis, multiclass classification analysis, enrichment analysis, cell communication analysis, and biomarker discovery.


# Installation
```{r echo=FALSE, results='hide', message=FALSE}
install.packages("devtools")
library(devtools)
install_github("chenchenxiong/Cycle")
```
## Quick example to use Cycle
For inferring cell type-specific lncRNA-mRNA regulatory networks, users should prepare the matched lncRNA and mRNA snRNA-seq expression data. Users can use the following scripts to infer cell type-specific lncRNA regulation. 
Load processed matched lncRNA and mRNA snRNA-seq data with ASD samples could be downloaded from https://drive.google.com/drive/folders/1_xjA_S77POIh28w49trPnuFUCTWagBzG

```{r echo=FALSE, results='hide', message=FALSE}
load('ASD_exp_3cell_types.rda')

microglia_Cycle_networks <- Cycle_network(ASD_Microglia_ncR_data[1:100,1:100], ASD_Microglia_mR_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)

ASTFB_Cycle_networks <- Cycle_network(ASD_ASTFB_ncR_data[1:100,1:100],ASD_ASTFB_mR_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)

Neumat_Cycle_networks <- Cycle_network(ASD_Neumat_ncR_data[1:100,1:100],ASD_Neumat_mR_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)

celltypes = c('Microglia', 'ASTFB', 'Neumat')
Cycle_networks_tmp = list(microglia_Cycle_networks, ASTFB_Cycle_networks, Neumat_Cycle_networks)

Cycle_networks = list()
for(i in 1:length(celltypes)){
  res_list <- Cycle_networks_tmp[[i]]
  cell_num <- length(res_list)
  overlap <- Overlap.net(net = res_list, overlap.num = round(cell_num*0.9),type = 'least')
  colnames(overlap) <- c('lncRNAs','mRNAs')
  Cycle_networks[[celltypes[i]]] <- overlap
}
```

## License
GPL-3
