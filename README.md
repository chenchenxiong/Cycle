# Cycle
**Modelling cell type-specific lncRNA regulatory network in autism with Cycle**

# Introduction

Autism spectrum disorder (ASD) is a class of complex neurodevelopment disorders with high genetic heterogeneity. Long non-coding RNAs (lncRNAs), as vital regulators, exert specific functions in different cell types and play pivotal roles in neurological diseases including ASD. Therefore, studying specific lncRNA regulation in various cell types is crucial for deciphering ASD molecular mechanisms. Existing computational methods utilize bulk transcriptomics data across all of cells or samples, which could reveal the commonalities, but ignore the specificity of lncRNA regulation across various cell types. Here, we present Cycle (Cell type-specific lncRNA regulatory network) to construct the landscape of cell type-specific lncRNA regulation in ASD. We find that each ASD cell type is unique in lncRNA regulation, and over one-third and all of cell type-specific lncRNA regulatory networks are scale-free and small-world networks, respectively. Across 17 ASD cell types, we discover that 19 rewired and 11 conserved modules, and eight rewired and three conserved hubs underlying in the identified cell type-specific lncRNA regulatory networks, are significantly enriched in ASD-related terms. Furthermore, more similar ASD cell types display stronger connection in the constructed cell similarity network. Finally, the comparison results demonstrate that Cycle is a potential method in uncovering cell type-specific lncRNA regulation.


## Flowchart illustration

A flowchart illustration of **Cycle** is shown in the following.

<p align="center">
  <img src="https://github.com/chenchenxiong/Cycle/blob/main/Fig.%201.%20Workflow%20of%20Cycle.jpg" alt="Cycle flowchart illustration" border="0.1">
</p>

Fig. 1. Workflow of Cycle. Firstly, Cycle extracts the matched lncRNA and mRNA expression data by using gene annotation information from HGNC (HUGO gene Nomenclature Committee), and further retain the highly expressed lncRNAs and mRNAs for each cell type. In total, we have obtained 17 cell type-specific expression data of highly expressed lncRNAs and mRNAs. Secondly, Cycle models cell type-specific lncRNA regulatory networks for 17 ASD cell types. Furthermore, Cycle identifies the rewired and conserved modules, and infers hubs based on the constructed cell type-specific lncRNA regulatory networks. Finally, Cycle conducts four types of downstream analyses, including modules identification, hub inference, network topological analysis, uniqueness analysis, cell similarity network construction, and enrichment analysis.

# Installation
```{r echo=FALSE, results='hide', message=FALSE}
install.packages("devtools")
library(devtools)
install_github("chenchenxiong/Cycle")
library(Cycle)
```
## Quick example to use Cycle
For inferring cell type-specific lncRNA-mRNA regulatory networks, users should prepare the matched lncRNA and mRNA snRNA-seq expression data. Users can use the following scripts to infer cell type-specific lncRNA regulation. 
Load processed matched lncRNA and mRNA snRNA-seq data with ASD samples could be downloaded from https://drive.google.com/drive/folders/1_xjA_S77POIh28w49trPnuFUCTWagBzG

```{r echo=FALSE, results='hide', message=FALSE}
data(ASD_exp_3cell_types.rda)

microglia_Cycle_networks <- Cycle_network(Microglia_ASD_lncRNAs_data[1:100,1:100], Microglia_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)

ASTFB_Cycle_networks <- Cycle_network(ASTFB_ASD_lncRNAs_data[1:100,1:100],ASTFB_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)

Neumat_Cycle_networks <- Cycle_network(Neumat_ASD_lncRNAs_data[1:100,1:100],Neumat_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)

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
