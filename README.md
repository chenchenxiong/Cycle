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
## Load required R packages
library(Cycle)
library(igraph)

## Load prepared datasets
load('ASD_exp_3cell_types.rda')

# Identificaton of cell type-specific lncRNA regulation
## lncRNA regulation of microglia
microglia_Cycle_networks <- Cycle_network(ASD_Microglia_ncR_data[1:100,1:100], ASD_Microglia_mR_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)

## lncRNA regulation of ASTFB
ASTFB_Cycle_networks <- Cycle_network(ASD_ASTFB_ncR_data[1:100,1:100],ASD_ASTFB_mR_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)

## lncRNA regulation of Neumat
Neumat_Cycle_networks <- Cycle_network(ASD_Neumat_ncR_data[1:100,1:100],ASD_Neumat_mR_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)

## Integrating cell type-specific lncRNA-mRNA regulatory network
celltypes = c('Microglia', 'ASTFB', 'Neumat')
Cycle_networks = list(microglia_Cycle_networks, ASTFB_Cycle_networks, Neumat_Cycle_networks)

# Identifying cell type-specific hub lncRNAs 
Cycle_hubs_lncRNAs <- hub_discovery(Cycle_networks)
names(Cycle_hubs_lncRNAs) = names(Cycle_networks)

# Discovering stable and rewired lncRNA regulation
## Given cell type-specific lncRNA-mRNA regulatory networks, we could discover stable and rewired lncRNA-mRNA regulatory networks.
Overlap_net.all <- Overlap.net(Cycle_networks, overlap.num = 1, type = "least") 
Overlap_net.1 <- Overlap.net(Cycle_networks, overlap.num = 1, type = "equal")
Overlap_net.2 <- Overlap.net(Cycle_networks, overlap.num = 2, type = "least")
Overlap_net.5 <- Overlap.net(Cycle_networks, overlap.num = 5, type = "least")
Overlap_net.9 <- Overlap.net(Cycle_networks, overlap.num = 9, type = "least")
Overlap_net.15 <- Overlap.net(Cycle_networks, overlap.num = 15, type = "least")
Overlap_net.17 <- Overlap.net(Cycle_networks, overlap.num = 17, type = "equal")

Overlap_Cycle_network_rewired <- Overlap_net.1
Overlap_Cycle_network_stable <- Overlap_net.15

## Given cell type-specific lncRNA-mRNA regulatory networks, we could discover stable and rewired hub lncRNAs.
Overlap.hub.all <- Overlap.hub(Cycle_hubs_lncRNAs, overlap.num = 1, type = "least") 
Overlap.hub.1 <- Overlap.hub(Cycle_hubs_lncRNAs, overlap.num = 1, type = "equal")
Overlap.hub.2 <- Overlap.hub(Cycle_hubs_lncRNAs, overlap.num = 2, type = "least")
Overlap.hub.5 <- Overlap.hub(Cycle_hubs_lncRNAs, overlap.num = 5, type = "least")
Overlap.hub.9 <- Overlap.hub(Cycle_hubs_lncRNAs, overlap.num = 9, type = "least")
Overlap.hub.15 <- Overlap.hub(Cycle_hubs_lncRNAs, overlap.num = 15, type = "least")
Overlap.hub.17 <- Overlap.hub(Cycle_hubs_lncRNAs, overlap.num = 17, type = "equal")

Overlap_Cycle_hubs_rewired <- Overlap.hub.1
Overlap_Cycle_hubs_stable <- Overlap.hub.15

# Uniqueness of cell type-specific lncRNA regulation
##　Uniqueness of　lncRNA-mRNA regulatory networks across cell types. 
Cycle_networks_sim <- Sim.network(Cycle_networks, Cycle_networks, directed = TRUE)
Cycle_networks_unique <- 1-Cycle_networks_sim

## Uniqueness of hub lncRNAs across cell types. 
Cycle_hub_Sim <- Sim.hub(Cycle_hubs_lncRNAs, Cycle_hubs_lncRNAs)
Cycle_hub_unique <- 1-Cycle_hub_Sim

#　Cell similarity network in terms of network similarity matrix between 17 cell types
## lncRNA-mRNA regulatory networks
Cycle_network_adjacency_matrix <- ifelse(Cycle_networks_sim > median(Cycle_networks_sim[lower.tri(Cycle_networks_sim)]), 1, 0)
diag(Cycle_network_adjacency_matrix) <- 0
colnames(Cycle_network_adjacency_matrix) <- rownames(Cycle_network_adjacency_matrix) <- names(Cycle_networks)
Cycle_network_adjacency_matrix_graph <- graph_from_adjacency_matrix(Cycle_network_adjacency_matrix, mode = "undirected")

## hub lncRNAs 
Cycle_hub_adjacency_matrix <- ifelse(Cycle_hub_Sim > median(Cycle_hub_Sim[lower.tri(Cycle_hub_Sim)]), 1, 0)
diag(Cycle_hub_adjacency_matrix) <- 0
colnames(Cycle_hub_adjacency_matrix) <- rownames(Cycle_hub_adjacency_matrix) <- names(Cycle_networks)
Cycle_hub_adjacency_matrix_graph <- graph_from_adjacency_matrix(Cycle_hub_adjacency_matrix, mode = "undirected")

```

## License
GPL-3
