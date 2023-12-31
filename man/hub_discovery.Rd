% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cycle.R
\name{hub_discovery}
\alias{hub_discovery}
\title{Inferring hub lncRNAs}
\usage{
hub_discovery(net)
}
\arguments{
\item{net}{List object, the list of networks.}
}
\value{
List object, the list of hub lncRNAs.
}
\description{
Inferring hub lncRNAs
}
\examples{
# NOT RUN
# load('ASD_exp_3cell_types.rda')
# microglia_Cycle_networks <- Cycle_network(ASD_Microglia_ncR_data[1:100,1:100], ASD_Microglia_mR_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
# microglia_hubs_lncRNAs <- hub_discovery(microglia_Cycle_networks)
}
