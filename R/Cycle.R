#' Define the upper and lower neighborhood of each plot
#'
#' @param gx  Gene expression values of gene x (a vector) in n cells
#' @param gy  Gene expression values of gene y (a vector) in n cells
#' @param boxsize Size of neighborhood (0.1 in default)
#' @importFrom pracma eps
#' @import stats
#' @return upperlower_result: a list of upper result and lower result
#' @export
#'
#' @references
#' Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.,22(1):578.
#' Dai H, Li L, Zeng T, Chen L. Cell-specific network constructed by single-cell RNA sequencing data. Nucleic Acids Res. 2019 Jun 20;47(11):e62. doi: 10.1093/nar/gkz172.
upperlower <-  function(gx, gy, boxsize = 0.1) {
  n <- length(gx)
  upper <- zeros(2, n)
  lower <- zeros(2, n)
  for (i in seq_len(2)) {
    g <- gx * (i == 1) + gy * (i == 2)
    s1 <- sort(g, index.return = TRUE)[[1]]
    s2 <- sort(g, index.return = TRUE)[[2]]
    n0 <- n - sum(sign(s1))
    h <- round(boxsize / 2 * sum(sign(s1)) + eps(1))
    k <- 1
    while (k <= n) {
      s <- 0
      while ((n >= k + s + 1) && (s1[k + s + 1] == s1[k])) {
        s <- s + 1
      }
      if (s >= h) {
        upper[i, s2[k:(k + s)]] <- g[s2[k]]
        lower[i, s2[k:(k + s)]] <- g[s2[k]]
      } else {
        upper[i, s2[k:(k + s)]] <- g[s2[min(n, k + s + h)]]
        lower[i, s2[k:(k + s)]] <- g[s2[max(n0 * (n0 > h) + 1, k - h)]]
      }
      k <- k + s + 1
    }
  }
  upperlower_result <- list(upper = upper, lower = lower)
  return(upperlower_result)
}


#' Define the upper and lower neighborhood of each plot using local standard deviation
#'
#' @param gene1 Gene expression values of gene1 (a vector) in n cells
#' @param gene2 Gene expression values of gene2 (a vector) in n cells
#' @param boxsize Size of neighborhood (0.1 in default)
#' @param iteration Whether need iterate (TRUE in default)
#' @param cell_id cell numbers, NULL means all cell to compute (NULL in default)
#' @param maxiter Maximum Number Of Iterations (20 in default)
#' @importFrom pracma eps
#' @import stats
#' @return upperlower_devresult: a list of upper result and lower result
#' @export
#' @references
#' Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.,22(1):578.
#' Dai H, Li L, Zeng T, Chen L. Cell-specific network constructed by single-cell RNA sequencing data. Nucleic Acids Res. 2019 Jun 20;47(11):e62. doi: 10.1093/nar/gkz172.
#' Wang X, Choi D, Roeder K. Constructing local cell-specific networks from single-cell data. Proc Natl Acad Sci U S A. 2021;118(51):e2113178118. doi:10.1073/pnas.2113178118
upperlower_dev <- function(gene1,
                          gene2,
                          boxsize = 0.1,
                          iteration = TRUE,
                          cell_id = NULL,
                          maxiter = 20) {
  if (length(gene1) == length(gene2)) {
    n1 <- 2
    n2 <- length(gene1)
    if (is.null(cell_id)) {
      cell_id <- seq_len(n2)
    }
    upperlower_result <-  upperlower(gene1, gene2, boxsize = 0.1)
    up_q <- upperlower_result$upper
    low_q <- upperlower_result$lower
    upper <- zeros(2, n2)
    lower <- zeros(2, n2)
    y <- 0
    md <- 1

    if (iteration) {
      maxiter <- maxiter
      for (k in cell_id) {
        if (gene1[k] * gene2[k] > 0) {
          d2_0 <- md * sd(gene2[gene1 <= up_q[1, k] & gene1 >= low_q[1, k]])
          d2_0 <- ifelse(is.na(d2_0), 0, d2_0)

          d1_0 <- md * sd(gene1[gene2 <= up_q[2, k] &
                                 gene2 >= low_q[2, k]])
          d1_0 <- ifelse(is.na(d1_0), 0, d1_0)

          len1_1 <- length(gene2[(gene1 <= gene1[k] + d1_0) &
                                  (gene1 >= gene1[k] - d1_0)])
          len1_2 <- length(gene1[gene2 <= (gene2[k] + d2_0) &
                                  gene2 >= (gene2[k] - d2_0)])


          d2_1 <- md * sd(gene2[gene1 <= (gene1[k] + d1_0) &
                                 gene1 >= (gene1[k] - d1_0)])
          d2_1 <- ifelse(is.na(d2_1), 0, d2_1)

          d1_1 <- md * sd(gene1[gene2 <= (gene2[k] + d2_0) &
                                 gene2 >= (gene2[k] - d2_0)])
          d1_1 <- ifelse(is.na(d1_1), 0, d1_1)

          enter_data <- sqrt((d2_0 - d2_1) ^ 2 + (d1_0 - d1_1) ^ 2)
          enter <- (sqrt((d2_0 - d2_1) ^ 2 + (d1_0 - d1_1) ^ 2) < 10 ^ (-5))
          count <- 0
          while ((sqrt((d2_0 - d2_1) ^ 2 + (d1_0 - d1_1) ^ 2) < 10 ^ (-5)) &&
                 (count < maxiter)) {
            d2_0 <- d2_1
            d1_0 <- d1_1
            len2_1 <- len1_1
            len2_2 <- len1_2
            d2_1 <- md * sd(gene2[gene1 <= (gene1[k] + d1_0) &
                                   gene1 >= (gene1[k] - d1_0)])
            d2_1 <- ifelse(is.na(d2_1), 0, d2_1)

            d1_1 <- md * sd(gene1[gene2 <= (gene2[k] + d2_0) &
                                   gene2 >= (gene2[k] - d2_0)])
            d1_1 <- ifelse(is.na(d1_1), 0, d1_1)

            len1_1 <- length(gene2[(gene1 <= gene1[k] + d1_0) &
                                    (gene1 >= gene1[k] - d1_0)])
            len1_2 <- length(gene1[gene2 <= (gene2[k] + d2_0) &
                                    gene2 >= (gene2[k] - d2_0)])
            if ((len2_1 == len1_1) || (len2_2 == len1_2)) {
              y <- y + 1
              break
            }
            count <- count + 1
          }
          if (count >= maxiter) {
            message('Iteration at cell ', k, ' exceeds ', maxiter)
            next
          }
          upper[1, k] <- gene1[k] + d1_1
          upper[2, k] <- gene2[k] + d2_1
          lower[1, k] <- gene1[k] - d1_1
          lower[2, k] <- gene2[k] - d2_1
        }
      }
    } else {
      for (k in cell_id) {
        if (gene1[k] * gene2[k] > 0) {
          d2 <- md * sd(gene2[gene1 <= up_q[1, k] & gene1 >= low_q[1, k]])
          d2 <- ifelse(is.na(d2), 0, d2)

          d1 <- md * sd(gene1[gene2 <= up_q[2, k] &
                               gene2 >= low_q[2, k]])
          d1 <- ifelse(is.na(d1), 0, d1)

          upper[1, k] <- gene1[k] + d1
          upper[2, k] <- gene2[k] + d2
          lower[1, k] <- gene1[k] - d1
          lower[2, k] <- gene2[k] - d2
        }
      }
    }
  }
  upperlower_devresult <- list(upper = upper, lower = lower)
  return(upperlower_devresult)
}



#' Cell-specific network of edge gx-gy
#'
#' @param gx  Gene expression values of gene x (a vector) in n cells
#' @param gy  Gene expression values of gene y (a vector) in n cells
#' @param boxsize Size of neighborhood (0.1 in default)
#' @param dev Whether to use local standard deviation (TRUE in default)
#' @param iteration Whether need iterate (TRUE in default)
#' @param cell_id cell numbers, NULL means all cell to compute (NULL in default)
#' @param maxiter Maximum Number Of Iterations (20 in default)
#' @importFrom pracma eps
#' @return res a vector, the normalized statistic of edge gx-gy in n cells
#' @export
#'
#' @references
#' Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.,22(1):578.
#' Dai H, Li L, Zeng T, Chen L. Cell-specific network constructed by single-cell RNA sequencing data. Nucleic Acids Res. 2019 Jun 20;47(11):e62. doi: 10.1093/nar/gkz172.
csn_edge <-
  function(gx,
           gy,
           boxsize = 0.1,
           dev = TRUE,
           iteration = TRUE,
           cell_id = NULL,
           maxiter = 20) {
    # Define the neighborhood of each plot
    if (dev) {
      upperlower_result <- upperlower_dev(gx,gy,boxsize = boxsize,iteration = iteration,cell_id = cell_id,maxiter = maxiter)
    } else {
      upperlower_result <- upperlower(gx, gy, boxsize = boxsize)
    }
    upper <- upperlower_result$upper
    lower <- upperlower_result$lower

    n <- length(gx)
    a <- zeros(2, n)
    B <- list()

    for (i in seq_len(2)) {
      g <- gx * (i == 1) + gy * (i == 2)
      upper_tmp <- upper[i, ]
      lower_tmp <- lower[i, ]
      B[[i]] <-
        (do.call(cbind, lapply(seq_len(n), function(i)
          g <= upper_tmp[i]))) &
        (do.call(cbind, lapply(seq_len(n), function(i)
          g >= lower_tmp[i])))
      a[i, ] <- colSums(B[[i]])
    }
    # Calculate the normalized statistic of edge gx-gy
    res <-(colSums(B[[1]] & B[[2]]) * n - a[1, ] * a[2, ]) / sqrt(a[1, ] * a[2, ] * (n - a[1, ]) * (n - a[2, ]) / (n - 1) + eps(1))
    return(res)
  }


#' Constructing cell type-specific lncRNA-mRNA regulatory network
#'
#' @param lncR Gene expression values of lncRNAs in single cells, rows are cells and columns are lncRNAs
#' @param mR Gene expression values of mRNAs in single cells, rows are cells and columns are mRNAs
#' @param boxsize Size of neighborhood (0.1 in default)
#' @param p.value.cutoff cutoff of p.value (0.05 in default)
#' @param num.cores number of parallel cores (2 in default)
#' @param dev Whether to use local standard deviation (TRUE in default)
#' @param iteration Whether need iterate (TRUE in default)
#' @param cell_id cell numbers, NULL means all cell to compute (NULL in default)
#' @param maxiter Maximum Number Of Iterations (20 in default)
#' @import foreach
#' @importFrom foreach %dopar%
#' @importFrom pracma eps
#' @import parallel
#' @import doParallel
#' @import stats
#'
#' @return res_list: a list of cell type-specific lncRNA-mRNA regulatory network
#' @export
#'
#' @examples
#' # NOT RUN
#' # load('ASD_exp_3cell_types.rda')
#' # microglia_Cycle_networks <- Cycle_network(Microglia_ASD_lncRNAs_data[1:100,1:100], Microglia_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # ASTFB_Cycle_networks <- Cycle_network(ASTFB_ASD_lncRNAs_data[1:100,1:100],ASTFB_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # Neumat_Cycle_networks <- Cycle_network(Neumat_ASD_lncRNAs_data[1:100,1:100],Neumat_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' @references
#' Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.,22(1):578.
Cycle_network <-
  function(lncR,
           mR,
           boxsize = 0.1,
           p.value.cutoff = 0.05,
           num.cores = 2,
           dev = TRUE,
           iteration = TRUE,
           cell_id = NULL,
           maxiter = 20) {

    lncR = t(lncR)
    mR = t(mR)
    lncRs_num <- ncol(lncR)
    mRs_num <- ncol(mR)
    cell_num <- nrow(lncR)
    int_num <- lncRs_num * mRs_num
    index <- matrix(NA, nrow = int_num, ncol = 2)

    for (i in seq(lncRs_num)) {
      for (j in seq(mRs_num)) {
        index[(i - 1) * mRs_num + j, 1] <- i
        index[(i - 1) * mRs_num + j, 2] <- j
      }
    }

    cl <- makeCluster(num.cores)
    registerDoParallel(cl)
    interin <- foreach::foreach(
      m = seq(1, int_num, 1),
      .packages = c("pracma"),
      .export = c("csn_edge", "upperlower_dev", "upperlower")
    ) %dopar% {
      csn_edge(
        lncR[, index[m, 1]],
        mR[, index[m, 2]],
        boxsize = boxsize,
        dev = dev,
        iteration = iteration,
        cell_id = cell_id,
        maxiter = maxiter
      )
    }
    # shut down the workers
    stopImplicitCluster()
    stopCluster(cl)

    interin <- do.call(rbind, interin)

    res <- matrix(NA, nrow = int_num, ncol = cell_num + 2)
    for (i1 in seq(int_num)) {
      res[i1, 1] <- colnames(lncR)[index[i1, 1]]
      res[i1, 2] <- colnames(mR)[index[i1, 2]]
    }
    res[, 3:(cell_num + 2)] <- interin
    q <- -qnorm(p.value.cutoff)
    res_list <-
      lapply(seq(cell_num), function(i2)
        res[which(as.numeric(res[, i2 + 2]) > q), seq(2)])

    return(res_list)
  }



#' Inferring hub lncRNAs
#'
#' @param net List object, the list of networks.
#'
#' @return List object, the list of hub lncRNAs.
#' @import igraph
#' @import stats
#' @export
#'
#' @examples
#' # NOT RUN
#' # load('ASD_exp_3cell_types.rda')
#' # microglia_Cycle_networks <- Cycle_network(Microglia_ASD_lncRNAs_data[1:100,1:100], Microglia_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # microglia_hubs_lncRNAs <- hub_discovery(microglia_Cycle_networks)
hub_discovery <- function(net){

  if(class(net)!= "list") {
    stop("Please check your input network! The input network should be list object! \n")
  }
  hubs_lncRNAs <- lapply(seq(net),function(i){

    network_graph <- make_graph(c(t(net[[i]][,1:2])),directed=TRUE)
    network_degree <- degree(network_graph,mode='out')
    network_degree <- network_degree[network_degree>0]

    n <- length(network_degree)
    m <- nrow(net[[i]])
    p1 <- m/(n^2)
    lamda <- n*p1
    p <- ppois(network_degree - 1, lambda = lamda, lower=FALSE)
    hub <- p[which(p < 0.05)]
    hub_names <- names(hub)
  })
  return(hubs_lncRNAs)
}


#' Identifying the overlap between multiple networks.
#'
#' @title Overlap.net
#' @param net List object, the list of networks.
#' @param overlap.num The minimum number of interactions existing in multiple networks.
#' @param type The overlapped interactions in overlap.num networks ("equal") or at least overlap.num networks ("least").
#' @importFrom stringr str_split_fixed
#' @export
#' @return Matrix object: The overlapped interactions.
#'
#' @examples
#' # NOT RUN
#' # load('ASD_exp_3cell_types.rda')
#' # microglia_Cycle_networks <- Cycle_network(Microglia_ASD_lncRNAs_data[1:100,1:100], Microglia_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # ASTFB_Cycle_networks <- Cycle_network(ASTFB_ASD_lncRNAs_data[1:100,1:100],ASTFB_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # Neumat_Cycle_networks <- Cycle_network(Neumat_ASD_lncRNAs_data[1:100,1:100],Neumat_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # Cycle_networks = list(microglia_Cycle_networks, ASTFB_Cycle_networks, Neumat_Cycle_networks)
#' # Overlap_net.all <- Overlap.net(Cycle_networks, overlap.num = 1, type = "least")
#' # Overlap_net.1 <- Overlap.net(Cycle_networks, overlap.num = 1, type = "equal")
#'
#' @references
#' Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.,22(1):578.
Overlap.net <- function(net,
                           overlap.num = 1,
                           type = c("equal", "least")){

  if(class(net)!="list") {
    stop("Please check your input network! The input network should be list object! \n")
  }

  net.transform <- unlist(lapply(seq(net), function(i) {net[[i]]=matrix(net[[i]],ncol=2);paste(net[[i]][, 1], net[[i]][, 2], sep = " & ")}))
  net.table <- table(net.transform)

  if(type == "least"){
    net.overlapped <- which(net.table >= overlap.num)
  } else if(type == "equal"){
    net.overlapped <- which(net.table == overlap.num)
  }

  overlapped <- names(net.overlapped)
  overlapped <- str_split_fixed(overlapped, " & ", 2)
  return(overlapped)
}

#' Identifying the overlap between multiple lists of hubs.
#'
#' @title Overlap.hub
#' @param hub List object, the list of hubs.
#' @param overlap.num The minimum number of hubs existing in multiple lists of hubs
#' @param type The overlapped hubs in overlap.num hub lists ("equal") or at least overlap.num hub lists ("least").
#' @export
#' @return A vector: The overlapped hubs.
#'
#' @examples
#' # NOT RUN
#' # load('ASD_exp_3cell_types.rda')
#' # microglia_Cycle_networks <- Cycle_network(Microglia_ASD_lncRNAs_data[1:100,1:100], Microglia_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # ASTFB_Cycle_networks <- Cycle_network(ASTFB_ASD_lncRNAs_data[1:100,1:100],ASTFB_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # Neumat_Cycle_networks <- Cycle_network(Neumat_ASD_lncRNAs_data[1:100,1:100],Neumat_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # Cycle_networks = list(microglia_Cycle_networks, ASTFB_Cycle_networks, Neumat_Cycle_networks)
#' # Cycle_hubs_lncRNAs <- hub_discovery(Cycle_networks)
#' # Overlap_hub.2 <- Overlap.hub(Cycle_hubs_lncRNAs, overlap.num = 2, type = "least")
#'
#' @references
#' Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.,22(1):578.
Overlap.hub <- function(hub,
                        overlap.num = 1,
                        type = c("equal", "least")){

  if(class(hub)!="list") {
    stop("Please check your input hub! The input hub should be list object! \n")
  }

  hub.transform <- unlist(hub)
  hub.table <- table(hub.transform)

  if(type == "least"){
    hub.overlapped <- which(hub.table >= overlap.num)
  } else if(type == "equal"){
    hub.overlapped <- which(hub.table == overlap.num)
  }

  overlapped <- names(hub.overlapped)
  return(overlapped)
}


#' Generate random networks in parallel and statistical significance p value of topological characteristics (path length and density) in biological networks.
#'
#' @param obser_path Observed characteristics path length of a biological network.
#' @param obser_density Observed density of a biological network.
#' @param nodes.num the number of nodes
#' @param edges.num the number of edges
#' @param perm The number of permutations for generating random networks.
#' @param directed  a logical object, TRUE or FALSE means directed or undirected network, respectively.
#' @param num.cores number of parallel cores (2 in default)
#' @import foreach
#' @importFrom foreach %dopar%
#' @importFrom pracma eps
#' @import parallel
#' @import doParallel
#' @import igraph
#' @importFrom stats pnorm
#' @importFrom stats sd
#'
#' @examples
#' #' # NOT RUN
#' # load('ASD_exp_3cell_types.rda')
#' # Cycle_networks <- Cycle_network(Neumat_ASD_lncRNAs_data[1:100,1:100], Neumat_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # g1 <- Cycle_networks
#' # g1 <- make_graph(c(t(g1[, 1:2])),directed = F)
#' # m <- gsize(g1)
#' # n <- vcount(g1)
#' # l <- mean_distance(g1)
#' # d <- edge_density(g1)
#' # pvalue_path_density <- Random_net_parallel(l, d, n, m, perm = 100, directed = FALSE, num.cores = 4)
#'
#' @return A List: The significance p-value and -log10(p-value) of topological characteristics (path length and density) of a biological network.
#' @export
#' @references Csardi G, Nepusz T. The igraph software package for complex network research.
#' InterJournal, Complex Systems. 2006, 1695.
Random_net_parallel <-
  function(obser_path,
           obser_density,
           nodes.num,
           edges.num,
           perm = 100,
           directed = FALSE,
           num.cores = 2) {
    set.seed(123)
    cores <- makeCluster(num.cores)
    registerDoParallel(cores)

    res <- foreach(i = seq(perm), .packages = "igraph") %dopar% {
      g <- sample_pa(n = nodes.num,
                     m = edges.num,
                     directed = directed)
      g <-
        delete_edges(g, sample(1:gsize(g), size = gsize(g) - edges.num))
      tmp_path <- mean_distance(g)
      tmp_density <- edge_density(g)
      return(c(tmp_path, tmp_density))
    }

    stopCluster(cores)

    res <- do.call(rbind, res)
    pvalue_density_raw <-
      pnorm(obser_density, mean(res[, 2]), sd(res[, 2]), lower.tail = FALSE)
    pvalue_density <-
      -log10(pnorm(obser_density, mean(res[, 2]), sd(res[, 2]), lower.tail = FALSE) + eps(1))

    random_res <- list(
      res <- res,
      pvalue_density_raw <- pvalue_density_raw,
      pvalue_density <- pvalue_density
    )

    return(random_res)
  }


#' Calculating similarity matrix between two list of cell type-specific lncRNA-mRNA networks
#'
#' @param net1 List object, the first list of network
#' @param net2 List object, the second list of network
#' @param directed Logical value, network directed (TRUE) or undirected (FALSE)
#' @importFrom igraph %s%
#' @import igraph
#'
#' @examples
#' # NOT RUN
#' # load('ASD_exp_3cell_types.rda')
#' # microglia_Cycle_networks <- Cycle_network(Microglia_ASD_lncRNAs_data[1:100,1:100], Microglia_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # ASTFB_Cycle_networks <- Cycle_network(ASTFB_ASD_lncRNAs_data[1:100,1:100],ASTFB_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # Neumat_Cycle_networks <- Cycle_network(Neumat_ASD_lncRNAs_data[1:100,1:100],Neumat_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # Cycle_networks = list(microglia_Cycle_networks, ASTFB_Cycle_networks, Neumat_Cycle_networks)
#' # Cycle_network_Sim <- Sim.network(Cycle_networks, Cycle_networks, directed = TRUE)
#'
#' @return A Matrix: a similarity matrix between two list of networks
#' @export
#'
#' @references
#' Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.,22(1):578.
Sim.network <- function(net1, net2, directed = TRUE) {
  if (class(net1) != "list" | class(net2) != "list") {
    stop("Please check your input network! The input network should be list object! \n")
  }

  m <- length(net1)
  n <- length(net2)
  Sim <- matrix(NA, m, n)
  for (i in seq(m)) {
    for (j in seq(n)) {
      net1_graph_interin <-
        make_graph(c(t(net1[[i]][, 1:2])), directed = directed)
      net2_graph_interin <-
        make_graph(c(t(net2[[j]][, 1:2])), directed = directed)
      overlap_interin <-
        nrow(as_data_frame(net1_graph_interin %s% net2_graph_interin))
      Sim[i, j] <-
        overlap_interin / min(nrow(net1[[i]]), nrow(net2[[j]]))
    }
  }

  return(Sim)
}


#' Calculating similarity matrix between two list of hubs
#'
#' @param hub1 List object, the first list of hub
#' @param hub2 List object, the second list of hub
#'
#' @return Sim: a similarity matrix between two list of hubs
#' @export
#'
#' @examples
#' # NOT RUN
#' # load('ASD_exp_3cell_types.rda')
#' # microglia_Cycle_networks <- Cycle_network(Microglia_ASD_lncRNAs_data[1:100,1:100], Microglia_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # ASTFB_Cycle_networks <- Cycle_network(ASTFB_ASD_lncRNAs_data[1:100,1:100],ASTFB_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # Neumat_Cycle_networks <- Cycle_network(Neumat_ASD_lncRNAs_data[1:100,1:100],Neumat_ASD_mRNAs_data[1:100,1:100], boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2, dev = TRUE, iteration = TRUE, cell_id = NULL, maxiter = 20)
#' # Cycle_networks = list(microglia_Cycle_networks, ASTFB_Cycle_networks, Neumat_Cycle_networks)
#' # Cycle_hubs_lncRNAs <- hub_discovery(Cycle_networks)
#' # Cycle_hub_Sim <- Sim.hub(Cycle_hubs_lncRNAs, Cycle_hubs_lncRNAs)
#'
#' @references
#' Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.,22(1):578.
Sim.hub <- function(hub1, hub2) {
  if (class(hub1) != "list" | class(hub2) != "list") {
    stop("Please check your input hub! The input hub should be list object! \n")
  }

  m <- length(hub1)
  n <- length(hub2)
  Sim <- matrix(NA, m, n)
  for (i in seq(m)) {
    for (j in seq(n)) {
      overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
      Sim[i, j] <-
        overlap_interin / min(length(hub1[[i]]), length(hub2[[j]]))
    }
  }
  return(Sim)
}
