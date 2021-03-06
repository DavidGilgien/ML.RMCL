
#' Synthetic Graph Generation
#'
#' This function creates a synthetic graph with truth clusters properties.
#' @param k integer: the number of clusters.
#' @param aid integer: average intra-cluster degree. Poisson parameters for number of neighbors within the cluster.
#' @param aed integer: average inter-cluster degree. Poisson parameters for number of neighbors outside the cluster.
#' @param n integer: number of vertices in the graph.
#' @param prob vector of size \code{k} with probabilities that a node belongs to a cluster. If omitted, all probabilities are equal \code{1/k}.
#' @param directed boolean value describing if the edges are directed.
#' @keywords graph generation
#' @return A graph object of class \code{\link{ELgraph}} with truth clusters.
#' @export
#' @examples
#' g = generateEdgeList(n=1000)
generateEdgeList = function(k = 2, aid = 100, aed = 10, n = 1000, prob = NULL, directed = FALSE){
  z = sample(k, size = n, replace = TRUE, prob = prob)
  el_list = list()
  for (i in seq(1,n,1)) {
    j_same = sample(which(z==z[i]), size = rpois(1,aid))
    j_diff = sample(which(z!=z[i]), size = rpois(1,aed))
    if(length(j_same) + length(j_diff) > 0) el_list[[i]] = data.table("i" = i, "j" = as.numeric(c(j_same, j_diff)))
  }
  el = rbindlist(el_list)
  if(!directed){
    el = el[el[, .I[1], by = list(pmin(i, j), pmax(i, j))]$V1]
  }
  return(ELgraph(el, clusters = z, directed = directed))
}


nCut = function(g, cluster, ignore.zero = FALSE){
  el = g$edgeList
  cut = el[xor(i %in% cluster, j %in% cluster), sum(x)]
  deg = el[i %in% cluster & j %in% cluster, sum(x)] * 2
  if(deg == 0){
    warning("One cluster has degree 0")
    if(ignore.zero) return(NA)
    else return(cut)
  }
  return(cut/deg)
  
}


#' Average Normalized N-Cut
#'
#' This function computes the average normalized n-cut of a clustering.
#' @param g graph object of class \code{\link{ELgraph}}, the base graph.
#' @param clusters vector of size \code{n}. Element \code{i} indicates the cluster of node \code{v_i}.
#' @param verb boolean value: if \code{TRUE}, the function prints the n-cut of all clusters.
#' @param ignore.zero boolean value: if \code{TRUE}, clusters with internal degree 0 are ignored.
#' @keywords average normalized ncut n-cut
#' @return integer: the average normalized n-cut of a clustering.
#' @export
#' @examples
#' g = generateEdgeList(n=1000)
#' res = ML_RMCL(g)
#' avgNcut(g, res, verb = TRUE)
avgNcut = function(g, clusters, verb = FALSE, ignore.zero=FALSE ){
  ncs = c()
  for(i in unique(clusters)){
    cl = g$nodeList[clusters == i]
    nc = nCut(g, cl, ignore.zero = ignore.zero)
    if(verb) print(paste0("cluster ", i, " has a ncut of ", nc))
    ncs= c(ncs, nc)
  }
  return(mean(ncs, na.rm = TRUE))
}


#' Clustering Evaluation Against Ground Truth
#'
#' This function evaluates the quality of a clustering when a ground truth clustering is available. It computes the
#' percentage of mismatched nodes when matching the clusters.
#' @param x vector of size \code{n} describing a clustering.
#' @param y vector of size \code{n} describing a second clustering that will be compared to \code{x}.
#' @param verb boolean value: if \code{TRUE}, the function prints a line for each cluster.
#' @param perf boolean value: if \code{TRUE}, the function stops when it detects a different number of clusters and return 0.
#' @keywords clustering evaluation ground truth checking
#' @return integer: percentage of correctly matched nodes when comparing the two clustering \code{x} and \code{y}.
#' @export
#' @examples
#' g = generateEdgeList(n=1000)
#' res = ML_RMCL(g)
#' clusterEval(g$clusters, res)
clusterEval = function(x, y, verb = TRUE, perf = FALSE){
  tot = 0
  if(length(unique(x)) != length(unique(y))){
    if(verb) print(paste("Different number of clusters:", length(unique(x)), "and",length(unique(y))))
    if(perf) return(0)
    else warning("Not perfect")
  }
  for(i in unique(x)){
    tab = table(y[x == i])
    j = names(which(tab == max(tab))[1])
    f = sum(tab[names(tab) != j])
    if(verb) print(paste("Cluster", i, "matches to", j, "with **", f , "** faults"))
    tot = tot + f
  }
  if(verb) print(paste("Total missmatch:", tot))
  if(verb) print(paste("Rate of success: ", 1-(tot/length(x))))
  return(1-(tot/length(x)))
}