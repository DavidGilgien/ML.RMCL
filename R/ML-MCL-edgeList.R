
HWM = function(g){
  m1 = c()
  m2 = c()
  r = sample(g$nodeList, length(g$nodeList))
  el = g$edgeList[i!=j]
  setkey(el, j)
  flag = rep(TRUE, max(g$nodeList))
  k = 1
  for(node in r){
    if(flag[node]){
      subset = el[.(node),]
      if(!is.na(subset[1,x])){
        subset = subset[flag[i]]
        if(nrow(subset)!=0){
          subset = subset[x == max(x), list(i, j)]
          match = subset[sample(nrow(subset), 1),c(i,j)]
          m1[k] = match[1]
          m2[k] = match[2]
          flag[match] = FALSE
          k = k+1
        }
      }
    }
  }
  return(data.table("m1" = m1, "m2" = m2))
}

## with weights
coarse = function(g, weights = TRUE){
  el = copy(g$edgeList)
  map = HWM(g)
  m = map[match(el$i, map$m1), m2]  #inverse m1/m2
  el[!is.na(m), i:=na.omit(m)]
  m = map[match(el$j, map$m1), m2]  #inverse m1/m2
  el[!is.na(m), j:=na.omit(m)]
  if(weights) el = el[, x:=sum(x), by=list(i, j)]
  el = unique(el)
  return(ELgraph(edgeList = el, nodeMap = map, directed = g$directed))
}


projectFlow = function(net, M_c){
  # Is called if the matrix was not sparse enough and was stored in a different format
  if(!(class(M_c) == "dgCMatrix")){
    M_c = as(M_c, "dgCMatrix")
  }
  el = data.table(summary(M_c))
  el$i = net$nodeList[el$i]
  el$j = net$nodeList[el$j]
  new_el = copy(el)
  new_nodes = net$nodeMap$m1  #newNodes m1
  m = net$nodeMap[match(el$j, net$nodeMap$m2), m1] #inverse m1/m2
  new_el[!is.na(m), j:=na.omit(m)]
  el = unique(rbind(el, new_el))
  M_r = getSparseMat(ELgraph(el, nodeList = sort(c(net$nodeList, new_nodes))))
  return(M_r)
}

MCL_iter = function(M, M_g, reg, r, eps){
  
  #Expansion
  if(reg){
    M_exp = M%*%M_g
  }else{
    M_exp = M%*%M
  }
  
  #Inflation
  M_infl = M_exp ^ r
  M = M_infl %*% Diagonal(x = 1/Matrix::colSums(M_infl))
  
  
  #Pruning
  if(!(class(M) == "dgCMatrix")) M = as(M, "dgCMatrix")
  el = data.table(summary(M))[, count:= .N, by = j]
  el = el[x > eps/count, list(i,j,x)]
  M = sparseMatrix(i = el$i,
                   j = el$j,
                   x = el$x,
                   dims = dim(M), 
                   dimnames = dimnames(M))
  M = round(M, 12) # Some values were not accurate due to precision points
  M = M %*% Diagonal(x = 1/Matrix::colSums(M))
  return(M)
}

interpret = function(M){
  dat = data.table(summary(M))[, m:= max(x) ,by = j]
  dat = dat[m == x, min(i), by = j]
  return(dat$V1)
}


weighTrans = function(graph){
  el = copy(graph$edgeList)
  outdeg = el[, sum(x), by = i]
  indeg = el[, sum(x), by = j]
  deg = merge(outdeg, indeg, by.x = "i", by.y = "j", all=TRUE)
  deg[ , deg :=rowSums(.SD, na.rm = TRUE), .SDcols = c("V1.x", "V1.y")]
  deg = deg[, list(i, deg)]
  if(!graph$directed) deg[, deg:= deg*2]
  el = merge(el, deg, by.x = "j", by.y = "i")
  el = merge(el, deg, by.x = "i", by.y = "i")
  el[, x := x/deg.x + x/deg.y]
  g2 = copy(graph)
  g2$edgeList = el[,list(i,j,x)]
  return(g2)
}

#' Multi-Level Regularized Markov Clustering
#'
#' This function implements the multi-level regularized markov clustering. It uses flow simulation of a graph to find clusters.
#' @param g graph object of class \code{\link{ELgraph}} on which to perform the clustering
#' @param r inflation parameter. A higher value yields more refined clusters.
#' @param reg boolean value. When \code{TRUE}, the MCL iterations in phase 2 are regularized.
#' @param iter.max integer: maximum number of MCL iterations during phase 3 before the matrix is interpreted as a clustering. The function gives a warning when reached.
#' @param n.threshold integer: number of nodes under which phase 1 stops the coarsening.
#' @param last.reg boolean value. When \code{TRUE}, the MCL iterations in phase 3 are regularized.
#' @param eps integer: epsilon threshold for the pruning operation. During pruning, values under \code{eps/m} are set to 0, where \code{m} is the number of non-zero entries in the column
#' @param save.p1 string: a name for the storing a the graph hierarchy after the coarsening. A file will be saved at "./phase1/\code{save.p1}.rds". The algorithm can then be used with this file as a strating point by giving the name as the \code{load.p1} parameter.
#' @param load.p1 string: a name for loading the graph hierarchy. When not \code{NULL}, \code{g} will not be used but the function will take the file "./phase1/\code{load.p1}.rds" as input. Can be used when phase 1 is time consuming and we want different tests with the same coarsening.
#' @param iter.curt integer: number of iteration of MCL during phase 2.
#' @param coarse.weight boolean value. When \code{TRUE}, the coarsening phase keeps the edges' weights as indications of the number of merged edges.
#' @param init.weight.trans boolean value. When \code{TRUE}, the algorithm performs a weight transformation step at the very beginning, before the coarsening.
#' @param curt.weight.trans boolean value. When \code{TRUE}, the algorithm performs a weight transformation step for each flow matrix initialisation during phase 2.
#' @param timer boolean value. When \code{TRUE}, the algorithm returns a vector of phases duration in addition of the clustering.
#' @keywords multilevel multi-level regularizd markov clustering flow simulation
#' @return A vector of size \code{n} describing the associated cluster of each node.
#' When \code{timer=TRUE}, a list of two elements: \code{res} the clustering vector and 
#' \code{time} a vector of size 3 with the duration of each phase in seconds.
#' @export
#' @examples
#' g = generateEdgeList(n = 1000)
#' clustering = ML_RMCL(g)
ML_RMCL = function(g, r = 2, reg = TRUE, iter.max = 200, n.threshold = 200,
                  last.reg = TRUE, eps = 0.001, save.p1 = NULL, load.p1 = NULL, 
                  iter.curt = 4, coarse.weight = TRUE, 
                  init.weight.trans = FALSE, curt.weight.trans = FALSE,
                  timer = FALSE){
  
  
  # Phase 1: coarsening
  p0 = as.numeric(Sys.time())
  if(is.null(load.p1)){
    if(init.weight.trans) g = weighTrans(g)
    graphs = list()
    graphs[[1]] = g
    i=2
    ratio = 0
    while(length(g$nodeList) > n.threshold & ratio < 0.9){
      g = coarse(g, weights = coarse.weight)
      graphs[[i]] = g
      ratio = length(g$nodeList)/length(graphs[[i-1]]$nodeList)
      i = i+1
    }
    if(!is.null(save.p1)) saveRDS(graphs, paste0("phase1/", save.p1, ".rds"))
  }else{
    graphs = readRDS(paste0("phase1/", load.p1, ".rds"))
    sizes = sapply(graphs, function(x) length(x$nodeList))
    if(n.threshold > min(sizes)) graphs = graphs[1:which(sizes < n.threshold)[1]]
  }
  p1 = as.numeric(Sys.time())
  # Phase 2: curtailed MCL:
  if(length(graphs) != 1){
    for(i in seq(length(graphs), 2, -1)){
      if(curt.weight.trans) M_g = getSparseMat(weighTrans(graphs[[i]]))
      else M_g = getSparseMat(graphs[[i]])
      diag(M_g) = 1
      M_g = M_g %*% Diagonal(x = 1/Matrix::colSums(M_g))
      if(i == length(graphs)) M = M_g
      for(j in seq(1,iter.curt)){
        M = MCL_iter(M, M_g, reg, r, eps)
      }
      M = projectFlow(graphs[[i]], M)
    }
  }
  # Phase 3: Run MCL on the original graph
  p2 = as.numeric(Sys.time())
  M_g = getSparseMat(graphs[[1]])
  diag(M_g) = 1
  M_g = M_g %*% Diagonal(x = 1/Matrix::colSums(M_g))
  if(length(graphs) == 1) M = M_g
  M_last = M
  i = 0
  repeat{
    M = MCL_iter(M, M_g, last.reg, r, eps)
    i = i+1
    if(i > iter.max){
      warning("Max iterations reached, tried to interpret the flow matrix anyways")
      if(timer) return(list("res" = interpret(M), "time" = c(p1-p0, p2-p1, as.numeric(Sys.time())-p2)))
      else return(interpret(M))
    }
    if(identical(M, M_last)){
      break
    } 
    M_last = M
  }
  p3 = as.numeric(Sys.time())
  if(timer) return(list("res" = interpret(M), "time" = c(p1-p0, p2-p1, p3-p2)))
  else return(interpret(M))
}





