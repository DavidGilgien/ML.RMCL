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
  M = M_infl %*% Diagonal(x = 1/colSums(M_infl))
  
  
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
  M = M %*% Diagonal(x = 1/colSums(M))
  return(M)
}


ML_MCL = function(g, r = 2, reg = TRUE, iter.max = 200, n.threshold = 200,
                  lastReg = TRUE, eps = 0.001, save.p1 = NULL, load.p1 = NULL, 
                  iter.curt = 4, coarse.weight = TRUE, 
                  init.weight.trans = FALSE, 
                  curt.weight.trans = FALSE){
  
  
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
  graphs <<- graphs
  # Phase 2: curtailed MCL:
  if(length(graphs) != 1){
    for(i in seq(length(graphs), 2, -1)){
      if(curt.weight.trans) M_g = getSparseMat(weighTrans(graphs[[i]]))
      else M_g = getSparseMat(graphs[[i]])
      #diag(M_g) = colMeans(M_g)
      #diag(M_g)[diag(M_g) == 0] = 1
      diag(M_g) = 1
      M_g = M_g %*% Diagonal(x = 1/colSums(M_g))
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
  #diag(M_g) = colMeans(M_g)
  diag(M_g) = 1
  #diag(M_g)[diag(M_g) == 0] = 1
  M_g = M_g %*% Diagonal(x = 1/colSums(M_g))
  if(length(graphs) == 1) M = M_g
  M_last = M
  i = 0
  repeat{
    M = MCL_iter(M, M_g, lastReg, r, eps)
    i = i+1
    if(i > iter.max){
      warning("Max iterations reached, tried to interpret the flow matrix anyways")
      return(list("res" = interpret(M), "time" = c(p1-p0, p2-p1, as.numeric(Sys.time())-p2)))
    }
    if(identical(M, M_last)){
      break
    } 
    M_last = M
  }
  p3 = as.numeric(Sys.time())
  return(list("res" = interpret(M), "time" = c(p1-p0, p2-p1, p3-p2)))
}





