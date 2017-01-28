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