#' Graph Creation From File
#'
#' This function creates a graph object from a file containing the edges
#' @param path string: the path to the file containing the edge list.
#' @param nskip integer: number of lines to skip before starting to read the edges.
#' @param directed boolean: value describing if the edges are directed.
#' @param norm boolean: if \code{TRUE}, translate the vertices names so that they are from 1 to n.
#' @param inverse.direction boolean: if \code{FALSE}, labels the first column \code{i} and the second \code{j}, otherwise the inverse.
#' this option is used when creating a graph for the ML-RMCL function, because it understand the edges going from \code{j} to \code{i}.
#' @keywords graph initialization file edge list
#' @return A graph object of class \code{\link{ELgraph}} from the edge list found at \code{path}
#' @export
#' @examples
#' g = createGraph("graphs/myEdges.txt")
createGraph = function(path, nskip = 0, directed = TRUE, norm = TRUE, inverse.direction = TRUE){
  el = data.table(read.table(path, skip = nskip))
  if(inverse.direction) colnames(el) = c("j", "i")
  else colnames(el) = c("i", "j")
  g = ELgraph(el, directed = directed)
  if(norm) g = normalizeNames(g)
  return(g)
}

normalizeNames = function(g){
  el = g$edgeList
  el[, i:=match(el$i, g$nodeList)]
  el[, j:=match(el$j, g$nodeList)]
  return(ELgraph(el, directed = g$directed))
}


subsetGraph = function(g, keepNodes){
  el = g$edgeList
  el = el[i %in% keepNodes, ]
  if(any(is.na(match(el$j,  keepNodes)))) stop("The keepNodes do not form an independant component")
  return(ELgraph(el, directed = g$directed))
}




