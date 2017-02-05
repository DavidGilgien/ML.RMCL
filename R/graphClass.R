#' A Graph Object
#'
#' This function allows you to create graph objects from an edge list.
#' @param edgeList \code{data.table} with a column i and a cloumn j of edges endpoints
#' @param nodeList vector of node names. if \code{NULL}, it is derived from the edge list
#' @param nodeMap table of nodes mapped. Can be used when the graph is the result of coarsening
#' @param clusters vector of size \code{length(nodeList)} of truth cluster belonging.
#' @param directed boolean value describing if the edges are directed.
#' @keywords graphs object
#' @return A graph object 
#' @import data.table
#' @export
#' @examples
#' g = ELgraph(data.table("i"=c(2,5,6,3,2), "j" =c(4,4,3,1,5)))
ELgraph <- function(edgeList, nodeList = NULL, nodeMap = NULL, clusters = NULL, directed = TRUE)
{
  
  if(is.null(nodeList)){
    nodeList = sort(unique(c(edgeList[,i], edgeList[,j])))
  }
  if(is.null(edgeList$x)){
    edgeList$x = 1
  }
  if(!all(nodeList != 0)){
    new = max(nodeList)+1
    edgeList[i == 0, i:= new]
    edgeList[j == 0, j:= new]
    nodeList = c(nodeList[-1], new)
  }
  me <- list(
    edgeList = edgeList,
    nodeList = nodeList,
    nodeMap = nodeMap,
    clusters = clusters,
    directed = directed
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"ELgraph")
  return(me)
}


#' Extract the Adjancency Sparse Matrix
#'
#' This function allows you to get the adjacency matrix from a graph object.
#' @param graph The graph from which to obtain the matrix
#' @keywords sparse matrix adjacency
#' @return A matrix of sparse format from package Matrix with entries corresponding to the edge list of the given graph
#' @export
#' @import Matrix
#' @examples
#' M = getSparseMat(g)
getSparseMat <- function(graph)
{
  M = sparseMatrix(i = match(graph$edgeList$i, graph$nodeList),
                   j = match(graph$edgeList$j, graph$nodeList),
                   x = graph$edgeList$x,
                   dims = rep(length(graph$nodeList), 2),
                   dimnames = list(graph$nodeList, graph$nodeList))
  if(!graph$directed){
    M = M + t(M)
    diag(M) = diag(M)/2
  }
  return(M)
}