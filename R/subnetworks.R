

#' Detect feedback loops in a network
#' 
#' Detects feedback loops in a network.
#' 
#' Note: this has not yet been tested and may not behave as intended.
#'
#' @param mat A matrix of edge weights
#' @param maxnodes,minnodes Maximum and minimum number of nodes 
#' permitted to be involved in the loop for it to be included in the output.
#'
#' @md
#' @return A list with:
#' 
#' * A list of character vectors, each being a loopderived from the network
#' * A \code{data.frame} with one loop per row, describing
#'     * a character description of the involved nodes and 
#'     the sign of the edges between them
#'     * The product of all edge weights, 
#'     indicating whether it is a positive or negative feedback loop
#'     * The number of nodes involved
#' * A list of matrices, each of the same size as \code{mat} and
#' retaining only the edges involved in the loop
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' 
#' 
#' 
findloops <- function(mat, maxnodes, minnodes=2){
  nodenames <- colnames(mat)
  is_symmetric <- all(mat==t(mat))

  find.cycles <- function(graph, k){
    ring <- make_ring(k, TRUE)
    subgraph_isomorphisms(ring, graph)
  }
  
  netloops <- list()
  loopdepictions <- list()

  g <- graph_from_adjacency_matrix(abs(mat), weighted=T)
  l <- unlist(lapply(minnodes:maxnodes, find.cycles, graph=g), recursive=FALSE)
  
  #extract the vertices in each cycle
  looplist <- Filter(Negate(is.null), lapply(l, function(e) {
    if (length(e) > 1L) {
      nm <- names(e)
      c(nm, nm[1L])
    }
  }))
  
  # Ensure the first node is always the same node; then remove identical nodes
  for(j in seq_along(looplist)){
    preferredinitial <- nodenames[which(nodenames %in% looplist[[j]])[1]]
    newinitial <- which(looplist[[j]] == preferredinitial)[1]
    if(newinitial != 1){
      looplist[[j]] <-
        c(looplist[[j]][c(newinitial:(length(looplist[[j]])-1),1:(newinitial-1))],
          looplist[[j]][newinitial])
    }
  }
  looplist <- unique(looplist)
  
  # Remove identical loops in opposite direction if undirected network
  if(is_symmetric){
    loopvecs <- sapply(looplist,paste,collapse=" ")
    revloopvecs <- sapply(looplist,\(x){paste(rev(x),collapse=" ")})
    matches <- cbind(seq_along(loopvecs),match(loopvecs,revloopvecs))
    matches <- apply(matches,1,sort,na.last=T)
    matches <- unique(matches)
    looplist <- looplist[matches[,1,drop=T]]
  }
    
  # Compute edge product per loop
  # Also produce submatrices per loop and textual descriptions
  loopprods <- numeric(length(looplist))
  looptexts <- looplist
  loopmats <- list()
  for(j in seq_along(looplist)){
    key <- cbind(looplist[[j]],lead(looplist[[j]]))[-length(looplist[[j]]),]
    powers <- mat[key]
    loopprods[j] <- prod(powers)
    
    looptexts[[j]] <- paste(looptexts[[j]],c(ifelse(sign(powers)==1,"+>","->"),"")) %>%
      paste(collapse=" ") %>% trimws()
    looptexts[[j]] <- paste0(ifelse(sign(prod(powers))==1,"(+) ","(-) "),looptexts[[j]])
    
    currloopmat <- mat
    currloopmat[] <- 0
    currloopmat[key] <- mat[key]
    loopmats[[j]] <- currloopmat
  }
  
  if(length(looplist)>0){
    looplengths <- sapply(looplist,length)-1
  }else{
    looplengths<-numeric(0)
  }
  
  loopdf <- data.frame(looptext=unlist(looptexts),
                       loopprod=loopprods,
                       looplength=looplengths)
  
  out<-list(loops=looplist,stats=loopdf,mats=loopmats)
  return(out)
}


#' Get all connected nodes in a network
#'
#' @param mat An adjacency matrix.
#' @param net A \code{qgraph} object. 
#' One of \code{mat} or \code{net} must be given.
#' @param node A node name. If given, the function only returns nodes that are
#' connected to this node, directly or indirectly.
#' @param ignore.selfloops Should self-loops be ignored when determining
#' if a node is unconnected?
#' @param ignore.negative Should negative edges be ignored when determining
#' if a node is unconnected?
#' @param components If \code{"largest"}, this function only returns nodes 
#' that are part of the largest component of the network. If \code{"all"},
#' no such subsetting is performed.
#'
#' @return A vector of node names that are connected to other nodes according to
#' the rules provided in this function's arguments.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' 
#' 
#' 
getConnectedNodes <- function(mat=NULL, net=NULL, node=NULL, 
                              ignore.selfloops=TRUE,
                              ignore.negative=FALSE,
                              components=c("all", "largest")){
  
  # Get both network and matrix form of the network
  if(is.null(net)){
    stopifnot(!is.null(mat))
    net <- qgraph::qgraph(mat, DoNotPlot=T, labels=colnames(mat))
  }else if(is.null(mat)){
    stopifnot(!is.null(net))
    mat <- getWmat(net)
  }
  
  # Extract node names from largest component if desired, else get all node names
  components <- match.arg(components)
  if(components == "largest"){
    coms <- components(igraph::graph_from_adjacency_matrix(mat))
    largest <- which.max(coms$csize)
    key <- colnames(mat)[coms$membership == largest]
  }else{
    key <- colnames(mat)
  }
  
  # Detect all connected nodes
  if(!is.null(node)){
    if(ignore.negative){
      negfun <- function(x){ x[x < 0] <- 0; x }
    }else{
      negfun <- abs
    }
    paths <- centrality(mat,posfun=negfun)$ShortestPathLengths
    reachfrom <- rownames(paths)[!is.infinite(paths[node,])]
    reachto <- colnames(paths)[!is.infinite(paths[,node])]
    key2 <- c(node,reachfrom,reachto) |> unique()
  }else{
    if(ignore.selfloops){
      diag(mat) <- 0
    }
    if(ignore.negative){
      mat[mat<0] <- 0
    }else{
      mat <- abs(mat)
    }
    key2 <- unique(c(rownames(mat)[rowSums(mat)!=0],
                     colnames(mat)[colSums(mat)!=0]))
  }
  
  # Return connected node names
  key <- intersect(key,key2)
  out <- net$Arguments$labels[net$Arguments$labels %in% key]
  return(out)
}

#' @rdname getConnectedNodes
#' @export
#' 
pruneUnconnectedNodes <- function(mat=NULL, node=NULL, 
                                  ignore.selfloops=TRUE,
                                  ignore.negative=FALSE,
                                  components=c("all", "largest")){
  components <- match.arg(components)
  connodes <- getConnectedNodes(mat=mat,
                                node=node,
                                ignore.selfloops=ignore.selfloops,
                                ignore.negative=ignore.negative,
                                components=components)
  return(mat[connodes,connodes])
}

#' Get all edges leading to/from a node
#'
#' @param mat An adjacency matrix
#' @param from,to Node names between which the shortest path must be computed.
#' At least one of these must be given. 
#' If only \code{from} is given, then this returns 
#' the shortest paths from this node to all other nodes, 
#' and if only \code{to} is given, then this returns 
#' the shortest paths to this node from all other nodes.
#' 
#' @return The same adjacency matrix, but the only non-null values are those edges
#' on the shortest path from the \code{from} node and/or to the \code{to} node.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' 
#' 
#' 
getNodePaths <- function(mat, from=NULL, to=NULL){
  
  # Compute shortest paths
  cents <- centrality(mat, all.shortest.paths=TRUE)
  
  # Extract relevant paths from centralities object
  if(!is.null(from) & is.null(to)){ 
    currpaths <- cents$ShortestPaths[from,] 
  }else if(is.null(from) & !is.null(to)){ 
    currpaths <- cents$ShortestPaths[,to] 
  }else if(!is.null(from) & !is.null(to)){ 
    currpaths <- cents$ShortestPaths[from,to] 
  }else{
    stop("Args from & to are both missing; supply at least one")
  }
  
  # takes a vector of IDs and turns it into a two-column matrix
  # indicating the path from each node to the next in an adjacency matrix
  ids2inds <- function(x){ 
    matrix(c(x[-length(x)],x[-1]), ncol=2, nrow=length(x)-1)
  }
  
  # Create a mask for the adjacency matrix and set all valid edges to 1
  maskmat <- matrix(0, nrow=nrow(mat), ncol=ncol(mat))
  for(i in seq_along(currpaths)){
    if(length(currpaths[i][[1]])>0){
      maskmat[ids2inds(currpaths[[i]][[1]])]<-1
    }
  }
  
  return(mat*maskmat)
}

#' Make an adjacency matrix symmetric
#' 
#' This makes a adjacency matrix symmetric by averaging each value with its transpose.
#'
#' @param mat An adjacency matrix.
#' @param rule When two nodes are not bidirectionally connected, what should be done?
#' \code{"AND"} (default) deletes the edge, 
#' while \code{"OR"} keeps the halved edge weight.
#' 
#' @details This follows the implementation of the IsingFit package.
#'
#' @return A symmetric adjacency matrix.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' 
#' 
symmetrizeNetwork <- function(mat, rule=c("AND","OR")){
  rule <- match.arg(rule)
  symmmat <- (mat + t(mat))/2
  if(rule == "AND"){
    nonzeromask <- mat!=0 & t(mat)!=0
    symmmat[!nonzeromask] <- 0
  }
  return(symmmat)
}

