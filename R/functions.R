
hamming_loop <- function(X) {
  d <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  for ( i in 1:(nrow(X) - 1) ) {
    for ( j in (i + 1):nrow(X) ) {
      d[i,j] <- d[j,i] <- sum(X[i,] != X[j,]) / length(X[i,])
    }
  }
}
hamming_binary <- function(X) {
  D <- t(1 - X) %*% X
  D + t(D)
}

doEdgesCluster<-function(X,nn=30,do.jaccard=TRUE,method="Louvain") {
  nearest<-nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
  print("Found nearest neighbors")
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1] #Convert to a similarity score
  nearest$nn.sim <- 1*(nearest$nn.dists >= 0 )
  edges <- melt(t(nearest$nn.idx))
  colnames(edges) <- c("B", "A", "C")
  edges <- edges[,c("A","B","C")]
  edges$B <- edges$C
  edges$C<-1
  
  #Remove repetitions
  edges <- unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
  if (do.jaccard){
    NN = nearest$nn.idx
    jaccard_dist <- apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
    edges$C <- jaccard_dist
    edges <- subset(edges, C != 0) ## skip no-overlapping
    edges$C <- edges$C/max(edges$C) ## convert to ratio
  }
  Adj <- matrix(0, nrow=nrow(X), ncol=nrow(X))
  rownames(Adj) <- rownames(X)
  colnames(Adj) < rownames(X)
  Adj[cbind(edges$A,edges$B)] <- edges$C ## up
  Adj[cbind(edges$B,edges$A)] <- edges$C ## down
  g<-graph.adjacency(Adj, mode = "undirected", weighted=TRUE)
  if (method=="Louvain") graph.out <- cluster_louvain(g)
  if (method=="Infomap") graph.out <- cluster_infomap(g)
  return(graph.out)
}
doSNNCluster<-function(X){
  # Run SNN graph
  grps <- buildSNNGraph(t(X),rand.seed = 1000)
  # extract memebership
  clusters <- cluster_fast_greedy(grps)
  memb <-  membership(clusters)
  return(memb)
}
plotTSneClustering <- function(X, tsne, methodname) {
  X1 <- tsne$Y[,1]
  X2 <- tsne$Y[,2]
  e <- doEdgesCluster(X,
                      nn = 30,
                      do.jaccard = TRUE,
                      method = methodname)
  mem <- e$membership
  d <-
    data.frame(
      'X1' = X1,
      'X2' = X2,
      'membership' = as.factor(e$membership)
    )
  md <- melt(d, id.vars = c("X1", "X2"))
  pl <- ggscatter(
    data = md,
    x = 'X1',
    y = 'X2',
    color = 'value',
    font.label = 8
  ) + theme_bw() + ggtitle(paste( methodname, ' Clustering'))
  return(list('plot' = pl, 'mem' = mem))
}