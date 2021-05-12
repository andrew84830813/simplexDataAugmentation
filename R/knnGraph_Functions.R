#' Computes a K-Farthest / Nearest Neighbor Graph
#'
#' Compute a K-Nearest/Farthest Neighbor (KFN / KNN) graph from a weighted adjacency matrix
#'
#' @param adj_mat pxp Adjacency Matrix
#' @param K K farthest neighbors
#' @param plot_TrueFalse should results be plotted
#' @param sim_ returns a K nearest neighbor graph
#'
#' @return A list containing:\tabular{ll}{
#'    \code{Graph} \tab the igraph graph object  \cr
#'    \tab \cr
#'    \code{adj} \tab The KFN or KNN weighted adjacency matrix \cr
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## compute pairwise distance
#' comp = sampleAddLogisticT(n1 = 30)[,-1]
#' dd = dist(compositions::clr(comp))
#' ## Compute 5-KNN Graph
#'  g = knn_graph(as.matrix(dd),K = 5,sim_ = T)
#' }
knn_graph <-
  function(adj_mat,K,plot_TrueFalse=FALSE,sim_=T){
    ## parms
    plt = plot_TrueFalse
    adj = adj_mat
    adj2 = adj_mat
    k = K

    #Similarity of Dissim
    if(sim_ == TRUE){
      diag(adj)=Inf# for smallest
    }else{
      diag(adj)=0
      adj = -adj
    }

    #Choose K
    k_nn = t(seq(1,k,by=1))

    #Replace with Adj Mat with Ranks
    for (r in 1:nrow(adj)){
      tt= as.matrix(rank(adj[r,],ties.method = "random"))
      tt= ifelse(tt %in% k_nn,1,0)
      tt=tt*adj2[r,]
      adj2[r,]=tt
    }

    diag(adj2)=0

    #Make Matric Symetric
    adj =  knnADJtoSYM(adj2)

    #Create Graph
    w.graph = igraph::graph.adjacency(adj,mode="undirected",weighted = TRUE)
    w.graph.simplified = igraph::simplify(w.graph, remove.loops = TRUE,
                                          edge.attr.comb = igraph::igraph_opt("edge.attr.comb"))
    #Plot if TRUE
    if (plt == TRUE){
      l = igraph::layout.fruchterman.reingold(w.graph.simplified)

      plot(w.graph.simplified,
           layout = l,
           main = "Test",
           vertex.size = 1,
           vertex.label = NA,
           edge.width =1.25,
           edge.color = "grey"
      )
    }


    list(Graph = w.graph.simplified,
         AdjMatrix = adj
    )
  }

#' Make Symetric
#'
#' Converts a KNN/KFN adjaceny matrix to a symetric matrix
#'
#' @param knnADJ_MAT knn / kfn dj. matrix
#'
#' @return a symetric knn/kfn weighted adj matrix
#'
knnADJtoSYM <-
  function(knnADJ_MAT){
    ut = NULL
    adj = knnADJ_MAT
    adjT = knnADJ_MAT
    adjT  =t(adj)
    lt = adj[lower.tri(adj)]
    uT = adjT[lower.tri(adjT)]
    df = data.frame(lt = lt,ut = uT)
    message("dplyr Start")
    df = dplyr::mutate(df,s = dplyr::if_else(lt==ut,lt,pmax(ut,lt)))
    message("dplyr Finished")

    adj[lower.tri(adj)]=df$s
    adj = t(adj)
    adj[lower.tri(adj)]=df$s

    return(adj)
  }
