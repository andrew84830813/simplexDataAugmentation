
#' Simulate / Augment Simplex Bound Data from Emipircal Distribution
#'
#' @importFrom foreach %dopar%
#'
#' @param comp empirical composition
#' @param labels class labels ( can be extracted later if needed)
#' @param k K nearest neighbors for topology matching
#' @param modthres modulairty threshold for simulating using community structure
#' @param ndirections number of random directions to sample from
#' @param fct scale factor for dispersion distance matching
#' @param alFactor alpha factor dirichlet concentration parameter vetcor
#' @param ptPerSample numer of points to agumnet per sample
#'
#' @return A list containing:\tabular{ll}{
#'    \code{simData} \tab simulated/augmented data  \cr
#'    \tab \cr
#'    \code{alphaData} \tab detailed metrics on distances and aplha power coefficents used to shift points \cr
#'     \tab \cr
#'    \code{ClassLabels} \tab Class labels of simulated/aug. data matching the class of origin
#' }
#' @export
#'
#' @examples
#'\dontrun{
#'add example here
#'}
simCompData = function(comp,labels = NULL,k = round(sqrt(nrow(comp))),modthres = .4,ndirections = 1000,fct = 1.5, alFactor = 1.5,ptPerSample = 3){

    ## compute pairwise distance
    dd = parallelDist::parallelDist(compositions::clr(comp))
    ## Compyte KNN Topology and communiyt structure
    g = knn_graph(as.matrix(dd),K = k,sim_ = T)
    #Extract KNN Adj Matrix
    adjacency_matrix <- igraph::as_adjacency_matrix(graph = g$Graph,sparse = F,
                                                    attr = "weight"
    )
    ## Define Community Structure
    partition <- leiden::leiden(object = adjacency_matrix)
    cls = data.frame(membership = partition)
    testMod = igraph::modularity(g$Graph,partition)
    plot(g$Graph,vertex.color =factor(cls$membership),vertex.label = NA,vertex.size = 5,edge.width = .5)


    ## if significant use community structure
    if(testMod<modthres){
      message("Not signif.", "; modularity = ",round(testMod,3))
      cls = data.frame(membership = rep(1,nrow(comp)))
      comIDs = 1
    }else{
      message("Signif.", "; modularity = ",round(testMod,3))
      ss = table(cls$membership)

      ## combine groups with single member
      ii = which(as.numeric(ss)==1)
      if(!rlang::is_empty(ii)){
        if(length(ii)==1){
          comb = as.numeric(names(ss)[ii])
          cls$membership[cls$membership %in%comb] = max(unique(cls$membership))
        }else{
          comb = as.numeric(names(ss)[ii])
          cls$membership[cls$membership %in%comb] = max(unique(cls$membership))+1
        }
      }
    }

    comIDs = unique(cls$membership)
    comIDs = sort(comIDs)

  ## Data Frames
  shiftsPoints = data.frame()
  reqAlpha = data.frame()
  allDat = data.frame()

  message("Compute Augmentation....")
  for(cid in comIDs){

    ## select community points
    phd = data.frame(comp[cls$membership==cid,])
    phd.labels = labels[cls$membership==cid]
    est = compositions::fitDirichlet(compositions::acomp(phd))
    alpha_ = est$alpha
    simDir_vec = compositions::rDirichlet.acomp(ndirections,alpha = alFactor*alpha_)#c(mean_comp)*200)

    compMean.fixed = compositions::mean.acomp(compositions::acomp(phd))
    dspLen.fixed = data.frame()
    for(i in 1:nrow(phd)){
      ph = compositions::clo(as.numeric(phd[i,])/as.numeric(compMean.fixed))
      ph = data.frame(pointNum = i,disCentr = sqrt(sum(compositions::clr(ph)^2)))
      dspLen.fixed = rbind(dspLen.fixed,ph)
    }
    d = dspLen.fixed$disCentr
    cv = stats::sd(d) / mean(d)
    mn_d = mean(d*fct*cv); sd_d = stats::sd(d*fct*cv)

    for(r in 1:ptPerSample){

      ad = foreach::foreach(i = 1:nrow(phd),.combine = rbind)%dopar%{

        ## seelct point
        pt = ((as.matrix(phd[i,])))
        ## randomly sample a direction
        rx = sample(x = 1:nrow(simDir_vec),size = 1)
        randDir = simDir_vec[rx,]
        ## define pertubation from point to new direction
        pertb = compositions::clo(as.numeric(pt) / as.numeric(randDir)) ## from point -> rand dir
        ## sample shift distance
        shiftDist = stats::rnorm(1,mn_d,sd_d)

        ## compute the required alpha
        alpha = shiftDist/sqrt(sum( compositions::clr(pertb)^2 ))

        ## Get alpha
        if(stats::runif(1)>.5){
          alpha = alpha*-1
        }else{
          alpha = alpha
        }


        ## define new point in random direction of length shift dist @ alpha n
        ph = compositions::clo( pt  * compositions::clo(pertb^alpha) )


        if(is.null(labels)){
          alphaData_debug = data.frame(rep = r,sampleNum = i,
                                       comm = cid,
                                       actAlpha = alpha,
                                       act_shiftDistance = sqrt( sum( compositions::clr( compositions::clo(ph/pt))^2 ) ),
                                       diff = abs(sqrt(sum( compositions::clr( compositions::clo(ph/pt))^2 ) ) - shiftDist),
                                       reqShift = shiftDist)
        }else{
          alphaData_debug = data.frame(classLabel = phd.labels[i],rep = r,sampleNum = i,
                                       comm = cid,
                                       actAlpha = alpha,
                                       act_shiftDistance = sqrt(sum( compositions::clr( compositions::clo(ph/pt))^2 ) ),
                                       diff = abs(sqrt(sum( compositions::clr( compositions::clo(ph/pt))^2 ) ) - shiftDist),
                                       reqShift = shiftDist)
        }




        ## Combine Data
        cbind(alphaData_debug,ph)

      }

      allDat = rbind(allDat,ad)
      message("Rep = ",r ," of ",ptPerSample,"; Comm =  ",cid," of ",length(comIDs))
    }

  }



 if(is.null(labels)){
   classLabel = data.frame()
   reqAlpha = allDat[,1:7]
   shiftsPoints = allDat[,-7:-1]
 }else{
   classLabel = allDat[,1]
   allDat = allDat[,-1]
   reqAlpha = allDat[,1:7]
   shiftsPoints = allDat[,-7:-1]
 }

  return(list(simData = shiftsPoints, alphaData = reqAlpha,ClassLabels = classLabel ))
}

