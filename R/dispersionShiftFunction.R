#' Alter Mulitvariate Dispersion of Point CLoud in SImplex
#'
#' @param tbl2 a closed and imputed data matrix of points cloud
#' @param dispFact Factor to increase (dispFact>0) or decrease (dispFact<0) dispersion by
#'
#' @return A list containing:\tabular{ll}{
#'    \code{combinedOutput} \tab xxx xx xx  \cr
#'    \tab \cr
#'    \code{dispAdjustedData} \tab xxx xxx \cr
#'     \tab \cr
#'    \code{dispDistr.shifted} \tab xxx \cr
#'     \tab \cr
#'    \code{dispDistr.adj} \tab xxx \cr
#' }
#' @export
dispersionShift = function(tbl2,dispFact = 1){



  ## centroids
  compMean.shifted = compositions::mean.acomp(compositions::acomp(tbl2))


  ### dispersion of the shifted distribution
  dspLen_shifted = data.frame()
  for(i in 1:nrow(tbl2)){
    ph = compositions::clo(as.numeric(tbl2[i,])/as.numeric(compMean.shifted))
    ph = data.frame(pointNum = i,disCentr = sqrt(sum( compositions::clr( ph)^2 ) ))
    dspLen_shifted = rbind(dspLen_shifted,ph)
  }


  #hist(dspLen.fixed$disCentr)
  d = dspLen_shifted$disCentr



  ########################################################
  shiftsPoints = foreach::foreach(i = 1:nrow(tbl2),.combine = rbind)%dopar%{

    ## Define base point
    pt = ((as.matrix(tbl2[i,])))
    ## compute linear shift direction
    pertb = compositions::clo(as.numeric(compMean.shifted) / as.numeric(pt))

    ## sample dispersion distance
    shiftDist = d[i]*stats::rnorm(1,dispFact,dispFact/3)

    # compute required alpha
    alpha = shiftDist/sqrt(sum( compositions::clr( pertb)^2 ))


    ## compute final shifted point at the solved alpha
    ph = compositions::clo( as.numeric(compMean.shifted)  * compositions::clo(pertb^alpha) )


    ph = data.frame(sampleNum = i,
                    reqAlpha = alpha,
                    actDist =  sqrt(sum(compositions::clr(compositions::clo(ph/as.numeric(compMean.shifted)))^2)),
                    requiredDist = shiftDist,
                    t(ph))

    ph

  }



  pointAtAlpha_n = shiftsPoints[,-4:-1]
  colnames(pointAtAlpha_n) = colnames(tbl2)

  ### New Dispersion Distr
  dspLen_after = data.frame()
  for(i in 1:nrow(pointAtAlpha_n)){
    ph = compositions::clo(as.numeric(pointAtAlpha_n[i,])/as.numeric(compMean.shifted))
    ph = data.frame(pointNum = i,disCentr = sqrt(sum( compositions::clr(ph)^2 )))
    dspLen_after = rbind(dspLen_after,ph)
  }




  #dispersion test
  simDat = rbind(
    data.frame(Type = "dispAdjusted",pointAtAlpha_n),
    data.frame(Type = "Base",tbl2))



  return(list(combinedOutput = simDat,
              dispAdjustedData = pointAtAlpha_n,
              dispDistr.shifted = dspLen_shifted,
              dispDistr.adj = dspLen_after   ))
}
