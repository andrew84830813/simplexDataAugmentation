#' Simulate Data from Additive Logistics Function
#'
#'Simulates data from the additive logistics normal distribution
#'
#'
#'
#' @param n1 number of samples to return
#' @param dms_ number of dimesions to simulate
#' @param className name of class from simulated data
#' @param seed random seed for simulation
#' @param varParm Controls the variance of simulated distribution
#' @param meanVec Sets the value of the mean vector for the multivariate distribution
#'
#' @return a data.frame with simulated data where the first column in the class name
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate 30 samples from the Dirichlet ditriubtion with default parms
#' simAddLogNormal(n1 = 30)
#' }
simAddLogNormal = function(n1 = 30, dms_ = 75,className = "S1",seed= 08272008,varParm = 2,meanVec = 0){
  set.seed(seed)
  mu_ = rep(meanVec,dms_)
  sigma_ = diag(dms_)
  U = matrix(stats::runif(dms_*dms_,0,varParm),nrow = dms_)
  U_ = sigma_ + U
  U_ = Matrix::nearPD(U_)$mat
  eg = min(eigen(as.matrix(sigma_))$values)
  sig = min(eigen(as.matrix(U_))$values)
  dd = min(eg,sig)+0.05
  sig1 = Matrix::nearPD( U_ + sigma_ + dd*diag(dms_) )$mat
  s1 = selEnergyPermR::sampleAddLogisticNormal(n = n1,dims_ = dms_,mu = mu_ ,sigma = sig1,sigmaScale = 1,sampleName = className)
  s1
}



#' Simulate Data additive logisitcs T
#'
#' @param n number of samples to return
#' @param df_ degrees of freedom
#' @param dms number of dimensions
#' @param sigmaScale scale covariance matrix
#' @param sampleName class name
#' @param seed random seed
#'
#' @return a data.frame with simulated data where the first column in the class name
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate 30 samples from the Dirichlet ditriubtion with default parms
#' sampleAddLogisticT(n1 = 30)
#' }
sampleAddLogisticT <-
  function(n, df_,dms,sigmaScale=1,sampleName = "S1",seed= 08272008){
    set.seed(seed)
    sigma = diag(dms)*sigmaScale
    y1 = mvtnorm::rmvt(n,sigma ,df = df_,type = "Kshirsagar")
    y1 = compositions::alrInv(y1)
    s1  = data.frame(Status = sampleName,y1)
    return(s1)
  }





#' Simulate Data From Dirichilet Distribution
#'
#'Simulate a single class of data from a dirichliet distribution with
#'
#' @param n1 number of samples
#' @param dms_ number of dimensions
#' @param seed random seed
#' @param concenParms parameters for dirichlet distribution
#' @param scale scale of concentration parms
#'
#' @return a data.frame with simulated data where the first column in the class name
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate 30 samples from the Dirichlet ditriubtion with default parms
#' simDirDistr(n1 = 30)
#' }
simDirDistr <-
  function( n1 = 30 ,dms_ = 75,seed= 08272008,concenParms = rep(1,dms_), scale = 1 ){

    ## define and scale alpa
    a = concenParms*scale
    set.seed(seed)
    s1 = selEnergyPermR::sampleDirichlet(n1 = n1,dims = dms_,sampleName = "S1",a1 = a)$Sample

    return(s1)

  }
