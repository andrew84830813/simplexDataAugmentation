% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simScenarios.R
\name{simDirDistr}
\alias{simDirDistr}
\title{Simulate Data From Dirichilet Distribution}
\usage{
simDirDistr(
  n1 = 30,
  dms_ = 75,
  seed = 8272008,
  concenParms = rep(1, dms_),
  scale = 1
)
}
\arguments{
\item{n1}{number of samples}

\item{dms_}{number of dimensions}

\item{seed}{random seed}

\item{concenParms}{parameters for dirichlet distribution}

\item{scale}{scale of concentration parms}
}
\value{
a data.frame with simulated data where the first column in the class name
}
\description{
Simulate a single class of data from a dirichliet distribution with
}
\examples{
\dontrun{
# simulate 30 samples from the Dirichlet ditriubtion with default parms
simDirDistr(n1 = 30)
}
}
