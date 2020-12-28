\name{TEMM}
\alias{TEMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit the Tensor Envelope Mixture Model (TEMM)}

\description{Fit the Tensor Envelope Mixture Model (TEMM)}

\usage{
TEMM(Xn, u, K, initial = "kmeans", iter.max = 500, 
stop = 1e-3, trueY = NULL, print = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xn}{The tensor for clustering, should be array tpye, the last dimension is the sample size \code{n}.}
  \item{u}{A vector of envelope dimension}
  \item{K}{Nuber of clusters, greater than or equal to \code{2}.}
  \item{initial}{Initialization methed for the regularized EM algorithm. Default value is "kmeans".}
  \item{iter.max}{Maximum number of iterations. Default value is \code{500}.}
  \item{stop}{Convergence threshold of relative change in cluster means. Default value is \code{1e-3}.}
  \item{trueY}{A vector of true cluster labels of each observation. Default value is NULL.}
  \item{print}{Whether to print information including current iteration number, ralative change in cluster means 
  and clustering error (\code{\%}) in each iteration.}
}


\value{
  \item{id}{A vector of estimated labels.}
  \item{pi}{Estimated cluster wight.}
  \item{eta}{A \code{n} by \code{K} matrix of estimated membership weights.}
  \item{Mu.est}{A list of estimated cluster means.}
  \item{SIG.est}{A list of estimated covariance matrices.}
  \item{Mm}{Estimation of \code{Mm} defined in paper.}
  \item{Nm}{Estimation of \code{Nm} defined in paper.}
  \item{Gamma.est}{A list of estimated envelope basis.}
  \item{PGamma.est}{A list of envelope projection matrices.}
}

\references{Deng, K. and Zhang, X. (2020) \emph{Tensor Envelope Mixture Model For Simultaneous Clustering and Multiway Dimension Reduction, Manuscript.}}

\author{Kai Deng, Xin Zhang}


\examples{
  A = array(c(rep(1,20),rep(2,20))+rnorm(40),dim=c(2,2,10))
  myfit = TEMM(A,u=c(2,2),K=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
