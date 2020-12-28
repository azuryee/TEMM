\name{TGMM}
\alias{TGMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit the Tensor Gaussian Mixture Model (TGMM)}

\description{Fit the Tensor Gaussian Mixture Model (TGMM)}

\usage{
TGMM(Xn, K, shape = "shared", initial = "kmeans", 
iter.max = 500, stop = 1e-3, trueY = NULL, print = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xn}{The tensor for clustering, should be array tpye, the last dimension is the sample size \code{n}.}
  \item{K}{Nuber of clusters, greater than or equal to \code{2}.}
  \item{shape}{"shared" if assume common covariance across mixtures, "distinct" if allow different covariance structures. Default value is "shared".}
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
}

\references{Deng, K. and Zhang, X. (2020) \emph{Tensor Envelope Mixture Model For Simultaneous Clustering and Multiway Dimension Reduction, Manuscript.}}

\author{Kai Deng, Xin Zhang}


\examples{
  A = array(c(rep(1,20),rep(2,20))+rnorm(40),dim=c(2,2,10))
  myfit = TGMM(A,K=2,shape="shared")
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

