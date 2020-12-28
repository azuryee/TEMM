\name{tune_u_joint}
\alias{tune_u_joint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tuning envelope dimension jointly by BIC in TEMM.}

\description{Tuning envelope dimension jointly by BIC in TEMM.}

\usage{
tune_u_joint(u_candi, K, X, iter.max = 500, stop = 0.001, trueY = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u_candi}{A list of length \code{M} containing candidate envelope dimension for each mode.}
  \item{K}{Nuber of clusters, greater than or equal to \code{2}.}
  \item{X}{The tensor for clustering, should be array tpye, the last dimension is the sample size \code{n}.}
  \item{iter.max}{Maximum number of iterations, default is \code{500}.}
  \item{stop}{Convergence threshold of relative change in cluster means, default is \code{1e-3}.}
  \item{trueY}{A vector of true cluster labels of each observation, default is NULL.}
}


\value{
  \item{opt.u}{Optimal envelope dimension selected.}
  \item{opt.id}{Estimated labels fitted by TEMM with the optimal envelope dimension.}
  \item{opt.Mu}{Estimated cluster means fitted by TEMM with the optimal envelope dimension.}
  \item{bic}{BIC value.}
}


\references{Deng, K. and Zhang, X. (2020) \emph{Tensor Envelope Mixture Model For Simultaneous Clustering and Multiway Dimension Reduction, Manuscript.}}

\author{Kai Deng, Xin Zhang}


\examples{
  A = array(c(rep(1,20),rep(2,20))+rnorm(40),dim=c(2,2,10))
  mytune = tune_u_joint(u_candi=list(1:2,1:2),K=2,A)
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

