\name{tune_u_sep}
\alias{tune_u_sep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tuning envelope dimension separately by BIC in TEMM.}

\description{Tuning envelope dimension separately by BIC in TEMM.}

\usage{
tune_u_sep(m, u_candi, K, X, C = 1, oneD = TRUE, 
iter.max = 500, stop = 0.001, trueY = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{The tensor mode to be tuned, can take value in \code{1,...,M}.}
  \item{u_candi}{A vector of candidate envelope dimension.}
  \item{K}{Nuber of clusters, greater than or equal to \code{2}.}
  \item{X}{The tensor for clustering, should be array tpye, the last dimension is the sample size \code{n}.}
  \item{C}{Constant in separate BIC creterion, default is \code{1}.}
  \item{oneD}{Whether to apply 1D-BIC tuning, default is TRUE.}
  \item{iter.max}{Maximum number of iterations, default is \code{500}.}
  \item{stop}{Convergence threshold of relative change in cluster means, default is \code{1e-3}.}
  \item{trueY}{A vector of true cluster labels of each observation, default is NULL.}
}

\value{
  \item{opt.u}{Optimal envelope dimension selected.}
  \item{bic}{BIC value.}
}


\references{Deng, K. and Zhang, X. (2020) \emph{Tensor Envelope Mixture Model For Simultaneous Clustering and Multiway Dimension Reduction, Manuscript.}}

\author{Kai Deng, Xin Zhang}


\examples{
  A = array(c(rep(1,20),rep(2,20))+rnorm(40),dim=c(2,2,10))
  mytune = tune_u_sep(1,1:2,K=2,A)
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
