% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiNorm-grad_l_mvn_generic.R
\name{grad_l_mvn_generic}
\alias{grad_l_mvn_generic}
\title{Gradient Vector of the Multivariate Normal Distribution - Generic}
\usage{
grad_l_mvn_generic(x, mu, sigmacap)
}
\arguments{
\item{x}{Numeric vector of length \code{k}.
The ith vector of observations.}

\item{mu}{Numeric vector.
Parameter.
Mean vector
\eqn{\boldsymbol{\mu}}.}

\item{sigmacap}{Numeric matrix.
Parameter.
Covariance matrix
\eqn{\boldsymbol{\Sigma}}.}
}
\value{
A vector.
}
\description{
Calculates gradient vector of the log of the likelihood function
of the multivariate normal distribution
for the ith observation.
}
\examples{
n <- 5
mu <- c(0, 0)
sigmacap <- matrix(
  data = c(
    1, 0.5, 0.5, 1
  ),
  nrow = 2
)

xcap <- as.data.frame(
  t(
    rmvn_chol(
      n = n,
      mu = mu,
      sigmacap = sigmacap
    )
  )
)

lapply(
  X = xcap,
  FUN = grad_l_mvn_generic,
  mu = mu,
  sigmacap = sigmacap
)
}
\seealso{
Other Multivariate Normal Distribution Functions: 
\code{\link{grad_l_mvn}()},
\code{\link{hess_l_mvn_generic}()},
\code{\link{hess_l_mvn}()},
\code{\link{l_mvn_generic}()},
\code{\link{l_mvn}()},
\code{\link{mvn_theta_helper}()},
\code{\link{negl_mvn}()},
\code{\link{rmvn_chol}()}
}
\concept{Multivariate Normal Distribution Functions}
\keyword{multiNorm}
