% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiNorm-mvn_theta_helper.R
\name{mvn_theta_helper}
\alias{mvn_theta_helper}
\title{Helper Function to Convert a Vector of Parameters
to the Mean Vector and Covariance Matrix}
\usage{
mvn_theta_helper(x)
}
\arguments{
\item{x}{Numeric vector.
Vector of parameters
\eqn{
    \boldsymbol{\theta}
    =
    \{
      \boldsymbol{\mu} ,
      \mathrm{vech} \left( \boldsymbol{\Sigma} \right)
    \}^{\prime}
  }.}
}
\value{
A list with elements \code{mu} and \code{sigmacap}.
}
\description{
Helper Function to Convert a Vector of Parameters
to the Mean Vector and Covariance Matrix
}
\examples{
mvn_theta_helper(
  x = c(0, 0, 1, 0.5, 1)
)
}
\seealso{
Other Multivariate Normal Distribution Functions: 
\code{\link{grad_l_mvn_generic}()},
\code{\link{grad_l_mvn}()},
\code{\link{hess_l_mvn_generic}()},
\code{\link{hess_l_mvn}()},
\code{\link{l_mvn_generic}()},
\code{\link{l_mvn}()},
\code{\link{negl_mvn}()},
\code{\link{rmvn_chol}()}
}
\concept{Multivariate Normal Distribution Functions}
\keyword{multiNorm}
