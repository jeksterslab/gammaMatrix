% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gammaMatrix-gammacap_mvnadj2.R
\name{gammacap_mvnadj2}
\alias{gammacap_mvnadj2}
\title{Asymptotic Covariance Matrix with Adjustment (Variant 2)}
\usage{
gammacap_mvnadj2(
  x,
  missing = FALSE,
  ml_cov = TRUE,
  drop_means = TRUE,
  names = TRUE,
  sep = "."
)
}
\arguments{
\item{x}{Numeric matrix, data frame, or vector.}

\item{missing}{Logical.
If \code{missing = TRUE},
the mean vector and the covariance matrix
will be estimated using the EM algorithm.
If \code{missing = FALSE},
all missing values will be dropped
and the mean vector and covariance matrix
will be estimated using \code{x}.}

\item{ml_cov}{Logical.
If \code{missing = FALSE} and \code{ml_cov = TRUE},
use maximum likelihood estimator of the covariance matrix.}

\item{drop_means}{Logical.
If \code{drop_means = TRUE},
remove the means from the output matrix.
If \code{drop_means = FALSE},
the first \code{ncol(x)} elements of the output matrix
will be names for the \code{ncol(x)} means.}

\item{names}{Logical.
Add names.}

\item{sep}{Character string.
Separator for variable names.}
}
\value{
A matrix.
}
\description{
Calculates the covariance matrix of the unique elements
of the covariance matrix with adjustment for nonnormality.
}
\section{Dependencies}{
\itemize{
\item \code{\link[=rmvn_chol]{rmvn_chol()}} (test)
}
}

\examples{
set.seed(42)
n <- 1000
k <- 2
z <- matrix(
  data = rnorm(n = n * k), nrow = n, ncol = k
)
q <- chol(
  matrix(
    data = c(1.0, 0.5, 0.5, 1.0),
    nrow = k, ncol = k
  )
)
x <- z \%*\% q

gammacap_mvnadj2(x)
gammacap_mvnadj2(x, drop_means = FALSE)
}
\references{
Add appropriate references here...
}
\seealso{
Other Gamma Matrix Functions: 
\code{\link{gammacap_adfnb}()},
\code{\link{gammacap_adf}()},
\code{\link{gammacap_gen}()},
\code{\link{gammacap_mvnadj1}()},
\code{\link{gammacap_mvn}()},
\code{\link{gammacap_nb}()},
\code{\link{gammacap_ols_generic}()},
\code{\link{gammacap_ols_hc_generic}()},
\code{\link{gammacap_ols_hc_qcap_generic}()},
\code{\link{gammacap_ols_hc_qcap}()},
\code{\link{gammacap_ols_hc}()},
\code{\link{gammacap_ols}()},
\code{\link{gammacapnames}()},
\code{\link{gammacap}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{Gamma Matrix Functions}
\keyword{gammaMatrix}
