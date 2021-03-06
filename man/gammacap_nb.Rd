% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gammaMatrix-gammacap_nb.R
\name{gammacap_nb}
\alias{gammacap_nb}
\title{Nonparametric Bootstrapped Covariance Matrix}
\usage{
gammacap_nb(x, bcap = 1000L, seed = NULL, names = TRUE, sep = ".")
}
\arguments{
\item{x}{Numeric matrix, data frame, or vector.}

\item{bcap}{Integer.
Number of bootstrap samples.}

\item{seed}{Integer.
Random number generation seed.}

\item{names}{Logical.
Add names.}

\item{sep}{Character string.
Separator for variable names.}
}
\value{
A matrix.
}
\description{
Calculates the nonparametric bootstrapped
asymptotic distribution-free (ADF)
covariance matrix of the unique elements
of the covariance matrix.
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

gammacap_adfnb(x)
}
\references{
Browne, M. W. (1984).
Asymptotically distribution-free methods
for the analysis of covariance structures.
British Journal of Mathematical and Statistical Psychology,
37, 62–83.
doi:\href{https://doi.org/10.1111/j.2044-8317.1984.tb00789.x}{10.1111/j.2044-8317.1984.tb00789.x}.

Yung, Y.-F., & Bentler, P. M. (1994).
Bootstrap-corrected ADF test statistics
in covariance structure analysis.
British Journal of Mathematical and Statistical Psychology,
47, 63–84.
doi:\href{https://doi.org/10.1111/j.2044-8317.1994.tb01025.x}{10.1111/j.2044-8317.1994.tb01025.x}.
}
\seealso{
Other Gamma Matrix Functions: 
\code{\link{gammacap_adfnb}()},
\code{\link{gammacap_adf}()},
\code{\link{gammacap_gen}()},
\code{\link{gammacap_mvnadj1}()},
\code{\link{gammacap_mvnadj2}()},
\code{\link{gammacap_mvn}()},
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
