% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gammaMatrix-gammacap_gen.R
\name{gammacap_gen}
\alias{gammacap_gen}
\title{Asymptotic Covariance Matrix (General)}
\usage{
gammacap_gen(x, names = TRUE, sep = ".")
}
\arguments{
\item{x}{Numeric matrix, data frame, or vector.}

\item{names}{Logical.
Add names.}

\item{sep}{Character string.
Separator for variable names.}
}
\value{
A matrix.
}
\description{
Calculates the asymptotic covariance matrix of the unique elements
of the covariance matrix.
}
\section{Dependencies}{
\itemize{
\item \code{\link[=vech]{vech()}}
\item \code{\link[=vechnames]{vechnames()}}
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

gammacap_gen(x)
}
\references{
Yuan, K.-H., & Hayashi, K. (2006).
Standard errors in covariance structure models:
Asymptotics versus bootstrap.
British Journal of Mathematical and Statistical Psychology,
59, 397–417.
doi:\href{https://doi.org/10.1348/000711005X85896}{10.1348/000711005X85896}.
}
\seealso{
Other Gamma Matrix Functions: 
\code{\link{gammacap_adfnb}()},
\code{\link{gammacap_adf}()},
\code{\link{gammacap_mvnadj1}()},
\code{\link{gammacap_mvnadj2}()},
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
