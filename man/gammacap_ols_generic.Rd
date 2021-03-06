% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gammaMatrix-gammacap_ols_generic.R
\name{gammacap_ols_generic}
\alias{gammacap_ols_generic}
\title{Asymptotic Covariance Matrix for Ordinary Least Squares Regression - Generic}
\usage{
gammacap_ols_generic(
  x,
  beta,
  sigmacap,
  sigmasq,
  ke,
  type = "gen",
  adf_unbiased = TRUE,
  bcap = 5000L,
  seed = NULL,
  yc = FALSE,
  names = TRUE,
  sep = "."
)
}
\arguments{
\item{x}{Numeric matrix, data frame, or vector.
Values for regressor variables.}

\item{beta}{Numeric vector.
Partial regression slopes.}

\item{sigmacap}{Numeric matrix.
Covariance matrix
\eqn{\boldsymbol{\Sigma}}
of
\eqn{\{y, x_1, \cdots, x_p \}^{\prime}}.}

\item{sigmasq}{Numeric.
Error variance
\eqn{\sigma^2}.}

\item{ke}{Numeric.
Kurtosis of errors.}

\item{type}{Character string.
Type of asymptotic covariance matrix of the regressors
If \code{type = "adf"},
calculate asymptotic distribution-free covariance matrix.
If \code{type = "adfnb"},
calculate nonparametric bootstrapped
asymptotic distribution-free covariance matrix.
If \code{type = "gen"},
calculate covariance matrix using the general formula.
If \code{type = "mvn"},
calculate covariance matrix with multivariate normal data.}

\item{adf_unbiased}{Logical.
If \code{adf_unbiased = TRUE},
use unbiased asymptotic distribution-free covariance matrix.
If \code{adf_unbiased = FALSE},
use consistent asymptotic distribution-free covariance matrix.}

\item{bcap}{Integer.
Number of bootstrap samples
when \code{type = "gen"} and \code{type = "adfnb"}.}

\item{seed}{Integer.
Random number generation seed
when \code{type = "gen"} and \code{type = "adfnb"}.}

\item{yc}{Logical.
If \code{yc = TRUE}, order the output
following Yuan and Chan (2011) page 674.
If \code{yc = FALSE}, the order of the output
following \eqn{\mathrm{vech} \left( \Sigma \right)}.}

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
of the covariance matrix using ordinary least squares estimates.
}
\section{Dependencies}{
\itemize{
\item \code{\link[=dcap]{dcap()}}
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
x <- as.data.frame(z \%*\% q)
colnames(x) <- c("y", "x")
obj <- lm(y ~ x, data = x)
beta <- coef(obj)[-1]
sigmacap <- cov(x)
sigmasq <- summary(obj)$sigma^2
ke <- 3

gammacap_ols_generic(
  x[, -1],
  beta = beta,
  sigmacap = sigmacap,
  sigmasq = sigmasq,
  ke = ke
)
}
\references{
Yuan, K.-H., & Chan, W. (2011).
Biases and standard errors of standardized regression coefficients.
Psychometrika,
76, 670–690.
doi:\href{https://doi.org/10.1007/S11336-011-9224-6}{10.1007/S11336-011-9224-6}.
}
\seealso{
Other Gamma Matrix Functions: 
\code{\link{gammacap_adfnb}()},
\code{\link{gammacap_adf}()},
\code{\link{gammacap_gen}()},
\code{\link{gammacap_mvnadj1}()},
\code{\link{gammacap_mvnadj2}()},
\code{\link{gammacap_mvn}()},
\code{\link{gammacap_nb}()},
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
