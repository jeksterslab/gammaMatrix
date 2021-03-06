% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linearAlgebra-kcap.R
\name{kcap}
\alias{kcap}
\title{Commutation Matrix}
\usage{
kcap(x)
}
\arguments{
\item{x}{Integer.
Dimension of the symmetric matrix.}
}
\value{
A matrix.
}
\description{
Create a commutation matrix.
}
\details{
The commutation matrix \eqn{K_{m}} is the
\eqn{
  \frac{m \left( m + 1 \right)}{2} \times m^2
}
matrix
for a given \eqn{m \times m} symmetric matrix \eqn{A}
where
\eqn{
  K_{m} \mathrm{vec} \left( A \right)
  =
  \mathrm{vech} \left( A^{\prime} \right)
}.
}
\examples{
kcap(3)
}
\references{
\href{https://en.wikipedia.org/wiki/Commutation_matrix}{Wikipedia: Commutation matrix}
}
\seealso{
Other Linear Algebra Functions: 
\code{\link{dcap}()},
\code{\link{sym_of_vech}()},
\code{\link{vechnames}()},
\code{\link{vech}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{Linear Algebra Functions}
\keyword{linearAlgebra}
\keyword{symmetric}
