#' Duplication Matrix
#'
#' Create a duplication matrix.
#'
#' The duplication matrix \eqn{D_{m}} is the
#' \eqn{
#'   m^2 \times \frac{m \left( m + 1 \right)}{2}
#' }
#' matrix
#' for a given \eqn{m \times m} symmetric matrix \eqn{A}
#' where
#' \eqn{
#'   D_{m} \mathrm{vech} \left( A \right)
#'   =
#'   \mathrm{vec} \left( A \right)
#' }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Integer.
#'   Dimension of the symmetric matrix.
#'
#' @references
#' [Wikipedia: Duplication matrix](https://en.wikipedia.org/wiki/Duplication_and_elimination_matrices#Duplication_matrix)
#'
#' @returns A matrix.
#'
#' @examples
#' dcap(3)
#' @export
#' @family Linear Algebra Functions
#' @keywords linearAlgebra symmetric
dcap <- function(x) {
  stopifnot(
    is.vector(x),
    length(x) == 1
  )
  m <- as.integer(x)
  sym <- diag(m)
  i <- seq_len(
    0.5 * m * (m + 1)
  )
  sym[lower.tri(sym, diag = TRUE)] <- i
  sym[upper.tri(sym)] <- t(sym)[upper.tri(sym)]
  return(
    outer(
      X = c(sym),
      Y = i,
      FUN = function(x, y) {
        ifelse(
          test = x == y,
          yes = 1,
          no = 0
        )
      }
    )
  )
}
