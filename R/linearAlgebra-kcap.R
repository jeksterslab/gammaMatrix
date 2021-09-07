#' Commutation Matrix
#'
#' Create a commutation matrix.
#'
#' The commutation matrix \eqn{K_{m}} is the
#' \eqn{
#'   \frac{m \left( m + 1 \right)}{2} \times m^2
#' }
#' matrix
#' for a given \eqn{m \times m} symmetric matrix \eqn{A}
#' where
#' \eqn{
#'   K_{m} \mathrm{vec} \left( A \right)
#'   =
#'   \mathrm{vech} \left( A^{\prime} \right)
#' }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Integer.
#'   Dimension of the symmetric matrix.
#'
#' @references
#' [Wikipedia: Commutation matrix](https://en.wikipedia.org/wiki/Commutation_matrix)
#'
#' @returns A matrix.
#'
#' @examples
#' kcap(3)
#' @export
#' @family Linear Algebra Functions
#' @keywords linearAlgebra symmetric
kcap <- function(x) {
  dcap <- dcap(x)
  return(
    tcrossprod(
      solve(
        crossprod(dcap)
      ),
      dcap
    )
  )
}
