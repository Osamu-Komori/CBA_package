#' Diagnose approximate normality from skewness and kurtosis
#'
#' Computes sample skewness and excess kurtosis, and checks whether they fall
#' in a simple acceptable range.
#'
#' @param s A numeric vector.
#'
#' @return A list with components:
#' \describe{
#'   \item{skewness}{Sample skewness.}
#'   \item{excess_kurtosis}{Sample excess kurtosis.}
#'   \item{acceptable}{Logical value indicating whether
#'   \code{-1 <= skewness <= 1} and \code{-2 <= excess_kurtosis <= 5}.}
#' }
#'
#' @export
normality_diag <- function(s) {
  s <- as.numeric(s)
  s <- s[!is.na(s)]

  if (length(s) < 2) {
    stop("s must contain at least two non-missing numeric values.")
  }

  m <- mean(s)
  sdv <- sd(s)

  if (sdv == 0) {
    stop("s must have non-zero standard deviation.")
  }

  skew <- mean(((s - m) / sdv)^3)
  exkurt <- mean(((s - m) / sdv)^4) - 3

  acceptable <- (-1 <= skew && skew <= 1 &&
                 -2 <= exkurt && exkurt <= 5)

  list(
    skewness = skew,
    excess_kurtosis = exkurt,
    acceptable = acceptable
  )
}