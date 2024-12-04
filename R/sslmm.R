#' Computing Summary-level Data from Individual-level Data
#'
#' @description sslmm can be used to compute the correlation-related summary statistics (summary-level data) for lmm function.
#'
#' @param X A design matrix for fixed effects, with rows corresponding to the columns of Y.
#' @param Y A features-by-samples matrix of responses (genes-by-cells matrix of gene expressions for scRNA-seq).
#' @param Z A design matrix for random effects, with rows corresponding to the columns of Y.
#'
#' @return A list of summary statistics:
#' XX = t(X)\%*\%X,
#' XY = t(X)\%*\%t(Y),
#' ZX = t(Z)\%*\%X,
#' ZY = t(Z)\%*\%t(Y),
#' ZZ = t(Z)\%*\%Z,
#' Ynorm = rowSums(Y*Y) and n = nrow(X).
#'
#' @examples
#' n <- 1e3
#' set.seed(2024)
#' p <- 2
#' X <- matrix(rnorm(p*n), n, p)
#' colnames(X) <- paste0("X", 1:ncol(X))
#' m <- 3
#' Y <- matrix(rnorm(m*n), m, n)
#' rownames(Y) <- paste0("Y", 1:nrow(Y))
#' q <- 4
#' Z <- gl(q, n/q, labels = letters[1:q])
#' Z <- model.matrix(~ 0 + Z)
#' sslmm(X, Y, Z)
#'
#' @import Matrix
#'
#' @export
sslmm <- function(X, Y, Z) {
stopifnot(!any(is.na(Y)), !any(is.na(X)), !any(is.na(Z)))
stopifnot(ncol(Y) == nrow(X), ncol(Y) == nrow(Z))

list(n = nrow(X),
	XX = t(X)%*%X,
	XY = t(Y%*%X),
	ZX = t(Z)%*%X,
	ZY = t(Y%*%Z),
	ZZ = t(Z)%*%Z,
	Ynorm = rowSums(Y*Y))
}
