#' Computing Summary-level Data from Individual-level Data
#'
#' @description sslmm can be used to compute the summary statistics (summary-level data) for \code{\link{lmm}} function, defined as
#' \itemize{
#' \item XX = t(X)\%*\%X
#' \item XY = t(X)\%*\%t(Y)
#' \item ZX = t(Z)\%*\%X
#' \item ZY = t(Z)\%*\%t(Y)
#' \item ZZ = t(Z)\%*\%Z
#' \item Ynorm = rowSums(Y*Y)
#' \item n = nrow(X)
#' }
#'
#' @param X A design matrix for fixed effects, with rows corresponding to the columns of Y.
#' @param Y A features-by-samples matrix of responses (genes-by-cells matrix of gene expressions for scRNA-seq).
#' @param Z A design matrix for random effects, with rows corresponding to the columns of Y.
#' @param nBlocks Number of the blocks, which a big data is subdivided into, used for reducing the storage in computing the summary statistics that are computed from a block of data. The default value may not be adequate. If encountering the error: vector memory limit reached, you should increase the nBlocks value to avoid the issue.
#'
#' @return A list of summary statistics:
#' XX, XY, ZX, ZY, ZZ, Ynorm and n.
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
#' s1 <- sslmm(X, Y, Z, nBlocks = 1)
#' s2 <- sslmm(X, Y, Z, nBlocks = 2)
#' s3 <- sslmm(X, Y, Z, nBlocks = 3)
#'
#' identical(s1, s2)
#' identical(s2, s3)
#'
#' @import Matrix
#'
#' @export
sslmm <- function(X, Y, Z, nBlocks = ceiling((ncol(Y)*1e-08)*nrow(Y))) {
stopifnot(!any(is.na(Y)), !any(is.na(X)), !any(is.na(Z)))
stopifnot(ncol(Y) == nrow(X), ncol(Y) == nrow(Z))

nr <- nrow(Y)
if (nBlocks > nr){
  message(paste0("Note: nBlocks=", nBlocks, " by default, changed to nrow(Y), i.e., nBlocks=", nr, "."))
  nBlocks <- nr
}
size <- round(nr/nBlocks)
if (nBlocks*size < nr) size <- round(nr/nBlocks) + 1

XY <- NULL
ZY <- NULL
Ynorm <- NULL
for (i in 1:nBlocks){
  j <- (1+(i-1)*size):(min(nr, i*size))
  XY <- cbind(XY, t(Y[j, , drop = FALSE]%*%X))
  ZY <- cbind(ZY, t(Y[j, , drop = FALSE ]%*%Z))
  Ynorm <- c(Ynorm, rowSums(Y[j, , drop = FALSE]*Y[j, , drop = FALSE]))
}

list(n = nrow(X),
	XX = t(X)%*%X,
	XY = XY,
	ZX = t(Z)%*%X,
	ZY = ZY,
	ZZ = t(Z)%*%Z,
	Ynorm = Ynorm)
}
