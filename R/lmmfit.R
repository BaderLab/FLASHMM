#' Fitting Linear Mixed-effects Models
#'
#' @description lmmfit, a wrapper function of lmm, fits linear mixed-effects models (LMM) by sample-level data. The LMM parameters are estimated by either restricted maximum likelihood (REML) or maximum likelihood (ML) method with Fisher scoring (FS) gradient descent algorithm.
#'
#' @param X A design matrix for fixed effects, with rows corresponding to the columns of Y.
#' @param Y A features-by-samples matrix of responses (genes-by-cells matrix of gene expressions for scRNA-seq).
#' @param Z A design matrix for random effects, with rows corresponding to the columns of Y. Z = [Z1, ..., Zk], and Zi, i=1,...,k, is the design matrix for the i-th type random factor.
#' @param d = (d1,...,dk), where di = ncol(Zi), number of columns in Zi. sum(d) = ncol(Z), number of columns in Z. For the model with only one random factor, d = ncol(Z).
#' @param theta0 A vector of initial values of the variance components, (s1, ...,sk, s_(k+1)), si = sigma_i^2, the variance component of the i-th type random effects. s_(k+1) = sigma^2, the variance component of model residual error.
#' @param nBlocks Number of the blocks, which a big data is subdivided into, used for reducing the storage in computing the summary statistics that are computed from a block of data. The default value may not be adequate. If encountering the error: vector memory limit reached, you should increase the nBlocks value to avoid the issue.
#' @param method The REML with Fisher scoring (FS) iterative algorithm, REML-FS.
#' @param max.iter The maximal number of iterations for the iterative algorithm.
#' @param epsilon Positive convergence tolerance. If the absolute value of the first partial derivative of log likelihood is less than epsilon, the iterations converge.
#' @param output.cov If TRUE, output the covariance matrices for the estimated coefficients, which are needed for testing contrasts.
#' @param output.RE If TRUE, output the best linear unbiased prediction (BLUP) of the random effects.
#'
#' @return A list containing the following components:
#'    \item{dlogL}{First partial derivatives of log-likelihoods for each feature.}
#'    \item{logLik}{Maximum log-likelihoods for ML method or maximum log-restricted-likelihood for REML method.}
#'    \item{niter}{Numbers of iterations for each feature.}
#'    \item{coef}{A matrix of estimated coefficients (fixed effects), each column corresponds to a feature and each row one covariate.}
#'    \item{se}{A matrix of standard errors of the estimated coefficients.}
#'    \item{t}{A matrix of t-values for the fixed effects, equal to coef/se.}
#'    \item{df}{Degrees of freedom for the t-statistics (values).}
#'    \item{p}{A matrix of two-sided p-values for the t-tests of the fixed effects.}
#'    \item{cov}{A array of covariance matrices of the estimated coefficients (fixed effects).}
#'    \item{theta}{A matrix of the estimated variance components, each column corresponds to a feature and each row one variance component. The last row is the variance component of the residual error.}
#'    \item{se.theta}{Standard errors of the estimated theta.}
#'    \item{RE}{A matrix of the best linear unbiased prediction (BLUP) of random effects.}
#'
#' @importFrom MASS ginv
#' @importFrom stats pt
#' @import Matrix
#'
#' @seealso \code{\link{lmm}}
#' @examples
#' #Generate data: X, Y, and Z.
#' set.seed(2024)
#'
#' n <- 1e3
#' m <- 10
#' Y <- matrix(rnorm(m*n), m, n)
#' rownames(Y) <- paste0("Gene", 1:nrow(Y))
#'
#' trt <- sample(c("A", "B"), n, replace = TRUE)
#' X <- model.matrix(~ 0 + trt)
#'
#' q <- 20
#' sam <- rep(NA, n)
#' sam[trt == "A"] <- paste0("A", sample.int(q/2, sum(trt == "A"), replace = TRUE))
#' sam[trt == "B"] <- paste0("B", sample.int(q/2, sum(trt == "B"), replace = TRUE))
#' Z <- model.matrix(~ 0 + sam)
#' d <- ncol(Z)
#'
#' #Fit LMM by the cell-level data
#' fit <- lmmfit(Y, X, Z, d = d)
#' str(fit)
#'
#' #Fit LMM by summary-level data
#' #Compute and store the summary-level data:
#' n <- nrow(X)
#' XX <- t(X)%*%X
#' XY <- t(Y%*%X)
#' ZX <- t(Z)%*%X
#' ZY <- t(Y%*%Z)
#' ZZ <- t(Z)%*%Z
#' Ynorm <- rowSums(Y*Y)
#' fitss <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d)
#'
#' identical(fit, fitss)
#'
#' #Hypothesis testing
#' lmmtest(fit)
#' lmmtest(fit, index = 2)
#' lmmtest(fit, contrast = cbind("B-A" = c(-1, 1)))
#'
#' @export
lmmfit <- function(Y, X, Z, d = ncol(Z), theta0 = NULL, nBlocks = ceiling((ncol(Y)*1e-08)*nrow(Y)), method = c("REML", "ML"), max.iter = 50, epsilon = 1e-5, output.cov = TRUE, output.RE = FALSE)
{
stopifnot(!any(is.na(Y)), !any(is.na(X)), !any(is.na(Z)))
if (is.vector(Y)) Y <- t(Y)
stopifnot(ncol(Y) == nrow(X), ncol(Y) == nrow(Z))

method <- match.arg(method)

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

n <- nrow(X)
XX <- t(X)%*%X
ZX <- t(Z)%*%X
ZZ <- t(Z)%*%Z

lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d, theta0 = theta0, method = method, max.iter = max.iter, epsilon = epsilon, output.cov = output.cov, output.RE = output.RE)
}
