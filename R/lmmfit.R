#' Fitting Linear Mixed-effects Models
#'
#' @description lmmfit, a wrapper function of lmm, fits linear mixed-effects models (LMM) by sample-level data. The LMM parameters are estimated by restricted maximum likelihood (REML) with Fisher scoring (FS) gradient descent algorithm.
#'
#' @param X A design matrix for fixed effects, with rows corresponding to the columns of Y.
#' @param Y A features-by-samples matrix of responses (genes-by-cells matrix of gene expressions for scRNA-seq).
#' @param Z A design matrix for random effects, with rows corresponding to the columns of Y. Z = [Z1, ..., Zk], and Zi, i=1,...,k, is the design matrix for the i-th type random factor.
#' @param d A vector of (m1,...,mk), mi = ncol(Zi), number of columns in Zi. m1 + ... + mk = ncol(Z), number of columns in Z.
#' @param theta0 A vector of initial values of the variance components, (s1, ...,sk, s_(k+1)), si = sigma_i^2, the variance component of the i-th type random effects. s_(k+1) = sigma^2, the variance component of model residual error.
#' @param nBlocks Number of blocks, used for blocking a big data to reduce the storage required in computing.
#' @param method The REML with Fisher scoring (FS) iterative algorithm, REML-FS.
#' @param max.iter The maximal number of iterations for the iterative algorithm.
#' @param epsilon Positive convergence tolerance. If the absolute value of the first partial derivative of log likelihood is less than epsilon, the iterations converge.
#' @param output.cov If TRUE, output the covariance matrices for the estimated coefficients, which are needed for testing contrasts.
#' @param output.RE If TRUE, output the best linear unbiased prediction (BLUP) of the random effects.
#'
#' @return A list containing the following components:
#' @return dlogL First partial derivatives of log-likelihoods for each feature (gene).
#' @return niter Nmbers of iterations for each feature (gene).
#' @return coef A matrix of estimated coefficients (fixed effects), each column corresponds to a feature (gene) and each row one covariate.
#' @return se A matrix of the standard errors of the estimated coefficients.
#' @return t A matrix of t-values for the fixed effects, equal to coef/se.
#' @return df Degrees of freedom.
#' @return p A matrix of two-sided p-values for the fixed effects.
#' @return cov A array of covariance matrices of the estimated coefficients (fixed effects).
#' @return theta A matrix of the estimated variance components, each column corresponds to a feature (gene) and each row one variance component. The last row is the variance component of the residual error.
#' @return se.theta Standard errors of the estimated theta.
#' @return RE A matrix of the best linear unbiased prediction (BLUP) of random effects.
#'
#' @importFrom MASS ginv
#' @importFrom stats pt
#' @import Matrix
#'
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
#' #Fitting LMM by the cell-level data
#' fit <- lmmfit(Y, X, Z, d = d)
#' str(fit)
#'
#' #Hypothesis tests
#' lmmtest(fit)
#' lmmtest(fit, index = 2)
#' lmmtest(fit, contrast = cbind("B-A" = c(-1, 1)))
#'
#' @export
lmmfit <- function(Y, X, Z, d, theta0 = NULL, nBlocks = 1, method = "REML-FS", max.iter = 50, epsilon = 1e-5, output.cov = TRUE, output.RE = FALSE)
{
stopifnot(!any(is.na(Y)), !any(is.na(X)), !any(is.na(Z)))
stopifnot(ncol(Y) == nrow(X), ncol(Y) == nrow(Z))

nr <- nrow(Y)
size <- round(nr/nBlocks)
if (nBlocks*size < nr) size <- round(nr/nBlocks) + 1

XY <- NULL
ZY <- NULL
Ynorm <- NULL
for (i in 1:nBlocks){
  j <- (1+(i-1)*size):(min(nr, i*size))
  XY <- cbind(XY, t(Y[j, ]%*%X))
  ZY <- cbind(ZY, t(Y[j, ]%*%Z))
  Ynorm <- c(Ynorm, rowSums(Y[j, ]*Y[j, ]))
}

n <- nrow(X)
XX <- t(X)%*%X
ZX <- t(Z)%*%X
ZZ <- t(Z)%*%Z

lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d, theta0 = theta0, method = method, max.iter = max.iter, epsilon = epsilon, output.cov = output.cov, output.RE = output.RE)
}
