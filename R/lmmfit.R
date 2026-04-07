#' Fitting Linear Mixed-effects Models
#'
#' @description lmmfit, a wrapper function of lmm, fits linear mixed-effects models (LMM) by sample-level data. The LMM parameters are estimated by either restricted maximum likelihood (REML) or maximum likelihood (ML) method with Fisher scoring (FS) gradient descent algorithm.
#'
#' @param Y a features-by-samples matrix of responses. For single-cell RNA sequencing (scRNA-seq) data, Y is a genes-by-cells gene expression matrix, whose entries can be raw counts, log-transformed values, or normalized values. If Y is a matrix of raw counts, it will be log-transformed as log2(1 + Y) for differential expression analysis.
#' @param X a design matrix for fixed effects, with rows corresponding to the columns of Y.
#' @param Z [Z1, ..., Zk] is a design matrix for k random components (factors), with rows corresponding to the columns of Y.
#' @param d (d1,...,dk), where di = ncol(Zi), the number of columns in Zi. sum(d) = ncol(Z), the number of columns in Z. For the model with only one random component, d = ncol(Z).
#' @param method a character string, either "REML" or "ML". Defaults to "REML". If ‘REML’, the model is fit using REML; otherwise, using ML.
#' @param nBlocks the number of blocks into which large datasets are subdivided to reduce memory usage when computing summary statistics. The default value, nBlocks = ceiling((ncol(Y)*1e-08)*nrow(Y)), may not be adequate. If encountering the error: vector memory limit reached, you should increase the nBlocks value to avoid the issue.
#' @param theta0 a vector of initial values of the variance components, (s1, ...,sk, s_(k+1)), si = sigma_i^2, the i-th variance component. s_(k+1) = sigma^2, the variance component of model residual error.
#' @param max.iter the maximal number of iterations for the iterative algorithm.
#' @param epsilon positive convergence tolerance. If the absolute value of the first partial derivative of log likelihood is less than epsilon, the iterations converge.
#' @param output.cov a logical scalar. If TRUE, output an array of the covariance matrices for the estimated coefficients, cov.beta.
#' @param output.FIM a logical scalar. If TRUE, output an array of the Fisher information matrices (FIM) for fixed effects, FIM.beta[,,i] = \eqn{X^T(Cov(Y[i,])^{-1}X} for the i-th feature, and variance components, FIM.theta[,,i] = \eqn{Var(\partial(logL_i)/\partial(theta)) = - E[\partial^2(logL_i)/{\partial(theta)\partial(theta^T)}]} for the i-th feature.
#' @param output.RE a logical scalar. If TRUE, output the best linear unbiased prediction (BLUP) of the random effects.
#' @param output.SS a logical scalar. If TURE, output list(XX, XY, ZX, ZY, ZZ, Ynorm, n), a list of the summary-level data (summary statistics), defined by XX=t(X)\%*\%X, XY=t(Y\%*\%X), ZX=t(Z)\%*\%X, ZY=t(Z)\%*\%Y, ZZ=t(Z)\%*\%Z, Ynorm=rowSums(Y*Y), and n=nrow(X), which can also be computed using the sslmm function.
#' @param verbose a logical scalar. If TRUE, print the number of features for which LMM fitting did not converge (abs(dlogL) > epsilon).
#'
#' @return A list containing the following components:
#'    \item{method}{the method, either REML or ML, for estimating the LMM parameters (fixed effects and variance components).}
#'    \item{epsilon}{convergence tolerance.}
#'    \item{dlogL}{first partial derivatives of log-likelihoods for each feature.}
#'    \item{logLik}{maximum log-likelihoods for ML method or maximum log-restricted-likelihood for REML method.}
#'    \item{niter}{numbers of iterations for each feature.}
#'    \item{coef}{a matrix of estimated coefficients (fixed effects or beta), each column corresponds to a feature and each row one covariate.}
#'    \item{se}{a matrix of standard errors of the estimated coefficients.}
#'    \item{t}{a matrix of t-values for the fixed effects, equal to coef/se.}
#'    \item{p}{a matrix of two-sided p-values for the t-tests of the fixed effects.}
#'    \item{df}{degrees of freedom for the t-statistics (values).}
#'    \item{cov}{a array of covariance matrices of the estimated coefficients (fixed effects).}
#'    \item{theta}{a matrix of the estimated variance components, each column corresponds to a feature and each row one variance component. The last row is the variance component of the residual error.}
#'    \item{se.theta}{standard errors of the estimated theta.}
#'	  \item{FIM.beta}{the Fisher information matrices (FIM) for fixed effects.}
#'	  \item{FIM.theta}{the Fisher information matrices (FIM) for variance components.}
#'    \item{RE}{a matrix of the best linear unbiased prediction (BLUP) of random effects.}
#'    \item{summary_data}{a list of the summary-level data (summary statistics).}
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
lmmfit <- function(Y, X, Z, d = ncol(Z), method = c("REML", "ML"), nBlocks = ceiling((ncol(Y)*1e-08)*nrow(Y)), theta0 = NULL, max.iter = 50, epsilon = 1e-5, output.cov = TRUE, output.FIM = FALSE, output.RE = FALSE, output.SS = FALSE, verbose = FALSE)
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
  ZY <- cbind(ZY, t(Y[j, , drop = FALSE]%*%Z))
  Ynorm <- c(Ynorm, rowSums(Y[j, , drop = FALSE]*Y[j, , drop = FALSE]))
}

n <- nrow(X)
XX <- t(X)%*%X
ZX <- t(Z)%*%X
ZZ <- t(Z)%*%Z

lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d, theta0 = theta0, method = method, max.iter = max.iter, epsilon = epsilon, output.cov = output.cov, output.FIM = output.FIM, output.RE = output.RE, output.SS = output.SS, verbose = verbose)
}
