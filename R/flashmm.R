#' Fitting Linear Mixed-effects Models
#'
#' @description flashmm, a wrapper function of lmm, fits linear mixed-effects models (LMM) by sample-level data. The LMM parameters are estimated by either restricted maximum likelihood (REML) or maximum likelihood (ML) method with Fisher scoring (FS) gradient descent algorithm.
#'
#' @param Y a features-by-samples matrix of responses. For single-cell RNA sequencing (scRNA-seq) data, Y is a genes-by-cells gene expression matrix, whose entries can be raw counts, log-transformed values, or normalized values. If Y is a matrix of raw counts, it will be log-transformed as log2(1 + Y) for differential expression analysis.
#' @param fixed a one-sided model formula used to create the fixed-effects design matrix by model.matrix function.
#' @param random a one-sided formula specifying a model with a single random-effects component or a list of one-sided formulas specifying multiple random-effects components, used to create the random-effects design matrix.
#' @param metadata a data frame containing the variables named in fixed and random.
#' @param is.counts a logical scalar indicating whether Y is a matrix of raw counts or processed values (e.g., log-transformed or normalized). If Y contains raw counts, it will be log-transformed as log2(1 + Y).
#' @param method a character string, either "REML" or "ML". Defaults to "REML". If ‘REML’, the model is fit using REML; otherwise, using ML.
#' @param nBlocks the number of blocks into which large datasets are subdivided to reduce memory usage when computing summary statistics. The default value, nBlocks = ceiling((ncol(Y) * 1e-08) * nrow(Y)), may not be adequate. If encountering the error: vector memory limit reached, you should increase the nBlocks value to avoid the issue.
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
#' @importFrom stats pt model.matrix
#' @import Matrix
#'
#' @seealso \code{\link{lmmfit}}, \code{\link{lmm}}
#' @examples
#' #Generate data
#' set.seed(2024)
#'
#' #Y: gene expression matrix
#' n <- 1e3
#' m <- 10
#' Y <- matrix(rnorm(m*n), m, n)
#' rownames(Y) <- paste0("Gene", 1:nrow(Y))
#'
#' #metadata: samples (sam) and treatments (trt)
#' metadata <- data.frame(trt = character(n), sam = character(n))
#' metadata$trt <- sample(c("A", "B"), n, replace = TRUE)
#' q <- 10
#' metadata$sam <- paste0("S", sample.int(q, n, replace = TRUE))
#'
#' #Model formulas
#' fixed <- ~ 0 + trt
#' random <- ~ 0 + sam
#'
#' #Fit LMM by flashmm with the model formulas based on the cell-level data
#' fit1 <- flashmm(Y, fixed, random, metadata = metadata, is.counts = FALSE)
#'
#' #Fit LMM by lmmfit with the model design matrices based on the cell-level data
#' #design matrices
#' X <- model.matrix(fixed, metadata)
#' Z <- model.matrix(random, metadata)
#' d <- ncol(Z)
#' fit2 <- lmmfit(Y, X, Z, d = d)
#'
#' identical(fit1, fit2)
#'
#' #Fit LMM by lmm based on summary-level data
#' #Compute and store the summary-level data:
#' n <- nrow(X)
#' XX <- t(X)%*%X
#' XY <- t(Y%*%X)
#' ZX <- t(Z)%*%X
#' ZY <- t(Y%*%Z)
#' ZZ <- t(Z)%*%Z
#' Ynorm <- rowSums(Y*Y)
#' fit3 <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d)
#'
#' identical(fit2, fit3)
#'
#' #Hypothesis testing
#' lmmtest(fit1)
#' lmmtest(fit1, index = 2)
#' lmmtest(fit1, contrast = cbind("B-A" = c(-1, 1)))
#'
#' @export
flashmm <- function(Y, fixed, random, metadata, is.counts = FALSE, method = c("REML", "ML"), nBlocks = NULL, theta0 = NULL, max.iter = 50, epsilon = 1e-5, output.cov = TRUE, output.FIM = FALSE, output.RE = FALSE, output.SS = FALSE, verbose = FALSE)
{
stopifnot(!any(is.na(Y)))
if (is.vector(Y)) Y <- t(Y)
stopifnot(ncol(Y) == nrow(metadata))

method <- match.arg(method)

##fixed effects desigm matrix
if (!inherits(fixed, "formula") || length(fixed) != 2) {
	stop("fixed-effects model must be a formula of the form \" ~ covariates\"")
	}
stopifnot(all(all.vars(fixed) %in% colnames(metadata)))

X <- model.matrix(fixed, data = metadata)

##random effects desigm matrix
if (is.list(random)) {
	Z <- NULL
	dcp <- NULL
	for (i in 1:length(random)) {
		if (!inherits(random[[i]], "formula") || length(random[[i]]) != 2) {
			stop("random-effects model must be a formula of the form \" ~ random_variables\"")
		}
		stopifnot(all(all.vars(random[[i]]) %in% colnames(metadata)))
	Z <- cbind(Z, model.matrix(random[[i]], data = metadata))
	dcp <- c(dcp, ncol(Z))
	}
	dcp <- dcp - c(0, dcp[-length(dcp)])
} else {
	if (!inherits(random, "formula") || length(random) != 2) {
		stop("random-effects model must be a formula of the form \" ~ random_variables\"")
	}
	stopifnot(all(all.vars(random) %in% colnames(metadata)))
	Z <- model.matrix(random, data = metadata)
	dcp <- ncol(Z)
}

stopifnot(sum(dcp) == ncol(Z))

##summary statistics
if (is.null(nBlocks)) nBlocks <- ceiling((ncol(Y)*1e-08)*nrow(Y))

nr <- nrow(Y)
if (nBlocks > nr){
	message(paste0("Note: nBlocks=", nBlocks, " by default, changed to nrow(Y), i.e., nBlocks=", nr, "."))
	nBlocks <- nr
}
blsize <- round(nr/nBlocks)
if (nBlocks*blsize < nr) blsize <- round(nr/nBlocks) + 1

XY <- NULL
ZY <- NULL
Ynorm <- NULL
for (i in 1:nBlocks){
  j <- (1+(i-1)*blsize):(min(nr, i*blsize))
  if (is.counts) {
  	Yj <- log2(1 + Y[j, , drop = FALSE])
  	XY <- cbind(XY, t(Yj%*%X))
  	ZY <- cbind(ZY, t(Yj%*%Z))
  	Ynorm <- c(Ynorm, rowSums(Yj*Yj))
  	rm(Yj)
  } else {
  	XY <- cbind(XY, t(Y[j, , drop = FALSE]%*%X))
  	ZY <- cbind(ZY, t(Y[j, , drop = FALSE]%*%Z))
  	Ynorm <- c(Ynorm, rowSums(Y[j, , drop = FALSE]*Y[j, , drop = FALSE]))
  }
}

n <- nrow(X)
XX <- t(X)%*%X
ZX <- t(Z)%*%X
ZZ <- t(Z)%*%Z

rm(X, Y, Z)
gc()

lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = dcp, theta0 = theta0, method = method, max.iter = max.iter, epsilon = epsilon, output.cov = output.cov, output.FIM = output.FIM, output.RE = output.RE, output.SS = output.SS, verbose = verbose)
}
