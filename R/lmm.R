#' Fitting Linear Mixed-effects Models
#'
#' @description lmm is used to fit a linear mixed-effects model (LMM) based on summary-level data. The LMM parameters are estimated by either restricted maximum likelihood (REML) or maximum likelihood (ML) method with Fisher scoring (FS) gradient descent algorithm.
#'
#' @param XX = t(X)\%*\%X, where X is a design matrix for fixed effects.
#' @param XY = t(Y\%*\%X), where Y is a features-by-samples matrix of observed responses (genes-by-cells expression matrix for scRNA-seq).
#' @param ZX = t(Z)\%*\%X, where Z = [Z1, ..., Zk] is a design matrix for k random components (factors).
#' @param ZY = t(Y\%*\%Z).
#' @param ZZ = t(Z)\%*\%Z.
#' @param Ynorm = rowSums(Y*Y), norms for features (each row).
#' @param n = nrow(X), number of samples (cells in scRNA-seq).
#' @param d = (d1,...,dk), where di = ncol(Zi), the number of columns in Zi. sum(d) = ncol(Z), the number of columns in Z. For the model with only one random component, d = ncol(Z).
#' @param method a character string, either "REML" or "ML". Defaults to "REML". If ‘REML’, the model is fit using REML; otherwise, using ML.
#' @param theta0 a vector of initial values of the variance components, (s1, ...,sk, s_(k+1)), si = sigma_i^2, the i-th variance component. s_(k+1) = sigma^2, the variance component of model residual error.
#' @param max.iter the maximal number of iterations for the iterative algorithm.
#' @param epsilon positive convergence tolerance. If the absolute value of the first partial derivative of log likelihood is less than epsilon, the iterations converge.
#' @param output.cov a logical scalar. If TRUE, output an array of the covariance matrices for the estimated coefficients, cov.beta.
#' @param output.FIM a logical scalar. If TRUE, output an array of the Fisher information matrices (FIM) for fixed effects, FIM.beta[,,i] = \eqn{X^T(Cov(Y[i,])^{-1}X} for the i-th feature, and variance components, FIM.theta[,,i] = \eqn{Var(\partial(logL_i)/\partial(theta)) = - E[\partial^2(logL_i)/{\partial(theta)\partial(theta^T)}]} for the i-th feature.
#' @param output.RE a logical scalar. If TRUE, output the best linear unbiased prediction (BLUP) of the random effects.
#' @param output.SS a logical scalar. If TURE, output list(XX, XY, ZX, ZY, ZZ, Ynorm, n), a list of the summary-level data (summary statistics), which can also be computed using the sslmm function.
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
#' #Fit LMM by summary-level data
#' #Compute and store the summary-level data:
#' n <- nrow(X)
#' XX <- t(X)%*%X
#' XY <- t(Y%*%X)
#' ZX <- t(Z)%*%X
#' ZY <- t(Y%*%Z)
#' ZZ <- t(Z)%*%Z
#' Ynorm <- rowSums(Y*Y)
#' fit <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d)
#' str(fit)
#'
#' @export
lmm <- function(XX, XY, ZX, ZY, ZZ, Ynorm, n, d = ncol(ZZ), method = c("REML", "ML"), theta0 = NULL, max.iter = 50, epsilon = 1e-5, output.cov = TRUE, output.FIM = FALSE, output.RE = FALSE, output.SS = FALSE, verbose = FALSE)
{
stopifnot(!any(is.na(XY)), !any(is.na(ZX)), !any(is.na(ZY)))
stopifnot(sum(d) == ncol(ZZ))

if (any(Ynorm <= 0)) {
	warning(sum(Ynorm <= 0), " feature(s) with only zero counts!")
}

method <- match.arg(method)
p <- ncol(ZX)
k <- length(d)

XXinv <- try(chol2inv(chol(XX)), silent = TRUE)
if (inherits(XXinv, "try-error")) {
	stop("XX is not positive-definite or X is not full column rank.")
	}

##zrz:= Z^TRZ, zry:= Z^TRy, yry:= [y^TRy],
##where R = I_n - X(X^TX)^{-1}X^T.
xxz <- XXinv%*%t(ZX)
zrz <- ZZ - ZX%*%(XXinv%*%t(ZX))
zry <- ZY - ZX%*%(XXinv%*%XY)
yry <- Ynorm - colSums(XY*(XXinv%*%XY))

if (method == "REML"){
	pres <- p
	ZZres <- zrz
} else {
	pres <- 0
	ZZres <- ZZ
}

niter <- rep(NA, ncol(XY))
names(niter) <- colnames(XY)
loglike <- niter
theta <- matrix(nrow = k + 1, ncol = ncol(XY), dimnames = list(paste0("var", c(1:k, 0)), colnames(XY)))
setheta <- theta
dlogL <- theta
beta <- matrix(nrow = nrow(XY), ncol = ncol(XY), dimnames = dimnames(XY))
sebeta <- beta
if (output.cov){
	covbeta <- array(dim = c(nrow(XY), nrow(XY), ncol(XY)),
		dimnames = list(rownames(XY), rownames(XY), colnames(XY)))
	} else covbeta <- NULL
if (output.FIM){
	FIMbeta <- array(dim = c(nrow(XY), nrow(XY), ncol(XY)),
		dimnames = list(rownames(XY), rownames(XY), colnames(XY)))
	FIMtheta <- array(dim = c(k+1, k+1, ncol(XY)),
		dimnames = list(rownames(theta), rownames(theta), colnames(XY)))
	} else {
	FIMbeta <- NULL
	FIMtheta <- NULL
	}
if (output.RE){
	RE <- matrix(nrow = nrow(ZY), ncol = ncol(XY), dimnames = dimnames(ZY))
	} else RE <- NULL

#for (jy in 1:ncol(ZY)) {
for (jy in which(Ynorm > 0)) {
	if (is.null(theta0)) {
		s <- c(rep(0, k), yry[jy]/(n-p))
	} else s <- theta0

	vest <- varest(ZZres, zrz, zryj = zry[, jy], yryj = yry[jy], n = n, d = d, s = s, pres = pres, max.iter = max.iter, epsilon = epsilon)

niter[jy] <- vest$iter
dlogL[, jy] <- vest$dl
#dl <- vest$dl
s <- vest$s
Minv <- vest$Minv
logdet <- vest$logdetM0
if (output.FIM) FIMtheta[,,jy] <- vest$FIM

#if (max(abs(dl)) > epsilon) {
#	warningText <- paste0("The first derivatives of log likelihood for Y", jy)
#	dlText <- paste0(ifelse(abs(dl) > 1e-3,
#	round(dl, 4), format(dl, digits = 3, scientific = TRUE)), collapse = ", ")
#	warning(paste0(warningText, ": ", dlText, ", doesn't reach epsilon ", epsilon))
#	}
#

##beta = (X^TV^{-1}X)^{-1}(X^TV^{-1}y)
##xvx:= (X^TV^{-1}X)^{-1} = [X^TX - (X^TZ)M0D(Z^TX)]^{-1}
##Use Sherman-Morrison-Woodbury formula if p >> q:
##xvx:= (X^TX)^{-1} + (X^TX)^{-1}X^TZ[I_q - M0DZ^TX(X^TX)^{-1}X^TZ]^{-1}M0DZ^TX(X^TX)^{-1}
##xvy:= X^TV^{-1}y = X^Ty - (X^TZ)M0D(Z^Ty)
##
##M:= M0 = (I_q + DZ^TZ)^{-1}
sr <- s[1:k]/s[k+1]
DZZ1 <- sweep(ZZ, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d))
M <- try(solve(DZZ1), silent = TRUE)
if (inherits(M, "try-error")) M <- ginv(DZZ1)
#qrM0 <- qr(M)
#logdet <- sum(log(abs(diag(qrM0$qr))))
##
##M0*D
M <- sweep(M, 2, STATS = rep(sr, times = d), FUN = "*")
##Use Sherman-Morrison-Woodbury formula if p >> q
##xvx <- XXinv + xxz%*%(ginv(diag(sum(d)) - M%*%(ZX%*%xxz))%*%(M%*%t(xxz)))
xvx <- XX - t(ZX)%*%M%*%ZX
if (output.FIM) FIMbeta[,,jy] <- (xvx + t(xvx))/(2*s[k+1])
#xvx <- ginv(XX - t(ZX)%*%M%*%ZX)
xvx <- ginv(xvx)
#if (all(diag(as.matrix(xvx)) >= 0))
	xvy <- XY[, jy] - t(ZX)%*%(M%*%ZY[, jy])
	b <- xvx%*%xvy
	if (output.cov) covbeta[,,jy] <- (xvx + t(xvx))*(s[k+1]/2)
	if (output.RE) RE[, jy] <- M%*%(ZY[, jy] - ZX%*%b)
	theta[, jy] <- s
	setheta[, jy] <- sqrt(diag(Minv))
	beta[, jy] <- b
	loglike[jy] <- -(n-pres)*(1+log(2*pi*s[k+1]))/2 + logdet/2
	sebeta[, jy] <- sqrt(diag(as.matrix(covbeta[,,jy])))
}

tval <- beta/sebeta
pval <- 2 * pt(-abs(tval), df = n-p)

nonconverge <- which(colSums(abs(dlogL) > epsilon) > 0)
if (verbose) {
	#warningText <- paste0("The first derivatives of log likelihood for Y", jy)
	#dlText <- paste0(ifelse(abs(dl) > 1e-3,
	#round(dl, 4), format(dl, digits = 3, scientific = TRUE)), collapse = ", ")
	#
	#dlrange <- range(apply(abs(dlogL[, nonconverge, drop = FALSE]), 2, max))
	#dltext <- " features (the rows of Y)
	#for which fitting LMM doesn't converge, i.e., abs(dlogL), "
	#warning(paste0(length(nonconverge), dltext,
	#dlrange[1], " ~ ", dlrange[2], ", > epsilon ", epsilon))
	dltext <- " feature(s) for which LMM fitting do not converge, i.e., abs(dlogL) > "
	message(paste0(length(nonconverge), dltext, "epsilon=", epsilon, "."))
	}

##summary-level data (summary statistics)
summary_data <- NULL
if (output.SS) summary_data <- list(n=n, XX=XX, XY=XY, ZX=ZX, ZY=ZY, ZZ=ZZ, Ynorm = Ynorm)

list(method = method, epsilon = epsilon, dlogL = dlogL, logLik = loglike, niter = niter, coef = beta, se = sebeta, t = tval, p = pval, df = n-p, cov = covbeta, theta = theta, se.theta = setheta, FIM.beta = FIMbeta, FIM.theta = FIMtheta, RE = RE, summary_data = summary_data)
}



#' A internal function to estimate variance components for one feature (gene).
#'
#' @description This function is used internally (inside lmm).
#'
#' @param ZZres Equal to t(Z)\%*\%Z for ML, otherwise, t(Z)\%*\%Z - ZX\%*\%(XXinv\%*\%t(ZX)) for REML.
#' @param zrz t(Z)\%*\%Z - ZX\%*\%(XXinv\%*\%t(ZX)).
#' @param zryj zry[, j], where zry = ZY - ZX\%*\%(XXinv\%*\%XY)
#' @param yryj yry[j], where yry = Ynorm - colSums(XY*(XXinv\%*\%XY))
#' @param n Numbers of samples (cells in scRNA-seq).
#' @param d A vector of (m1,...,mk), mi = ncol(Zi), number of columns in Zi.
#' @param s A vector of initial values of the variance components, (s1, ...,sk, s_(k+1)).
#' @param pres Equal to 0 for ML, otherwise, ncol(X) for REML.
#' @param max.iter The maximal number of iterations.
#' @param epsilon Positive convergence tolerance.
#'
#' @return A list consisting of
#' estimates of variance components (s),
#' first partial derivatives of log-likehood (dl),
#' number of iterations (iter),
#' inverse of Fisher information matrix (Minv), and
#' Fisher information matrix (FIM) for the variance components.
#'
#' @importFrom MASS ginv
#'
#' @keywords internal
#'
#' @noRd
varest <- function(ZZres, zrz, zryj, yryj, n, d, s, pres, max.iter = 50, epsilon = 1e-5)
{
k <- length(d)
sr <- s[1:k]/s[k+1]
M <- solve(sweep(zrz, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d)))
if (pres > 0){
	M0 <- M
} else {
	M0inv <- sweep(ZZres, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d))
	M0 <- try(solve(M0inv), silent = TRUE)
	if (inherits(M0, "try-error")) M0 <- ginv(M0inv)
    }

dl <- 100
iter <- 0
while ((max(abs(dl)) > epsilon)	& (iter < max.iter)){
	iter <- iter + 1
    fs <- matrix(NA, k+1, k+1)
    dl <- rep(NA, k+1)

    yRZ <- t(zryj)%*%M
    ZVZ <- ZZres%*%M0
    ZV2Z <- ZVZ%*%M0

    mi <- 0
    for (i in 1:k){
    	ik <- (mi+1):(mi+d[i])
    	dl[i] <- (sum((yRZ[ik])^2)/s[k+1]^2 - sum(diag(ZVZ[ik, ik, drop = FALSE]))/s[k+1])/2
    	mj <- 0
    	for (j in 1:i){
    		ji <- (mj+1):(mj+d[j])
    		fs[i, j] <- sum((ZVZ[ji, ik])^2)/s[k+1]^2/2
    		fs[j, i] <- fs[i, j]
    		mj <- mj + d[j]
    	}
    	j <- k+1
    	fs[i, j] <- sum(diag(ZV2Z[ik, ik, drop = FALSE]))/s[k+1]^2/2
    	fs[j, i] <- fs[i, j]
    	mi <- mi + d[i]
    }

    i <- k+1
    fs[i, i] <- (n - pres - sum(d) + sum(t(M0)*M0))/s[k+1]^2/2
    yR2y <- yryj - sum(((t(M) + diag(sum(d)))%*%zryj)*(M%*%(rep(sr, times = d)*zryj)))
    dl[i] <-  (yR2y/s[k+1]^2 - (n - pres - sum(d) + sum(diag(M0)))/s[k+1])/2

    Minv <- ginv(fs)
    s <- s + Minv%*%dl

    sr <- s[1:k]/s[k+1]
    M <- solve(sweep(zrz, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d)))
    if (pres > 0){
    	M0 <- M
    } else {
    	M0inv <- sweep(ZZres, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d))
    	M0 <- try(solve(M0inv), silent = TRUE)
    	if (inherits(M0, "try-error")) M0 <- ginv(M0inv)
    }
  }

##Recompute sigma^2, i.e., s[K+1].
##s[K+1] = y^T(R_gamma)y/n > 0 for ML
##s[K+1] = y^T(R_gamma)y/(n-p) > 0 for REML
##y^T(R_gamma)y = y^TRy - (y^TRZ)MD(Z^TRy)
##M:=M*D, D = diag(s[i]/s[k+1])
M <- sweep(M, 2, STATS = rep(sr, times = d), FUN = "*")
s[k+1] <- (yryj - sum(zryj*(M%*%zryj)))/(n - pres)
M <- NULL

##update fs and dl using new s.
    fs <- matrix(NA, k+1, k+1)
    dl <- rep(NA, k+1)

    sr <- s[1:k]/s[k+1]
    M <- solve(sweep(zrz, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d)))
    if (pres > 0){
    	M0 <- M
    } else {
    	M0inv <- sweep(ZZres, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d))
    	M0 <- try(solve(M0inv), silent = TRUE)
    	if (inherits(M0, "try-error")) M0 <- ginv(M0inv)
    }

    yRZ <- t(zryj)%*%M
    ZVZ <- ZZres%*%M0
    ZV2Z <- ZVZ%*%M0

    mi <- 0
    for (i in 1:k){
    	ik <- (mi+1):(mi+d[i])
    	dl[i] <- (sum((yRZ[ik])^2)/s[k+1]^2 - sum(diag(ZVZ[ik, ik, drop = FALSE]))/s[k+1])/2
    	mj <- 0
    	for (j in 1:i){
    		ji <- (mj+1):(mj+d[j])
    		fs[i, j] <- sum((ZVZ[ji, ik])^2)/s[k+1]^2/2
    		fs[j, i] <- fs[i, j]
    		mj <- mj + d[j]
    	}
    	j <- k+1
    	fs[i, j] <- sum(diag(ZV2Z[ik, ik, drop = FALSE]))/s[k+1]^2/2
    	fs[j, i] <- fs[i, j]
    	mi <- mi + d[i]
    }

    i <- k+1
    fs[i, i] <- (n - pres - sum(d) + sum(t(M0)*M0))/s[k+1]^2/2
    yR2y <- yryj - sum(((t(M) + diag(sum(d)))%*%zryj)*(M%*%(rep(sr, times = d)*zryj)))
    dl[i] <-  (yR2y/s[k+1]^2 - (n - pres - sum(d) + sum(diag(M0)))/s[k+1])/2

    Minv <- ginv(fs)

##log(det(M0))
qrM0 <- qr(M0)
logdetM0 <- sum(log(abs(diag(qrM0$qr))))

list(s = c(s), dl = dl, iter = iter, Minv = Minv, FIM = fs, logdetM0 = logdetM0)
}
