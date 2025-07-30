#' Fitting Linear Mixed-effects Models
#'
#' @description lmm is used to fit linear mixed-effects models (LMM) based on summary-level data. The LMM parameters are estimated by either restricted maximum likelihood (REML) or maximum likelihood (ML) method with Fisher scoring (FS) gradient descent algorithm.
#'
#' @param XX = t(X)\%*\%X, where X is a design matrix for fixed effects.
#' @param XY = t(Y\%*\%X), where Y is a features-by-samples matrix of observed responses (genes-by-cells expression matrix for scRNA-seq).
#' @param ZX = t(Z)\%*\%X, where Z = [Z1, ..., Zk], a design matrix for k random factors (variables or random components).
#' @param ZY = t(Y\%*\%Z).
#' @param ZZ = t(Z)\%*\%Z.
#' @param Ynorm = rowSums(Y*Y), norms for features (each row in Y).
#' @param n = nrow(X), number of samples (cells in scRNA-seq).
#' @param d = (d1,...,dk), where di = ncol(Zi), number of columns in Zi. sum(d) = ncol(Z), number of columns in Z. For the model with only one random factor, d = ncol(Z).
#' @param theta0 A vector of initial values of the variance components, (s1, ...,sk, s_(k+1)), si = sigma_i^2, the variance component of the i-th type random effects. s_(k+1) = sigma^2, the variance component of model residual error.
#' @param method Either REML or ML with Fisher scoring (FS) iterative algorithm.
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
lmm <- function(XX, XY, ZX, ZY, ZZ, Ynorm, n, d = ncol(ZZ), theta0 = NULL, method = c("REML", "ML"), max.iter = 50, epsilon = 1e-5, output.cov = TRUE, output.RE = FALSE)
{
stopifnot(!any(is.na(XY)), !any(is.na(ZX)), !any(is.na(ZY)))
stopifnot(sum(d) == ncol(ZZ))

method <- match.arg(method)
p <- ncol(ZX)
k <- length(d)

XXinv <- try(chol2inv(chol(XX)), silent = TRUE)
if (inherits(XXinv, "try-error")) {
	stop("XX is not positive-definite or X is not full column rank.")
	}

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

niter <- NULL
dlogL <- NULL
loglike <- NULL
theta <- matrix(nrow = k + 1, ncol = ncol(XY), dimnames = list(paste0("var", c(1:k, 0)), colnames(XY)))
setheta <- theta
RE <- matrix(nrow = nrow(ZY), ncol = ncol(XY), dimnames = dimnames(ZY))
beta <- matrix(nrow = nrow(XY), ncol = ncol(XY), dimnames = dimnames(XY))
sebeta <- beta
covbeta <- array(dim = c(nrow(XY), nrow(XY), ncol(XY)),
	dimnames = list(rownames(XY), rownames(XY), colnames(XY)))

for (jy in 1:ncol(ZY)) {
	if (is.null(theta0)) {
		s <- c(rep(0, k), yry[jy]/(n-p))
	} else s <- theta0

vest <- varest(ZZres, zrz, zryj = zry[, jy], yryj = yry[jy], n = n, d = d, s = s, pres = pres, max.iter = max.iter, epsilon = epsilon)
s <- vest$s
dl <- vest$dl
iter <- vest$iter
Minv <- vest$Minv
logdet <- vest$logdetM0

#if (max(abs(dl)) > epsilon) {
#	warningText <- paste0("The first derivatives of log likelihood for Y", jy)
#	dlText <- paste0(ifelse(abs(dl) > 1e-3,
#	round(dl, 4), format(dl, digits = 3, scientific = TRUE)), collapse = ", ")
#	warning(paste0(warningText, ": ", dlText, ", doesn't reach epsilon ", epsilon))
#	}
#

sr <- s[1:k]/s[k+1]
Dtheta <- sweep(ZZ, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d))
M <- try(solve(Dtheta), silent = TRUE)
if (inherits(M, "try-error")) M <- ginv(Dtheta)
#qrM0 <- qr(M)
#logdet <- sum(log(abs(diag(qrM0$qr))))

M <- sweep(M, 2, STATS = rep(sr, times = d), FUN = "*")
xvx <- XXinv + xxz%*%(ginv(diag(sum(d)) - M%*%(ZX%*%xxz))%*%(M%*%t(xxz)))
xvy <- XY[, jy] - t(ZX)%*%(M%*%ZY[, jy])
b <- xvx%*%xvy
covbeta[,,jy] <- (xvx + t(xvx))*(s[k+1]/2)

RE[, jy] <- M%*%(ZY[, jy] - ZX%*%b)

niter <- c(niter, iter)
theta[, jy] <- s
setheta[, jy] <- sqrt(diag(Minv))
beta[, jy] <- b
dlogL <- cbind(dlogL, dl)
loglike <- c(loglike, -(n-pres)*(1+log(2*pi*s[k+1]))/2 + logdet/2)
sebeta[, jy] <- sqrt(diag(as.matrix(covbeta[,,jy])))
}

tval <- beta/sebeta
pval <- 2 * pt(-abs(tval), df = n-p)

nonconverge <- which(colSums(abs(dlogL) > epsilon) > 0)
if (length(nonconverge) > 0) {
	#warningText <- paste0("The first derivatives of log likelihood for Y", jy)
	#dlText <- paste0(ifelse(abs(dl) > 1e-3,
	#round(dl, 4), format(dl, digits = 3, scientific = TRUE)), collapse = ", ")
	#
	#dlrange <- range(apply(abs(dlogL[, nonconverge, drop = FALSE]), 2, max))
	#dltext <- " features (the rows of Y)
	#for which fitting LMM doesn't converge, i.e., abs(dlogL), "
	#warning(paste0(length(nonconverge), dltext,
	#dlrange[1], " ~ ", dlrange[2], ", > epsilon ", epsilon))
	dltext <- " features (the rows of Y) for which fitting LMM doesn't converge with abs(dlogL) > "
	warning(paste0(length(nonconverge), dltext, "epsilon, ", epsilon))
	}

if (!output.cov) covbeta <- NULL
if (!output.RE) RE <- NULL

list(method = method, dlogL = dlogL, logLik = loglike, niter = niter, coef = beta, se = sebeta, t = tval, p = pval, cov = covbeta, df = n-p, theta = theta, se.theta = setheta, RE = RE)
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
#' @param p Equal to 0 for ML, otherwise, ncol(X) for REML.
#' @param max.iter The maximal number of iterations.
#' @param epsilon Positive convergence tolerance.
#'
#' @return A list consisting of
#' estimates of variance components (s),
#' first partial derivatives of log-likehood (dl),
#' number of iterations (iter), and
#' inverse of Fisher information matrix (Minv).
#'
#' @importFrom MASS ginv
#'
#' @keywords internal
#'
#' @noRd
varest <- function(ZZres, zrz, zryj, yryj, n, d, s, pres, max.iter = 50, epsilon = 1e-5)
{
  k <- length(d)

  dl <- 100
  iter <- 0
  while ((max(abs(dl)) > epsilon)	& (iter < max.iter)){
    iter <- iter + 1

    fs <- matrix(NA, k+1, k+1)
    dl <- rep(NA, k+1)

    sr <- s[1:k]/s[k+1]
    M <- solve(sweep(zrz, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d)))
    yRZ <- t(zryj)%*%M

    if (pres > 0){
    	M0 <- M
    } else {
    	M0inv <- sweep(ZZres, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d))
    	M0 <- try(solve(M0inv), silent = TRUE)
    	if (inherits(M0, "try-error")) M0 <- ginv(M0inv)
    }
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
  }

qrM0 <- qr(M0)
logdetM0 <- sum(log(abs(diag(qrM0$qr))))
list(s = c(s), dl = dl, iter = iter, Minv = Minv, logdetM0 = logdetM0)
}
