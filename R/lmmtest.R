#' Testing Fixed Effects and Contrasts of the Fixed Effects
#'
#' @description lmmtest is used to test fixed effects or contrasts of fixed effects by t-statistic.

#' @param fit Output of \code{\link{lmmfit}} or \code{\link{lmm}}, which contains
#' coef (estimates of fixed effects), a matrix with rows representing the fixed effects and columns the different response variables in the model,
#' cov (covariance matrix of the fixed effects), an array of three dimensions for different response variables in the model,
#' df (residual degree of freedom in the linear model).
#' @param index A vector of integers or characters indicating which fixed effects are to be tested. By default index consists of all of the fixed effects. Ignored if contrast is not NULL.
#' @param contrast A matrix with columns corresponding to contrasts of the fixed effects to be tested.
#' @param alternative A character string specifying the alternative hypothesis, one of "two.sided", "greater" or "less".
#'
#' @return A matrix of coefficients, t-values and p-values, in which the rows correspond to the features and the columns the fixed effects (covariates). .
#'
#' @importFrom stats pt
#'
#' @export
lmmtest <- function(fit, index, contrast = NULL, alternative = c("two.sided", "less", "greater")){
alternative <- match.arg(alternative)
	if (is.null(contrast)){
		if (missing(index)) index <- 1:nrow(fit$coef)
		contrast <- diag(nrow(fit$coef))
		colnames(contrast) <- rownames(fit$coef)
		contrast <- contrast[, index, drop = FALSE]
	}

	eff <- t(contrast)%*%fit$coef
	tval <- matrix(nrow = nrow(eff), ncol = ncol(eff), dimnames = dimnames(eff))
	for (j in 1:ncol(fit$coef)){
		tval[, j] <- eff[, j]/sqrt(diag(t(contrast)%*%fit$cov[,,j]%*%contrast))
		}

	df <- fit$df
	if (alternative == "less") {
        pval <- pt(tval, df)
    } else if (alternative == "greater") {
        pval <- pt(tval, df, lower.tail = FALSE)
    } else pval <- 2 * pt(-abs(tval), df)

	rownames(eff) <- paste0(rownames(eff), "_coef")
	rownames(tval) <- paste0(rownames(tval), "_t")
	rownames(pval) <- paste0(rownames(pval), "_p")

cbind(t(eff), t(tval), t(pval))
}
