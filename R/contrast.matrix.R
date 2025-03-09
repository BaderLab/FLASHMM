#' Construct Contrast Matrix
#'
#' @description Construct the contrast matrix to make various comparsions of different treatments.

#' @param contrast A vector of character strings specifying the various comparisons, which are the expressions constituted by model.matrix.names.
#' @param model.matrix.names Column names of model (design) matrix.
#'
#' @return Matrix which columns correspond to contrasts.
#'
#' @examples
#' model_variables <- c("A", "B", "C", "D")
#' contrast <- c("AvsB" = "A-B", "AvsC" = "A-C", 'AvsB.C.D'= "A-(B+C+D)/3")
#' contrast.matrix(contrast, model_variables)
#'
#' @export
contrast.matrix <- function(contrast, model.matrix.names){
	stopifnot(!any(duplicated(model.matrix.names)))

	cm <- matrix(nrow = length(model.matrix.names), ncol = length(contrast))
	rownames(cm) <- model.matrix.names
	colnames(cm) <- contrast
	nmct <- names(contrast)
	if (!is.null(nmct)){
		k <- (nmct != "") & (!is.na(nmct))
		colnames(cm)[k] <- nmct[k]
	}

	covars <- model.matrix.names
	#Special variables with bracket characters
	index <- union(grep("\\(", covars), grep("\\)", covars))
	if (length(index) > 0){
		for (i in index){
			xi <- paste0("SVar", i)
			vi <- gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", covars[i]))
			contrast <- gsub(vi, xi, contrast)
			covars[i] <- xi
		}
	}

	#With colon punctuation mark
	covars <- gsub(":", "cl", covars)
	contrast <- gsub(":", "cl", contrast)

	mdcovars <- diag(length(covars))
	colnames(mdcovars) <- covars
	mdcovars <- as.data.frame(mdcovars)
	for (i in 1:length(contrast)){
		cm[, i] <- with(data = mdcovars, eval(parse(text = contrast[i])))
	}

return(cm)
}
