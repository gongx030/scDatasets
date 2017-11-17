library(GenomicRanges)

exprs <- function(x, ...) UseMethod('exprs', x)
exprs.default <- Biobase::exprs

#' Extract the expression matrix
# 
#' @author Wuming Gong, \email{gongx030@umn.edu}
# 
#' @importFrom GenomicRanges values
exprs.Expression <- function(gene) GenomicRanges::values(gene$gr)[['expression']]

#' Process the scRNA-seq data matrix
preprocess2 <- function(X, min.expressed.gene = 2000, min.expressed.cell = 2, max.expressed.ratio = 0.9, normalize.by.size.effect = TRUE){

	M0 <- ncol(X) 
	N0 <- nrow(X)

	cat(sprintf('[%s] number of input cells(nrow(X))=%.d\n', Sys.time(), N0))
	cat(sprintf('[%s] number of input genes(ncol(X))=%.d\n', Sys.time(), M0))
	X <- X[, Matrix::colSums(X > 1) >= min.expressed.gene, drop = FALSE]
	M <- ncol(X)
	cat(sprintf('[%s] number of input cells that express at least %d genes=%.d\n', Sys.time(), min.expressed.gene, M))
	X <- X[Matrix::rowSums(X > 1) <= max.expressed.ratio * M & Matrix::rowSums(X > 1) >= min.expressed.cell, , drop = FALSE]
	N <- nrow(X)
	cat(sprintf('[%s] number of input genes that are expressed in at least %d cells and at most %.0f%% cells=%.d\n', Sys.time(), min.expressed.cell, max.expressed.ratio * 100, N))
	cat(sprintf('[%s] sparsity of expression matrix=%.1f%%\n', Sys.time(), 100 * (N * M - sum(X > 0)) / (N * M)))

	if (normalize.by.size.effect){
		cat(sprintf('[%s] scaling raw read counts by size factor\n', Sys.time()))
	  sf <- apply((X + 1) / exp(Matrix::rowMeans(log(X + 1))), 2, median)
		X <- t(t(X) / sf)
	}

	as.matrix(X)

} # end of preprocess2

