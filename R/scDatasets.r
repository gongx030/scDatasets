library(SummarizedExperiment)

#' scDatasets: Collections of processed single cell sequencing datasets
#'
#' Collections of processed single cell sequencing datasets
#'
#' @import SummarizedExperiment
#' @docType package
#' @name scDatasets
NULL
# > NULL

#' guo
#'
#' Resolution of Cell Fate Decisions Revealed by Single-Cell Gene Expression Analysis 
#' from Zygote to Blastocyst
#' @docType data
#' @keywords datasets
#' @name guo
#' @usage data(guo)
#' @examples
#' data(guo)
NULL

#' loh
#'
#' 
#' Mapping the Pairwise Choices Leading from Pluripotency to Human Bone, Heart, and Other Mesoderm Cell Types
#' @docType data
#' @keywords datasets
#' @name loh 
#' @usage data(loh)
#' @examples
#' data(loh)
NULL

#' trapnell
#'
#'
#' Single cell RNA-seq of human myoblast differentiation
#' @docType data
#' @keywords datasets
#' @name trapnell
#' @usage data(trapnell)
#' @examples
#' data(trapnell)
NULL

#' guo2
#'
#'
#' Single cell RNA-seq of the development of human primordial germ cells (PGC) and neighboring somatic cells from weeks 4 to 19 post-gestation
#' @docType data
#' @keywords datasets
#' @name guo2
#' @usage data(guo2)
#' @examples
#' data(guo2)
NULL

#' pollen
#'
#'
#' Low-coverage single-cell mRNA sequencing reveals cellular heterogeneity 
#' and activated signaling pathways in developing cerebral cortex
#' @docType data
#' @keywords datasets
#' @name pollen
#' @usage data(pollen)
#' @examples
#' data(pollen)
NULL

#' Process the scRNA-seq data matrix
#'
#' @export
#'
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

