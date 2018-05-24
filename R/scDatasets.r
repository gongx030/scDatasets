#' scDatasets: Collections of processed single cell sequencing datasets
#'
#' Collections of processed single cell sequencing datasets
#'
#' @import SummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @docType package
#' @name scDatasets
NULL

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

#' usoskin
#'
#'
#' RNA-Seq of single cells from the mouse lumbar dorsal root ganglion
#' @docType data
#' @keywords datasets
#' @name usoskin
#' @usage data(usoskin)
#' @examples
#' data(usoskin)
NULL

#' blakeley
#'
#' Single-Cell RNA-seq Defines the Three Cell Lineages of the Human Blastocyst
#' 
#' @docType data
#' @keywords datasets
#' @name usoskin
#' @usage data(blakeley)
#' @examples
#' data(blakeley)
NULL


#' preprocess
#'
#' A function for preprocessing gene expression matrix.
#'
#' @param X Gene expression matrix (Gene by Cell). 
#' @param min.expressed.gene Cell level filtering criteria. For a given cell, if the number of expressed genes are less than min.expressed.gene, we filter it out.  
#' @param min.expressed.cell Gene level filtering criteria. For a given gene, if the number of expressed cells are less than min.expressed.cell, we filter it out.  
#' @param max.expressed.ratio Gene level filtering criteria. For a given gene, if the ratio of expressed cells are larger than max.expressed.ratio, we filter it out.
#' @param normalize.by.size.effect Normaize using size factor.
#' @return Filtered gene expression matrix
#'
#' @export
#'
#' @author Wuming Gong and Il-Youp Kwak
#'
#' @references
#' Wuming Gong, Il-Youp Kwak, Pruthvi Pota, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
#' 
#' @examples
#'
#' library(scDatasets)
#' library(SummarizedExperiment)
#' data(usoskin)
#' X <- assays(usoskin)$count
#' X <- preprocess(X, min.expressed.gene = 0)
#' 
#' @seealso \code{\link{DrImpute}}
preprocess <- function(X, min.expressed.gene = 2000, min.expressed.cell = 2, max.expressed.ratio = 1, normalize.by.size.effect = FALSE){

	M0 <- ncol(X) 
	N0 <- nrow(X)

	cat(sprintf('[%s] number of input genes(nrow(X))=%.d\n', Sys.time(), N0))
	cat(sprintf('[%s] number of input cells(ncol(X))=%.d\n', Sys.time(), M0))
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

} # end of preprocess

