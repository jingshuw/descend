#' A subset of the scRNA-seq data from Zeisel et.al. (2015)
#'
#' This dataset is used to demostrate the usage of \code{descend}. The original data contains 3005 cells from 7 cell types and 12234 genes. This dataset contains 1464 cells from 3 of the cell types. The elements in the datasets are as follows:
#'
#' \itemize{
#' \item{\code{count.matrix.small}}{a smaller count matrix containing 100 genes}
#' \item{\code{count.matrix.large}}{a larger count matrix containing 1000 genes}
#' \item{\code{ercc.matrix}}{the UMI counts of the ERCC spike-in genes. Only 57 ERCC spike-ins are kept.}
#' \item{\code{library.size}}{the library sizes of each cell}
#' \item{\code{trueMol}}{the true number of molecules of the ERCC spike-ins}
#' \item{\code{cell.size}}{the cell sizes of each cell. Calculated as \code{library.size/cell efficiency}}
#' \item{\code{labels}}{the cell type labels of each cell}
#' }
#'
#' @name zeisel
#' @docType data
#' @keywords data
NULL
