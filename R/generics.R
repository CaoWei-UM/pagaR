#' Calculate PAGA connectivities
#'
#' Computes the PAGA connectivities for a given KNN adjacency matrix of cells.
#' Calculateconnectivities between every cell groups and its nearest neighbors.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return This function can return a adjacency matrix with the paga
#' connectivities information of given cell groups . When running on
#' a \code{\link{Seurat}} object, this function returns the
#' \code{\link{Seurat}} object with the Graphs stored in their respective
#' slots. Names of the Graph can be found with \code{\link{Graphs}}.
#'
#' @examples
#' library(igraph)
#' pbmc<-CalculatePAGA(object=pbmc, group.by='celltype')
#' plot(graph_from_adjacency_matrix(pbmc@graphs$RNA_paga))
#'
#' @rdname CalculatePAGA
#' @export CalculatePAGA
#'
CalculatePAGA <- function(object, ...) {
  UseMethod(generic = 'CalculatePAGA', object = object)
}
