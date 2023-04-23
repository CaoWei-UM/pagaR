#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param group A vector of the groups of each cell that specifies
#' the mapping. Its elements correspond to the KNN matrix.
#'
#' @importFrom igraph graph_from_adjacency_matrix contract simplify as.undirected strength E V
#'
#' @rdname CalculatePAGA
#' @export
#' @method CalculatePAGA default
#'
CalculatePAGA.default <- function(object, group){
  group <- factor(group)
  g <- graph_from_adjacency_matrix(object, mode="directed", weighted=T)
  vc <- contract(g,group)
  V(vc)$name <- levels(group)
  vc <- simplify(vc,remove.loops = F,edge.attr.comb=list(weight="sum"))
  es_inner_cluster <- E(vc)$weight[which_loop(vc)]
  es = strength(vc,mode = 'out',loops = T)
  inter_es = as.undirected(simplify(vc, remove.loops=T), mode='collapse', edge.attr.comb=list(weight="sum"))
  ns = table(group)
  v1 = unlist(head_of(inter_es,E(inter_es))$name)
  v2 = unlist(tail_of(inter_es,E(inter_es))$name)
  expected_random_null <- (es[v1] * ns[v2] + es[v2] * ns[v1])/ (length(group) - 1)
  scaled_value <- E(inter_es)$weight/expected_random_null
  scaled_value[which(expected_random_null==0|scaled_value>1)] <- 1
  E(inter_es)$connectivities <- scaled_value
  mat <- as_adjacency_matrix(inter_es, type="both", attr='connectivities')
  return(mat)
}

#' @param assay Assay used in calculate KNN matrix
#' @param group.by Name of one or more metadata columns to group
#' cells by (for example, orig.ident)
#' @param graph.name Name of graph to use for the paga algorithm
#'
#' @importFrom Seurat DefaultAssay Idents LogSeuratCommand
#'
#' @rdname CalculatePAGA
#' @export
#' @method CalculatePAGA Seurat
#'
CalculatePAGA.Seurat <- function(object, group.by = NULL, assay = NULL, graph.name=NULL){
  assay <- assay %||% DefaultAssay(object = object)
  graph.name <- graph.name %||% paste0(assay, "_nn")
  if (!graph.name %in% names(x = object)) {
    stop(paste0("Provided graph.name: ",graph.name," not present in Seurat object"))
  }
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  paga.results <- CalculatePAGA(object=object[[graph.name]], group=object@meta.data[[group.by]])
  paga.results <- as.Graph(paga.results)
  DefaultAssay(object=paga.results) <- assay
  object@graphs[[paste(assay, group.by, 'paga', sep='_')]] <- paga.results
  object <- LogSeuratCommand(object = object)
  return(object)
}

