
#' Plot PAGA with cell embedding
#'
#' Plot the PAGA connectivities from connectivities matrix of cell groups.
#' If parameter show.cell=F, only shows cell groups, each cell group is a
#' point, point size indicate cell count in group, point position is the
#' center of the position of cells based on the cell embeddings determined
#' by the reduction technique. otherwise a 2D scatter plot where each point
#' is a cell with same embeddings with cell groups will be put under cell
#' group points. Cells are colored by their group (can be changed with the
#' feature parameter).
#'
#' @param object Seurat object
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param threshold Threshold for connectivities between cell groups, default is 0.01
#' @param cell.plot Weature to plot cell points, default is not.
#' @param assay Assay used to calculate cell group connectivities
#' @param slot Which slot to pull expression data from? Default is data slot
#' @param feature character of a feature to plot. Feature can come from:
#' \itemize{
#'     \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
#'     \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
#'     \item A column name from a \code{DimReduc} object corresponding to the cell embedding values
#'     (e.g. the PC 1 scores - "PC_1")
#' }
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for feature
#' @param feature.cols The multiple colors to form the gradient over. Provide as string vector with
#' the first color corresponding to low values, the second to high.
#' @param pt.size Adjust cell point size for plotting, default is 1.
#' @param group.cols Vector of colors, each color corresponds to an identity class. By default, hue_pal is used.
#' @param label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param raster.dpi Cells convert points to raster format, Pixel resolution for rasterized plot.
#'
#' @return A ggplot object with PAGA and a feature in reduction.
#'
#' @importFrom SeuratObject DefaultAssay Idents DefaultDimReduc Key FetchData
#' @importFrom reshape2 melt
#' @importFrom ggplot2 geom_segment scale_linewidth geom_point scale_size_continuous scale_color_manual scale_fill_manual scale_color_gradientn ggtitle theme_bw
#' @importFrom ggrastr geom_point_rast
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'
#' @examples
#' PlotPAGA(object=pbmc)
#' PlotPAGA(object=pbmc, feature='PTPRC')
#'
#'
PlotPAGA<-function(object,
                   group.by=NULL,
                   reduction=NULL,
                   dims=c(1, 2),
                   threshold=0.01,
                   cell.plot=F,
                   assay=NULL,
                   slot='data',
                   feature=NULL,
                   min.cutoff=NULL,
                   max.cutoff=NULL,
                   feature.cols=c('grey','red'),
                   pt.size=1,
                   group.cols=NULL,
                   label=F,
                   label.size=5,
                   raster.dpi=500){
  # Get data to use
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  graph.name<-paste(assay,group.by,'paga',sep='_')
  group.by <- group.by %||% 'ident'
  object[['ident']] <- Idents(object=object)
  group<-factor(object@meta.data[[group.by]])
  if(!(graph.name %in% names(object@graphs))){stop(paste(graph.name, 'not found'))}
  reduction <- reduction %||% DefaultDimReduc(object = object)
  # Get cell plotting data
  if(length(x=dims)!=2 || !is.numeric(x=dims)){stop("'dims' must be a two-length integer vector")}
  dims <- paste0(Key(object=object[[reduction]]), dims)
  data <- FetchData(object=object, dims)
  # Get group plotting data
  group.data<-as.data.frame(apply(data,2,function(x){tapply(x,group,mean)}))
  colnames(group.data)<-dims
  group.data$group <- levels(group)
  group.data$weight <- as.integer(table(group))
  if(is.null(group.cols)){group.cols<-setNames(scales::hue_pal()(length(levels(group))),levels(group))}
  # Get connect plotting data
  connect.data <- as.matrix(object[[graph.name]])
  connect.data[connect.data < threshold] <- 0
  connect.data[upper.tri(connect.data,diag=T)] <- NA
  connect.data <- reshape2::melt(connect.data, na.rm=T)
  connect.data<-cbind.data.frame(group.data[match(connect.data$Var1,group.data$group),dims],
                   group.data[match(connect.data$Var2,group.data$group),dims],
                   connect.data$value)
  colnames(connect.data)<-c('xstart', 'ystart', 'xend', 'yend', 'connectivities')
  #plot
  if((!cell.plot) & is.null(feature)){
    plot <- ggplot() +
      geom_segment(data=connect.data, aes(x=xstart, y=ystart, xend=xend, yend=yend, linewidth=connectivities), color='black') +
      scale_linewidth(range = c(0.1, 4)) +
      geom_point(data=group.data, aes_string(x=dims[1], y=dims[2], color='group', size='weight')) +
      scale_size_continuous(range = c(4, 10))+
      scale_color_manual(values=group.cols)
  }else if(cell.plot & is.null(feature)){
    data[[group.by]]<-group
    plot <- ggplot() +
      ggrastr::geom_point_rast(data=data, aes_string(x=dims[1], y=dims[2], color=group.by), size=pt.size, raster.dpi=raster.dpi) +
      scale_color_manual(values=group.cols)+
      geom_segment(data=connect.data, aes(x=xstart, y=ystart, xend=xend, yend=yend, linewidth=connectivities), color='black') +
      scale_linewidth(range = c(0.1, 4))+
      geom_point(data=group.data, aes_string(x=dims[1], y=dims[2], size='weight'), color='black') +
      scale_size_continuous(range = c(4, 10))
  }else if(cell.plot & (!is.null(feature))){
    data.feature <- as.numeric(x=FetchData(object=object, vars=feature, slot=slot)[[feature]])
    # Determine cutoffs
    if(is.null(min.cutoff)){min.cutoff<-min(data.feature)}
    if(is.null(max.cutoff)){max.cutoff<-max(data.feature)}
    # Apply cutoffs
    data.feature[data.feature < min.cutoff] <- min.cutoff
    data.feature[data.feature > max.cutoff] <- max.cutoff
    data[[feature]]<-data.feature
    plot <- ggplot() +
      ggrastr::geom_point_rast(data=data, aes_string(x=dims[1], y=dims[2], color=feature), size=pt.size, raster.dpi=raster.dpi) +
      scale_color_gradientn(colours=feature.cols) +
      geom_segment(data=connect.data, aes(x=xstart, y=ystart, xend=xend, yend=yend, linewidth=connectivities), color='black') +
      scale_linewidth(range = c(0.1, 4)) +
      geom_point(data=group.data, aes_string(x=dims[1], y=dims[2], fill='group', size='weight'),color='black', shape = 21) +
      scale_size_continuous(range = c(4, 10)) +
      scale_fill_manual(values=group.cols) +
      ggtitle(feature)
  }else{stop('cell.plot must be TRUE if you want to show feature')}
  if(label){
    plot <- plot + ggrepel::geom_text_repel(data=group.data,aes_string(x=dims[1], y=dims[2],label='group'), size=label.size)
  }
  plot <- plot + theme_bw()
  return(plot)
}
