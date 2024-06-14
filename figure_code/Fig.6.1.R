library(igraph)
library(tidyverse)
library(tidyr)

# Load data
count_net <- read.csv("XXX/source_target_cellphonedb.csv")

# Data preprocessing
count_net <- count_net[, -1]
colnames(count_net) <- c("SOURCE", "TARGET", "count")
count_net <- count_net[!grepl('Plasma', count_net$SOURCE) & !grepl('Plasma', count_net$TARGET), ]
count_net <- count_net[!grepl('Cycling', count_net$SOURCE) & !grepl('Cycling', count_net$TARGET), ]

# Transform data
count_inter <- spread(count_net, TARGET, count)
rownames(count_inter) <- count_inter$SOURCE
count_inter <- count_inter[, -1]
count_inter <- as.matrix(count_inter)

# Visualize network
netVisual_circle(count_inter, weight.scale = TRUE)

# Visualize individual matrices
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Define color palette function
scPalette <- function(n) {
  colorSpace <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#F29403', '#F781BF', '#BC9DCC', '#A65628', '#54B0E4', '#222F75', '#1B9E77', '#B2DF8A',
                           '#E3BE00', '#FB9A99', '#E7298A', '#910241', '#00CDD1', '#A6CEE3', '#CE1261', '#5E4FA2', '#8CA77B', '#00441B', '#DEDC00', '#B3DE69', '#8DD3C7', '#999999')
                           if (n <= length(colorSpace)) {
                             colors <- colorSpace[1:n]
                           } else {
                             colors <- grDevices::colorRampPalette(colorSpace)(n)
                           }
  return(colors)
}

# Define network visualization function
netVisual_circle <- function(net, color.use = NULL, title.name = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, top = 1,
                             weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex = 1, vertex.label.color = "black",
                             edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, label.edge = FALSE, edge.label.color = 'black', edge.label.cex = 0.8,
                             edge.curved = 0.2, shape = 'circle', layout = in_circle(), margin = 0.2, vertex.size = NULL,
                             arrow.width = 1, arrow.size = 0.2) {
  
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  
  if (is.null(vertex.size.max)) {
    vertex.size.max <- if (length(unique(vertex.weight)) == 1) 5 else 15
  }
  
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1 - top)
  net[net < thresh] <- 0
  
  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  
  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = TRUE)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  coords_scale <- if (nrow(coords) != 1) scale(coords) else coords
  color.use <- if (is.null(color.use)) scPalette(length(igraph::V(g))) else color.use
  vertex.weight.max <- if (is.null(vertex.weight.max)) max(vertex.weight) else vertex.weight.max
  vertex.weight <- vertex.weight / vertex.weight.max * vertex.size.max + 5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 2] / coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 2] / coords_scale[igraph::V(g), 1]))
  
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- round(igraph::E(g)$weight, digits = 1)
  }
  edge.weight.max <- if (is.null(edge.weight.max)) max(igraph::E(g)$weight) else edge.weight.max
  igraph::E(g)$width <- if (weight.scale) 0.3 + igraph::E(g)$weight / edge.weight.max * edge.width.max else 0.3 + edge.width.max * igraph::E(g)$weight
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::V(g)$color[edge.start[, 1]], alpha.edge)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
  
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 1])] <- loop.angle[edge.start[which(edge.start[, 2] == edge.start[, 1]), 1]]
  }
  
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)), direction = -1, start = 0)
  label.dist <- vertex.weight / max(vertex.weight) + 2
  
  plot(g, edge.curved = edge.curved, vertex.shape = shape, layout = coords_scale, margin = margin, vertex.label.dist = label.dist,
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}