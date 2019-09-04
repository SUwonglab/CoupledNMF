library(uwot)
library(RColorBrewer)
library(Rtsne)
defaultPal <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"))

getColors <- function(n, palette){
	if(!is.null(palette) && n <= length(palette)){
		return(c(palette[1:n], "#B3B3B3"))
	}
	if(!is.null(palette) && n > length(palette)){
		warning("Number of clusters exceeds number of colors in palette. Default palettes used.")
	}

	if(n <= length(defaultPal)){
		return(c(defaultPal[1:n], "#B3B3B3"))
	}
	else{
		return(c(rainbow(n), "B3B3B3"))
	}
}

thresh <- function(x) {
 return(ifelse(length(which(x >= .5)) > 0, which(x >= .5), NA))
}
vis <- function(H1 = "H1.txt", H2 = "H2.txt", filename = "coupledNMF", palette = NULL, save_vars = FALSE ) {

H1 <- read.table(H1, sep = "\t")
H2 <- read.table(H2, sep = "\t")
H1H2 <- cbind(H1, H2)
h1h2colorvec <- rep(c(1, 2), c(ncol(H1), ncol(H2)))

clustlabfactor <- as.factor(apply(H1H2, 2, which.max))
thresholdclustlabs <- as.factor(apply(scale(H1H2, center=FALSE, scale=colSums(H1H2)), 2, thresh))
pal <- getColors(nrow(H1H2), palette)


na_pointvec <- replace(rep(16, length(thresholdclustlabs)), which(is.na(thresholdclustlabs)), 1)

H1H2 <- t(H1H2)


# H1H2 is samples x clusters
pca <- prcomp(H1H2, rank = 2)

# umap unsupervised
umap_unsup <- umap(H1H2)

# umap semisupervised
umap_semisup <- umap(H1H2, y = thresholdclustlabs)

# tsne
tsne <- Rtsne(H1H2)

pdf(paste(filename, ".pdf", sep = ""))

par(xpd = TRUE, mar = par()$mar + c(0,0,0,5))

# pca clust:
plot(pca$x, col = pal[clustlabfactor], pch = 16, axes = FALSE, frame.plot = TRUE, main  = "PCA", xlab = "", ylab = "") 
title(xlab = paste("PC1: ", round(summary(pca)[[1]][1], 2), "% of variance"), ylab = paste("PC2: ", round(summary(pca)[[1]][2], 2), "% of variance"), line = 0)
legend("topright", inset=c(-0.2,0), levels(clustlabfactor), fill = pal, title = "Cluster")

# pca h1h2: 
plot(pca$x, col = brewer.pal(8, "Set2")[h1h2colorvec], pch = 16, axes = FALSE, frame.plot = TRUE, main  = "PCA", xlab = "", ylab = "")
title(xlab = paste("PC1: ", round(summary(pca)[[1]][1], 2), "% of variance"), ylab = paste("PC2: ", round(summary(pca)[[1]][2], 2), "% of variance"), line = 0)
legend("topright", inset=c(-0.2,0), c("ATAC", "RNA"), fill = brewer.pal(8, "Set2"), title = "Data type")

# uu clust:
plot(umap_unsup, col = pal[clustlabfactor], pch = 16, axes = FALSE, ann = FALSE, frame.plot = TRUE) 
title(main  = "Unsupervised UMAP")
legend("topright", inset=c(-0.2,0), levels(clustlabfactor), fill = pal, title = "Cluster")

# uu h1h2:
plot(umap_unsup, col = brewer.pal(8, "Set2")[h1h2colorvec], pch = 16, axes = FALSE, ann = FALSE, frame.plot = TRUE)
title(main  = "Unsupervised UMAP")
legend("topright", inset=c(-0.2,0), c("ATAC", "RNA"), fill = brewer.pal(8, "Set2"), title = "Data type")

# us clust
plot(umap_semisup, col = pal[clustlabfactor], pch = na_pointvec, axes = FALSE, ann = FALSE, frame.plot = TRUE) 
title(main = "Semisupervised UMAP (unclustered values unfilled)")
legend("topright", inset=c(-0.2,0), levels(clustlabfactor), fill = pal, title = "Cluster")

# us h1h2
plot(umap_semisup, col = brewer.pal(8, "Set2")[h1h2colorvec], pch = na_pointvec, axes = FALSE, ann = FALSE, frame.plot = TRUE)
title(main = "Semisupervised UMAP (unclustered values unfilled)")
legend("topright", inset=c(-0.2,0), c("ATAC", "RNA"), fill = brewer.pal(8, "Set2"), title = "Data type")

# tsne clust
plot(tsne$Y, col = pal[clustlabfactor], pch = 16, axes = FALSE, ann = FALSE, frame.plot = TRUE) 
title(main = "t-SNE")
legend("topright", inset=c(-0.2,0), levels(clustlabfactor), fill = pal, title = "Cluster")

# tsne h1h2
plot(tsne$Y, col = brewer.pal(8, "Set2")[h1h2colorvec], pch = 16, axes = FALSE, ann = FALSE, frame.plot = TRUE)
title(main = "t-SNE")
legend("topright", inset=c(-0.2,0), c("ATAC", "RNA"), fill = brewer.pal(8, "Set2"), title = "Data type")

dev.off()


if(save_vars == TRUE){
	list2env(mget(ls()), envir = .GlobalEnv)
}

}


