# vis() for visualization of coupledNMF results

## Introduction
visualization.R presents a function vis() which can be used to generate relevant graphs for viewing coupledNMF results using different dimensionality reduction techniques. The following graphs are generated into one pdf file:
- PCA, colored by cluster label
- PCA, colored by data type (ATAC vs. RNAseq)
- unsupervised UMAP, colored by cluster label
- unsupervised UMAP, colored by data type (ATAC vs. RNAseq)
- semisupervised UMAP, colored by cluster label, with unclustered cells as empty circles
- semisupervised UMAP, colored by data type, with unclustered cells as empty circles
- t-SNE, colored by cluster label
- t-SNE, colored by data type (ATAC vs. RNAseq)

## Running vis()
> vis() function call:
>`vis(H1 = "H1.txt", H2 = "H2.txt", filename = "coupledNMF", palette = NULL, save_vars = FALSE)`

- `H1 = "H1.txt"` is the filename of the H1 matrix with dimensions [number of clusters] x [number of ATACseq samples]; saved as a text file with tab separation (sep = "\t")
- `H2 = "H2.txt"` is the filename of the H2 matrix with dimensions [number of clusters] x [number of RNAseq samples]; saved as a text file with tab separation (sep = "\t")
- `filename = "coupledNMF"` changes the name of the pdf file
- `palette = NULL` can be changed so that graphs are plotted using a specific color palette
- `save_vars = FALSE` can be set to TRUE to save all variables used in plotting (cluster labels, appended H1H2 matrix, PCA/UMAP/t-SNE results, etc. to the global environment)

#### Output:
vis() saves a pdf file called [filename].pdf (defaults to coupledNMF.pdf) to the local directory, containing all graphs.
 
## Requirements
- uwot
- RColorBrewer
- Rtsne

**For questions about vis(), please contact mirandal@stanford.edu.**