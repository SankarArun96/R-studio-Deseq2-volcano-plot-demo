# Installing DESEQ-2
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")

# Install ggplot2 and plotly if you do not have them
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
if (!requireNamespace("htmlwidgets", quietly = TRUE)) install.packages("htmlwidgets")

#For gene set enrichment analysis and heatmap generation
# Install Bioconductor packages for GSEA and Heatmaps
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "pheatmap"))
if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")

#For Pathway analysis using the enriched data
if (!requireNamespace("pathview", quietly = TRUE)) BiocManager::install("pathview")
BiocManager::install(c("pathview", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"))
install.packages("magick")