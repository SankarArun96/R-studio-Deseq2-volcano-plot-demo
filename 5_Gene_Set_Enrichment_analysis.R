library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(DESeq2)
library(plotly)

###########################################################################################
#------------------------------------------------------------------------------------------
#Gene Set Enrichment analysis
#------------------------------------------------------------------------------------------

# Read the DESeq2 results CSV
res_df <- read.csv("DESeq2_BoneMarrow_vs_Adipose_Results.csv", row.names = 1)

# Drop rows with missing log2FoldChange values
res_df <- res_df[!is.na(res_df$log2FoldChange), ]

# Create a ranked vector of genes based on Log2 Fold Change
gene_list <- res_df$log2FoldChange
names(gene_list) <- rownames(res_df)
gene_list <- sort(gene_list, decreasing = TRUE)

# Run GSEA targeting Biological Processes (BP)
gsea_res <- gseGO(
  geneList = gene_list,
  ont = "BP",
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  pvalueCutoff = 0.05,
  eps = 0 # Prevents p-value underflow errors in highly significant datasets
)

# View the top enriched pathways
head(as.data.frame(gsea_res)[, c("ID", "Description", "enrichmentScore", "p.adjust")])

# ---------------------------------------------------------
# EXPORT AND INTERACTIVE VISUALIZATION
# ---------------------------------------------------------

# 1. Export the full GSEA results to a CSV file
write.csv(as.data.frame(gsea_res), file = "GSEA_Biological_Processes_Results.csv", row.names = FALSE)

# 2a. Generate a static ggplot2 dotplot of the top 20 enriched pathways
# Note: dotplot requires the 'enrichplot' package loaded above
static_dotplot <- dotplot(gsea_res, showCategory = 20, title = "Top 20 Enriched Biological Processes")

# 2b.Save as a TIFF (often required for publication)
ggsave("GSEA_Dotplot.tiff", plot = static_dotplot, device = "tiff", width = 10, height = 8, dpi = 300, compression = "lzw")

# 3. Convert the ggplot object to an interactive WebGL-accelerated Plotly object
interactive_gsea <- ggplotly(static_dotplot) %>%
  toWebGL()


# 4. Save the interactive plot as a standalone HTML file
htmlwidgets::saveWidget(interactive_gsea, "Interactive_GSEA_Dotplot.html")

###########################################################################################