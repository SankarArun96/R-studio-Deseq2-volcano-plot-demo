# Read the DESeq2 results from the generated CSV file
res_df <- read.csv("DESeq2_BoneMarrow_vs_Adipose_Results.csv", row.names = 1)

# Drop genes with NA in the adjusted p-value column
res_df <- res_df[!is.na(res_df$padj), ]

# Create a new column categorizing each gene's significance
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up in Bone Marrow"
res_df$Significance[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Up in Adipose"

# Add a 'Gene' column from the row names for the interactive tooltip
res_df$Gene <- rownames(res_df)

######################################################################################################

#Preparing Basic Volcano Plot

library(ggplot2)
library(plotly)

# Build the ggplot object
basic_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), 
                                    color = Significance, text = Gene)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up in Adipose" = "blue", 
                                "Not Significant" = "grey80", 
                                "Up in Bone Marrow" = "red")) +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(title = "Volcano Plot: Bone Marrow vs Adipose",
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)") +
  theme(legend.position = "bottom")

# View the basic plot in the RStudio 'Plots' pane
print(basic_volcano)

# Save the plot as a high-resolution TIFF file
ggsave("Volcano_Plot_BM_vs_Adipose.tiff", plot = basic_volcano, 
       device = "tiff", width = 8, height = 6, dpi = 300, compression = "lzw")

#########################################################################################################

#Best for small data points

# Convert the basic ggplot into an interactive HTML widget
interactive_volcano <- ggplotly(basic_volcano, tooltip = c("text", "x", "y"))

# View the interactive plot in the RStudio 'Viewer' pane
interactive_volcano

#View the plot in your browser
htmlwidgets::saveWidget(interactive_volcano, "Interactive_Volcano_Plot.html")

##########################################################################################################

#Best for Large data points

# Convert to an interactive plot using WebGL for high performance
interactive_volcano <- toWebGL(ggplotly(basic_volcano, tooltip = c("text", "x", "y")))

# Save the hardware-accelerated widget
htmlwidgets::saveWidget(interactive_volcano, "Interactive_Volcano_WebGL.html")