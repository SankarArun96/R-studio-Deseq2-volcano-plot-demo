# Construct the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ condition
)

# Filter out lowly expressed genes (optional but speeds up computation)
# Keeps genes with at least 10 reads across at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run the standard DESeq2 pipeline
dds <- DESeq(dds)

# Extract the results, contrasting Bone Marrow against Adipose
res <- results(dds, contrast = c("condition", "Bone_Marrow", "Adipose"))

# Sort the results by adjusted p-value to see the most significant genes at the top
res_sorted <- res[order(res$padj), ]
head(res_sorted)

# Export the sorted differential expression results to a CSV file
write.csv(as.data.frame(res_sorted), file = "DESeq2_BoneMarrow_vs_Adipose_Results.csv", row.names = TRUE)

###############################################################################################################################


