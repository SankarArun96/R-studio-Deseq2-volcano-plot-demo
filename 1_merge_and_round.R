# Load the DESeq2 package
library(DESeq2)

# Load the Adipose and Bone Marrow count files
# Ensure the files are in your current working directory
adipose_df <- read.delim("ADP.tabular", row.names = 1, check.names = FALSE)
bm_df <- read.delim("BM.tabular", row.names = 1, check.names = FALSE)

# Merge the matrices by gene names
merged_df <- merge(adipose_df, bm_df, by = "row.names", all = FALSE)

# Set the row names to the gene names and remove the temporary Row.names column
rownames(merged_df) <- merged_df$Row.names
merged_df$Row.names <- NULL

# Convert the dataframe to a matrix and round the estimated counts to integers
count_matrix <- round(as.matrix(merged_df))

# Export the rounded count matrix to a CSV file
write.csv(count_matrix, file = "merged_rounded_counts.csv", row.names = TRUE)