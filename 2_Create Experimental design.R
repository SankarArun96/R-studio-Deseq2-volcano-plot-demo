# Create a factor defining the two groups in the exact order of your count matrix columns
condition <- factor(c(
  "Adipose", "Adipose", "Adipose",         # SRR86+87, SRR96+97, SRR06+07
  "Bone_Marrow", "Bone_Marrow", "Bone_Marrow" # SRR27603585, SRR27603586, SRR27603587
))

# Create the colData dataframe
col_data <- data.frame(
  row.names = colnames(count_matrix),
  condition = condition
)