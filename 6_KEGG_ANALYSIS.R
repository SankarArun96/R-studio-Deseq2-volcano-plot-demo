# ============================================================
# Generate TIFF maps for ALL statistically significant KEGG
# pathways using DESeq2 log2 fold changes
# Comparison: Bone Marrow vs Adipose
#
# Colour scheme: DOWNREGULATED = RED, UPREGULATED = BLUE
# Every image carries (a) pathview's colour-scale bar and
# (b) a caption explaining the colours and the scale range.
# ============================================================

# ----------------------------
# 0. Load packages
# ----------------------------

# List every package the script depends on
required_packages <- c(
  "pathview",         # draws KEGG pathway maps coloured by gene values
  "clusterProfiler",  # provides bitr() for gene ID conversion
  "org.Hs.eg.db",     # human gene annotation database (SYMBOL <-> ENTREZ)
  "AnnotationDbi",    # backend used by org.Hs.eg.db
  "magick"            # image processing (PNG -> TIFF, adding captions)
)

# Loop over each required package name
for (pkg in required_packages) {
  # Check whether the package is installed (without loading it)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Stop with a clear message if a package is missing
    stop(
      paste0(
        "Package '", pkg,
        "' is missing. Install it before running this script."
      )
    )
  }
}

library(pathview)         # load pathview into the session
library(clusterProfiler)  # load clusterProfiler into the session
library(org.Hs.eg.db)     # load the human annotation database
library(AnnotationDbi)    # load the annotation interface
library(magick)           # load the image library

# ----------------------------
# USER SETTINGS
# ----------------------------

de_file   <- "DESeq2_BoneMarrow_vs_Adipose_Results.csv"  # DESeq2 results table (gene symbols as row names)
gsea_file <- "GSEA_KEGG_Pathways_Results.csv"            # GSEA/KEGG enrichment results table

output_dir <- "KEGG_Pathway_TIFFs"  # folder where all TIFFs and the summary are saved

padj_cutoff  <- 0.05   # a pathway is "significant" if its adjusted p-value is below this
species_code <- "hsa"  # KEGG organism code (hsa = Homo sapiens)

# Colour-scale range (log2 fold change). Colours saturate at +/- this value.
# NA = choose automatically (99th percentile of |log2FC|, so a few extreme
# genes don't wash out the colours of everything else).
# Set a number such as 2 to force the same scale, e.g. -2 to +2.
fc_limit_user <- NA

# Wording used in the caption written onto every image
caption_comparison <- "Bone Marrow vs Adipose"

# ----------------------------
# 1. Read DESeq2 results
# ----------------------------

# Make the paths absolute NOW, because the script later changes working directory
de_file   <- normalizePath(de_file,   mustWork = FALSE)  # absolute path to DE file
gsea_file <- normalizePath(gsea_file, mustWork = FALSE)  # absolute path to GSEA file

# Stop early if the DESeq2 file cannot be found
if (!file.exists(de_file)) {
  stop("DESeq2 results file not found: ", de_file)
}

# Read the DESeq2 table; first column (gene symbols) becomes row names
res_df <- read.csv(
  de_file,
  row.names = 1,        # use the first column as row names
  check.names = FALSE   # keep column names exactly as written in the file
)

# The fold-change column is essential; stop if it is absent
if (!"log2FoldChange" %in% names(res_df)) {
  stop("DESeq2 results must contain a log2FoldChange column.")
}

# Copy the row names (gene symbols) into an explicit SYMBOL column
res_df$SYMBOL <- rownames(res_df)

# Keep only rows that have a real, finite log2 fold change
res_df <- res_df[
  !is.na(res_df$log2FoldChange) &       # drop NA values
    is.finite(res_df$log2FoldChange),   # drop Inf / -Inf values
  ,
  drop = FALSE                          # keep data-frame structure
]

# Stop if nothing is left after filtering
if (nrow(res_df) == 0) {
  stop("No valid log2FoldChange values were found.")
}

# ----------------------------
# 2. Convert SYMBOL to ENTREZID
# ----------------------------

# KEGG/pathview needs Entrez IDs, so convert the gene symbols
gene_df <- clusterProfiler::bitr(
  unique(res_df$SYMBOL),   # each symbol once
  fromType = "SYMBOL",     # input ID type
  toType   = "ENTREZID",   # output ID type
  OrgDb    = org.Hs.eg.db, # human annotation database
  drop     = TRUE          # drop symbols that could not be mapped
)

# Stop if no symbols mapped (usually means wrong species / ID type)
if (is.null(gene_df) || nrow(gene_df) == 0) {
  stop(
    "No gene symbols mapped to Entrez IDs. ",
    "Check whether the input row names are human gene symbols."
  )
}

# Attach each gene's log2 fold change to its Entrez ID
mapped_df <- merge(
  gene_df,                                   # SYMBOL + ENTREZID table
  res_df[, c("SYMBOL", "log2FoldChange")],   # SYMBOL + fold change table
  by = "SYMBOL"                              # join on the symbol
)

# Remove any rows with missing Entrez ID or non-finite fold change
mapped_df <- mapped_df[
  !is.na(mapped_df$ENTREZID) &
    is.finite(mapped_df$log2FoldChange),
  ,
  drop = FALSE
]

# Sort so the largest |log2FC| comes first (used to resolve duplicates)
mapped_df <- mapped_df[
  order(abs(mapped_df$log2FoldChange), decreasing = TRUE),
  ,
  drop = FALSE
]

# If several symbols share one Entrez ID, keep only the strongest one
mapped_df <- mapped_df[
  !duplicated(mapped_df$ENTREZID),
  ,
  drop = FALSE
]

# Build the named numeric vector pathview expects:
# values = log2FC, names = Entrez IDs
fold_changes <- mapped_df$log2FoldChange
names(fold_changes) <- as.character(mapped_df$ENTREZID)

# Stop if no usable pairs remain
if (length(fold_changes) == 0) {
  stop("No usable Entrez ID / fold-change pairs remain.")
}

# Report how many genes made it through
message("Genes mapped to Entrez IDs: ", length(fold_changes))

# ----------------------------
# 2b. Decide the colour-scale limit (same for every pathway)
# ----------------------------

if (!is.na(fc_limit_user)) {
  # User supplied a limit: use it
  fc_limit <- fc_limit_user
} else {
  # Otherwise use the 99th percentile of absolute fold changes
  fc_limit <- as.numeric(quantile(abs(fold_changes), 0.99, na.rm = TRUE))
}

# Round up to 1 decimal so the scale labels are tidy
fc_limit <- ceiling(fc_limit * 10) / 10

# Safety net: if the limit is invalid or zero, fall back to 1
if (!is.finite(fc_limit) || fc_limit <= 0) {
  fc_limit <- 1
}

# Report the scale that will be used
message("Colour scale: -", fc_limit, " to +", fc_limit, " (log2 fold change)")

# ----------------------------
# 3. Read and filter KEGG results
# ----------------------------

# Stop early if the GSEA/KEGG file cannot be found
if (!file.exists(gsea_file)) {
  stop("KEGG results file not found: ", gsea_file)
}

# Read the enrichment results table
kegg_results <- read.csv(
  gsea_file,
  check.names = FALSE,       # keep original column names
  stringsAsFactors = FALSE   # keep text as plain strings
)

# Pathway IDs are required to know which maps to draw
if (!"ID" %in% names(kegg_results)) {
  stop("KEGG results must contain an ID column.")
}

# Possible names for the adjusted p-value column
padj_candidates <- c(
  "p.adjust", "padj", "FDR",
  "qvalue", "qvalues", "p.adjusted"
)

# Pick the first candidate that exists in the table
padj_col <- padj_candidates[
  padj_candidates %in% names(kegg_results)
][1]

# Stop if none of the candidates is present
if (is.na(padj_col)) {
  stop(
    "No adjusted p-value column found. Expected one of: ",
    paste(padj_candidates, collapse = ", ")
  )
}

# Possible names for the pathway description column
description_candidates <- c(
  "Description", "description", "Pathway",
  "pathway", "Name", "name"
)

# Pick the first description column that exists
description_col <- description_candidates[
  description_candidates %in% names(kegg_results)
][1]

# If there is none, mark it as NA (pathway ID will be used as the name)
if (is.na(description_col)) {
  description_col <- NA_character_
}

# Convert the adjusted p-values to numbers (non-numeric -> NA, silently)
kegg_results$adjusted_p <- suppressWarnings(
  as.numeric(as.character(kegg_results[[padj_col]]))
)

# Keep only significant pathways with a valid ID
significant <- kegg_results[
  !is.na(kegg_results$adjusted_p) &            # p-value present
    is.finite(kegg_results$adjusted_p) &       # p-value is a real number
    kegg_results$adjusted_p < padj_cutoff &    # passes the significance cutoff
    !is.na(kegg_results$ID) &                  # ID present
    nzchar(as.character(kegg_results$ID)),     # ID not empty
  ,
  drop = FALSE
]

# Sort by adjusted p-value (most significant first)
significant <- significant[
  order(significant$adjusted_p),
  ,
  drop = FALSE
]

# Remove duplicated pathway IDs (keeps the most significant copy)
significant <- significant[
  !duplicated(significant$ID),
  ,
  drop = FALSE
]

# Stop if no pathway passed the cutoff
if (nrow(significant) == 0) {
  stop(
    "No significant pathways found at adjusted p-value < ",
    padj_cutoff,
    ". Check the input table and significance threshold."
  )
}

# Report how many pathways will be drawn
message(
  "Significant pathways to process: ",
  nrow(significant)
)

# ----------------------------
# 4. Prepare output folder
# ----------------------------

# Create the output folder (and parents) if it does not exist
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Convert to an absolute path so it still works after setwd()
output_dir <- normalizePath(output_dir)

# Helper: turn any text into a safe filename
safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)  # replace unsafe characters with "_"
  x <- gsub("_+", "_", x)               # collapse repeated underscores
  x <- gsub("^_|_$", "", x)             # trim underscores at start/end
  
  # If nothing useful is left, use a placeholder name
  if (is.na(x) || !nzchar(x)) {
    x <- "Unnamed_Pathway"
  }
  
  substr(x, 1, 120)  # limit the length to 120 characters
}

# ----------------------------
# 5. Render every significant pathway
# ----------------------------

# Pre-allocate a list that will hold one result row per pathway
summary_rows <- vector("list", nrow(significant))

# IMPORTANT: pathview saves its output PNG in the CURRENT working directory,
# so we move into the output folder while drawing and restore it afterwards.
original_wd <- getwd()                 # remember where we started
setwd(output_dir)                      # work inside the output folder
on.exit(setwd(original_wd), add = TRUE)  # restore automatically if the script errors

# Loop over every significant pathway
for (i in seq_len(nrow(significant))) {
  
  # Get this pathway's ID as text
  pathway_id <- as.character(significant$ID[i])
  
  # Strip a leading "hsa" so both "hsa00010" and "00010" are accepted
  pathway_id <- sub("^hsa", "", pathway_id)
  
  # A KEGG ID must be exactly 5 digits; otherwise skip it
  if (!grepl("^[0-9]{5}$", pathway_id)) {
    warning(
      "Skipping unsupported KEGG pathway ID: ",
      significant$ID[i]
    )
    next  # move on to the next pathway
  }
  
  # Keep the 5-digit number for pathview, and the full ID for naming
  pathway_num <- pathway_id                       # e.g. "00010"
  pathway_id  <- paste0("hsa", pathway_num)       # e.g. "hsa00010"
  
  # Use the human-readable description if the table has one
  if (!is.na(description_col)) {
    pathway_name <- as.character(
      significant[[description_col]][i]
    )
  } else {
    pathway_name <- pathway_id  # otherwise fall back to the ID
  }
  
  # If the description is missing/empty, use the ID instead
  if (is.na(pathway_name) || !nzchar(pathway_name)) {
    pathway_name <- pathway_id
  }
  
  # Build the base filename: pathway name + ID (prevents collisions)
  file_stem <- paste0(
    safe_filename(pathway_name),
    "_",
    pathway_id
  )
  
  # Print progress to the console
  message(
    "[", i, "/", nrow(significant), "] ",
    pathway_name, " (", pathway_id, ")"
  )
  
  # Suffix pathview adds to its output PNG name; used to find the file again
  suffix <- paste0(
    "BoneMarrow_vs_Adipose_",
    file_stem
  )
  
  # tryCatch so one failing pathway doesn't stop the whole run
  result <- tryCatch({
    
    # ---- Draw the KEGG map with pathview ----
    pathview::pathview(
      gene.data   = fold_changes,           # named vector of log2FC (names = Entrez IDs)
      pathway.id  = pathway_num,            # 5-digit pathway number
      species     = species_code,           # "hsa"
      gene.idtype = "entrez",               # our IDs are Entrez IDs
      kegg.native = TRUE,                   # draw on the native KEGG image (PNG)
      kegg.dir    = ".",                    # download/cache KEGG templates here
      limit       = list(gene = fc_limit),  # scale runs from -fc_limit to +fc_limit
      low         = list(gene = "red"),     # most NEGATIVE values (downregulated) = RED
      mid         = list(gene = "white"),   # values near zero = WHITE
      high        = list(gene = "blue"),    # most POSITIVE values (upregulated) = BLUE
      na.col      = "grey",                 # genes with no data = grey
      bins        = list(gene = 20),        # 20 colour steps (smooth gradient)
      plot.col.key = TRUE,                  # draw the colour-scale bar on the image
      key.pos     = "bottomright",          # put the scale bar at the bottom right
      out.suffix  = suffix                  # unique suffix added to the output file name
    )
    
    # ---- Locate the PNG that pathview just created ----
    # Pathview names it "<pathway_id>.<suffix>.png" in the working directory
    png_file <- paste0(pathway_id, ".", suffix, ".png")
    
    # Confirm it exists; otherwise raise an error (caught below)
    if (!file.exists(png_file)) {
      stop("Pathview PNG not found for this pathway.")
    }
    
    # ---- Read the PNG into magick ----
    img <- magick::image_read(png_file)
    
    # ---- Add a caption explaining colours and scale ----
    caption_text <- paste0(
      "log2 fold change (", caption_comparison, "):  ",
      "RED = downregulated   |   WHITE = no change   |   BLUE = upregulated   |   ",
      "scale: -", fc_limit, " to +", fc_limit
    )
    
    # Get the image width so the caption bar spans the full image
    img_width <- magick::image_info(img)$width
    
    # Create a white strip (same width, 60 px tall) to hold the caption
    caption_bar <- magick::image_blank(
      width  = img_width,
      height = 60,
      color  = "white"
    )
    
    # Write the caption text into the white strip
    caption_bar <- magick::image_annotate(
      caption_bar,
      text     = caption_text,
      gravity  = "center",   # centre the text in the strip
      size     = max(14, round(img_width / 70)),  # text size scales with image width
      color    = "black"
    )
    
    # Stack the pathway image on top and the caption strip below it
    img <- magick::image_append(c(img, caption_bar), stack = TRUE)
    
    # ---- Write the final TIFF ----
    tiff_file <- file.path(
      output_dir,
      paste0(file_stem, ".tiff")
    )
    
    magick::image_write(
      img,
      path    = tiff_file,    # destination file
      format  = "tiff",       # output format
      density = "300x300"     # store 300 dpi in the file metadata
    )
    
    # Confirm the TIFF really exists on disk
    if (!file.exists(tiff_file)) {
      stop("TIFF output was not created.")
    }
    
    # Return a one-row "success" record for the summary table
    data.frame(
      pathway_id       = pathway_id,
      pathway_name     = pathway_name,
      adjusted_p_value = significant$adjusted_p[i],
      status           = "Success",
      tiff_file        = tiff_file,
      message          = "",
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    
    # On error: print a warning with the reason
    warning(
      "Failed for ", pathway_name, ": ",
      conditionMessage(e)
    )
    
    # Return a one-row "failed" record for the summary table
    data.frame(
      pathway_id       = pathway_id,
      pathway_name     = pathway_name,
      adjusted_p_value = significant$adjusted_p[i],
      status           = "Failed",
      tiff_file        = NA_character_,
      message          = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  
  # Store this pathway's result row in the list
  summary_rows[[i]] <- result
}

# Go back to the original working directory now that drawing is finished
setwd(original_wd)

# ----------------------------
# 6. Save processing summary
# ----------------------------

# Combine all result rows (skipping pathways that were skipped = NULL)
summary_df <- do.call(
  rbind,
  summary_rows[!vapply(summary_rows, is.null, logical(1))]
)

# Path of the summary CSV inside the output folder
summary_file <- file.path(
  output_dir,
  "KEGG_Pathway_Rendering_Summary.csv"
)

# Save the summary table
write.csv(
  summary_df,
  summary_file,
  row.names = FALSE
)

# Final report in the console
message("========================================")
message("Pathway rendering finished.")
message("Successful TIFFs: ", sum(summary_df$status == "Success"))
message("Failed pathways: ", sum(summary_df$status == "Failed"))
message("Output directory: ", output_dir)
message("Summary: ", summary_file)
message("========================================")