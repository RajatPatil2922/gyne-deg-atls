
library(DESeq2)
library(GEOquery)

setwd("D:/Rajat_atls/gyne-deg-atls")

print("--> Step 1: Parsing metadata and assembling count matrix...")
gse216609 <- getGEO("GSE216609", GSEMatrix=TRUE, getGPL=FALSE)[[1]]
meta216609 <- pData(gse216609)

# Build a mapping data frame between GSM IDs, file names, and disease groups
meta_map <- data.frame(
  gsm = rownames(meta216609),
  title = meta216609$title,
  group = ifelse(grepl("PCOS", meta216609$title), "case", "control"),
  stringsAsFactors = FALSE
)

# Find all extracted sample files
extracted_dir <- "data/expr/GSE216609/extracted"
sample_files <- list.files(extracted_dir, pattern = "\\.tsv\.gz$", full.names = TRUE)

# Read the first file to initialize our master matrix structure
first_sample <- read.delim(sample_files[1], header=TRUE, stringsAsFactors=FALSE)
# Assume column 1 is Gene and column 2 is the raw count values
gene_names <- first_sample[, 1]
counts_list <- list()

# Loop through all files to extract raw counts cleanly
for (f in sample_files) {
  fname <- basename(f)
  gsm_id <- strsplit(fname, "_")[[1]][1]
  
  # Match the file back to its metadata entry
  match_idx <- match(gsm_id, meta_map$gsm)
  if(!is.na(match_idx)) {
    sdata <- read.delim(f, header=TRUE, stringsAsFactors=FALSE)
    counts_list[[gsm_id]] <- sdata[, 2]
  }
}

# Bind into a solid matrix layout
counts_matrix <- do.call(cbind, counts_list)
rownames(counts_matrix) <- gene_names

# Align the metadata groups perfectly with the column matrix order
meta_final <- meta_map[match(colnames(counts_matrix), meta_map$gsm), ]
grp_factor <- factor(meta_final$group, levels=c("control", "case"))

print("--- Final Verified Sample Breakdown Summary ---")
print(table(grp_factor))

print("--> Step 2: Running statistical estimation via DESeq2...")
dds <- DESeqDataSetFromMatrix(round(as.matrix(counts_matrix)), data.frame(group=grp_factor), ~group)
dds <- dds[rowSums(counts(dds)) >= 10, ] 
dds_executed <- DESeq(dds)
res_data <- as.data.frame(results(dds_executed, contrast=c("group", "case", "control")))

# Step 3: Construct the precise 5-column output format requested
final_table <- data.frame(
  gene   = rownames(res_data),  
  log2FC = res_data$log2FoldChange,
  p      = res_data$pvalue,
  padj   = res_data$padj,
  n      = ncol(counts_matrix)
)

# Sort strictly by significance
final_table <- final_table[order(final_table$p), ]

# Export to your standardized results folder
dir.create("results/deg", recursive=TRUE, showWarnings=FALSE)
write.table(final_table, "results/deg/GSE216609.RNAseq.deg.tsv", sep="\\t", quote=FALSE, row.names=FALSE)

print("--> SUCCESS: GSE216609.RNAseq.deg.tsv completed and formatted perfectly!")
