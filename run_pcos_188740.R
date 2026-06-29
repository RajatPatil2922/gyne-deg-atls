
library(DESeq2)

setwd("D:/Rajat_atls/gyne-deg-atls")

# Load the verified count matrix path
counts_file <- "data/expr/GSE188740/GSE188740_counts.tsv.gz"
counts188740 <- read.delim(counts_file, row.names=1, check.names=FALSE)

# Isolate only numerical data columns
counts188740 <- counts188740[, sapply(counts188740, is.numeric)]

# Classify groups using the first letter prefix rule
sample_names <- colnames(counts188740)
grp_labels <- ifelse(substr(sample_names, 1, 1) == "P", "case", "control")
grp188740 <- factor(grp_labels, levels=c("control","case"))

print("--- GSE188740 Sample Breakdown Summary ---")
print(table(grp188740))

# Execute the DESeq2 pipeline
dds188740 <- DESeqDataSetFromMatrix(round(as.matrix(counts188740)), data.frame(group=grp188740), ~group)
dds188740 <- dds188740[rowSums(counts(dds188740)) >= 10, ] 
dds_executed <- DESeq(dds188740)
res_data <- as.data.frame(results(dds_executed, contrast=c("group","case","control")))

# Format the 5-column layout template
final_table <- data.frame(
  gene   = sub("\\..*", "", rownames(res_data)),  
  log2FC = res_data$log2FoldChange,
  p      = res_data$pvalue,
  padj   = res_data$padj,
  n      = ncol(counts188740)
)

# Sort strictly by raw p-value
final_table <- final_table[order(final_table$p), ]

# Export the standardized matrix file to your results directory
dir.create("results/deg", recursive=TRUE, showWarnings=FALSE)
write.table(final_table, "results/deg/GSE188740.RNAseq.deg.tsv", sep="\\t", quote=FALSE, row.names=FALSE)

print("--> SUCCESS: GSE188740.RNAseq.deg.tsv completed and saved successfully!")
