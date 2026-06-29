
library(DESeq2)

setwd("D:/Rajat_atls/gyne-deg-atls")

downloaded_files <- list.files("data/expr", pattern = "counts.tsv.gz", recursive = TRUE, full.names = TRUE)
counts199225 <- read.delim(downloaded_files[1], row.names=1, check.names=FALSE)
counts199225 <- counts199225[, sapply(counts199225, is.numeric)]

sample_names <- colnames(counts199225)
grp_labels <- ifelse(substr(sample_names, 1, 1) == "P", "case", "control")
grp199225 <- factor(grp_labels, levels=c("control","case"))

dds199225 <- DESeqDataSetFromMatrix(round(as.matrix(counts199225)), data.frame(group=grp199225), ~group)
dds199225 <- dds199225[rowSums(counts(dds199225)) >= 10, ] 
dds_executed <- DESeq(dds199225)
res_data <- as.data.frame(results(dds_executed, contrast=c("group","case","control")))

final_table <- data.frame(
  gene   = sub("\\..*", "", rownames(res_data)),  
  log2FC = res_data$log2FoldChange,
  p      = res_data$pvalue,
  padj   = res_data$padj,
  n      = ncol(counts199225)
)

final_table <- final_table[order(final_table$p), ]

dir.create("results/deg", recursive=TRUE, showWarnings=FALSE)

# --- Updated filename to include RNAseq ---
write.table(final_table, "results/deg/GSE199225.RNAseq.deg.tsv", sep="\\t", quote=FALSE, row.names=FALSE)

print("--> SUCCESS: GSE199225.RNAseq.deg.tsv written to disk!")
