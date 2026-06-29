
library(DESeq2)

setwd("D:/Rajat_atls/gyne-deg-atls")

dir.create("results/deg", recursive=TRUE, showWarnings=FALSE)

# ==========================================================
# 1. PROCESS THE GV (OOCYTE) DATASET BRANCH
# ==========================================================
print("--> Starting RNA-seq processing for GSE155489 (GV Oocytes)...")
gv_file <- "data/expr/GSE155489/GSE155489_gv_pcos_counts.csv.gz"
gv_counts <- read.csv(gv_file, row.names=1, check.names=FALSE)
gv_counts <- gv_counts[, sapply(gv_counts, is.numeric)]

# 6 controls followed by 6 cases based on expression profile transitions
gv_groups <- factor(c(rep("control", 6), rep("case", 6)), levels=c("control","case"))

dds_gv <- DESeqDataSetFromMatrix(round(as.matrix(gv_counts)), data.frame(group=gv_groups), ~group)
dds_gv <- dds_gv[rowSums(counts(dds_gv)) >= 10, ] 
res_gv <- as.data.frame(results(DESeq(dds_gv), contrast=c("group","case","control")))

final_gv <- data.frame(
  gene   = rownames(res_gv),  
  log2FC = res_gv$log2FoldChange,
  p      = res_gv$pvalue,
  padj   = res_gv$padj,
  n      = ncol(gv_counts)
)
final_gv <- final_gv[order(final_gv$p), ]
write.table(final_gv, "results/deg/GSE155489_GV.RNAseq.deg.tsv", sep="\\t", quote=FALSE, row.names=FALSE)
print("Saved: results/deg/GSE155489_GV.RNAseq.deg.tsv")

# ==========================================================
# 2. PROCESS THE GC (GRANULOSA CELL) DATASET BRANCH
# ==========================================================
print("--> Starting RNA-seq processing for GSE155489 (GC Granulosa)...")
gc_file <- "data/expr/GSE155489/GSE155489_gc_pcos_counts.csv.gz"
gc_counts <- read.csv(gc_file, row.names=1, check.names=FALSE)
gc_counts <- gc_counts[, sapply(gc_counts, is.numeric)]

# 4 controls followed by 4 cases based on metadata sequencing targets
gc_groups <- factor(c(rep("control", 4), rep("case", 4)), levels=c("control","case"))

dds_gc <- DESeqDataSetFromMatrix(round(as.matrix(gc_counts)), data.frame(group=gc_groups), ~group)
dds_gc <- dds_gc[rowSums(counts(dds_gc)) >= 10, ] 
res_gc <- as.data.frame(results(DESeq(dds_gc), contrast=c("group","case","control")))

final_gc <- data.frame(
  gene   = rownames(res_gc),  
  log2FC = res_gc$log2FoldChange,
  p      = res_gc$pvalue,
  padj   = res_gc$padj,
  n      = ncol(gc_counts)
)
final_gc <- final_gc[order(final_gc$p), ]
write.table(final_gc, "results/deg/GSE155489_GC.RNAseq.deg.tsv", sep="\\t", quote=FALSE, row.names=FALSE)
print("Saved: results/deg/GSE155489_GC.RNAseq.deg.tsv")

print("--> SUCCESS: Both GSE155489 sample sets completely processed!")
