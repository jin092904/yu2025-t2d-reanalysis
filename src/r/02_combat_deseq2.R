# ComBat-Seq 로 배치를 보정하고 DESeq2 로 OT2D vs OH 차등발현을 구한다.
# 사전 명세(2026-11-09 고정): batch = GEO 시리즈, group 없음, 유전자 사전 거르기 없음,
#   DESeq2 ~condition 기본값, LFC 축소 없음, DEG = padj < 0.05 & |log2FC| >= 2.
# 실행: Rscript src/r/02_combat_deseq2.R  (프로젝트 폴더에서)
# 탐색용 인자: "group" = ComBat-Seq 에 group = condition, "pc" = NCBI 주석표의 protein-coding 유전자만 쓴다.
suppressPackageStartupMessages({library(sva); library(DESeq2)})

args <- commandArgs(trailingOnly = TRUE)
use_group <- "group" %in% args
use_pc <- "pc" %in% args
tag <- paste0("R", if (use_group) "_group" else "", if (use_pc) "_pc" else "")
annot <- read.delim(gzfile("data/raw/Human.GRCh38.p13.annot.tsv.gz"), quote = "")[, c("GeneID", "Symbol", "GeneType")]

samples <- read.delim("data/processed/samples.tsv")
counts <- as.matrix(read.delim(gzfile("data/processed/counts_main.tsv.gz"), row.names = 1, check.names = FALSE))
meta <- samples[match(colnames(counts), samples$gsm), ]
meta$condition <- factor(meta$paper_group, levels = c("OH", "OT2D"))
stopifnot(identical(meta$gsm, colnames(counts)))
if (use_pc) counts <- counts[rownames(counts) %in% annot$GeneID[annot$GeneType == "protein-coding"], ]

adj <- ComBat_seq(counts, batch = meta$series, group = if (use_group) meta$condition else NULL)

dds <- DESeqDataSetFromMatrix(adj, colData = meta, design = ~condition)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("condition", "OT2D", "OH")))

res <- cbind(GeneID = as.integer(rownames(res)), res)
res$Symbol <- annot$Symbol[match(res$GeneID, annot$GeneID)]
res$deg <- !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 2

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.table(res, gzfile(sprintf("results/tables/02_deseq2_%s.tsv.gz", tag)), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(adj, gzfile(sprintf("data/processed/combat_counts_%s.tsv.gz", tag)), sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf("[%s] 시료 %d, 유전자 %d, DEG %d (상향 %d, 하향 %d)\n", tag, ncol(adj), nrow(adj), sum(res$deg),
            sum(res$deg & res$log2FoldChange > 0), sum(res$deg & res$log2FoldChange < 0)))
