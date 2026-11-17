# 기계학습 입력 행렬(시료 × 논문 DEG 210개)을 만든다.
# 원문에 입력값 종류가 없어 후보 4가지를 만든다. 그림 2D 가 "vst-transformed data" 라서 사전 명세는 vst(기본값 blind = TRUE).
#   vst      : DESeq2 vst(blind = TRUE)      vst_nb : vst(blind = FALSE)
#   norm     : 크기 인자로 정규화한 count     lognorm: log2(정규화 count + 1)
# 보정 count 는 02_combat_deseq2.R 의 재현 설정(ComBat-Seq 군 공변량 + protein-coding)이다.
# 실행: Rscript src/r/03_ml_inputs.R  (프로젝트 폴더에서)
suppressPackageStartupMessages({library(DESeq2); library(readxl)})

adj <- as.matrix(read.delim(gzfile("data/processed/combat_counts_R_group_pc.tsv.gz"), row.names = 1, check.names = FALSE))
samples <- read.delim("data/processed/samples.tsv")
meta <- samples[match(colnames(adj), samples$gsm), ]
meta$condition <- factor(meta$paper_group, levels = c("OH", "OT2D"))
dds <- DESeq(DESeqDataSetFromMatrix(adj, meta, ~condition), quiet = TRUE)

xl <- "paper/41598_2025_29141_MOESM2_ESM.xlsx"
deg <- c(read_excel(xl, "Table S1", skip = 2)[[1]], read_excel(xl, "Table S2", skip = 2)[[1]])
annot <- read.delim(gzfile("data/raw/Human.GRCh38.p13.annot.tsv.gz"), quote = "")[, c("GeneID", "Symbol")]
ids <- as.character(annot$GeneID[match(deg, annot$Symbol)])
stopifnot(length(deg) == 210, !anyNA(ids), all(ids %in% rownames(adj)))

norm <- counts(dds, normalized = TRUE)
inputs <- list(vst = assay(vst(dds, blind = TRUE)), vst_nb = assay(vst(dds, blind = FALSE)),
               norm = norm, lognorm = log2(norm + 1))
for (k in names(inputs)) {
  x <- t(inputs[[k]][ids, ])
  colnames(x) <- deg
  write.table(data.frame(gsm = rownames(x), condition = meta$condition, x, check.names = FALSE),
              gzfile(sprintf("data/processed/ml_input_%s.tsv.gz", k)), sep = "\t", quote = FALSE, row.names = FALSE)
}
cat("입력 4종:", names(inputs), "/ 시료", ncol(adj), "× 유전자", length(ids), "\n")
