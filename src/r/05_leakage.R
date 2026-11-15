# Stage 5 ②: 유전자 고르기에 검증 11명이 들어간 효과(누수)를 본다.
# 사전 명세(2026-11-09 고정):
#   분할 = Stage 3·4 와 같은 시드 1~100 (data/processed/04_splits_100.tsv).
#   누수 없는 선택 = 학습 20명만으로 논문 설정(ComBat-Seq 군 공변량 + protein-coding → DESeq2, padj < 0.05 & |log2FC| ≥ 2).
#   누수 있는 선택 = 논문 DEG 210개(31명 전체로 고름).
#   평가 = 검증 11명의 단일 유전자 AUC. 검증 값은 라벨을 쓰지 않는 log2(정규화 count + 1)(크기 인자는 31명 전체, 배치 보정 없음)이고
#          두 선택에 똑같이 쓴다 → 차이는 "고를 때 검증 11명을 썼나" 하나만 남는다. 방향은 각 선택의 log2FC 부호.
# 실행: Rscript src/r/05_leakage.R  (프로젝트 폴더에서, 8코어)
suppressPackageStartupMessages({library(sva); library(DESeq2); library(readxl); library(parallel)})

annot <- read.delim(gzfile("data/raw/Human.GRCh38.p13.annot.tsv.gz"), quote = "")[, c("GeneID", "Symbol", "GeneType")]
samples <- read.delim("data/processed/samples.tsv")
counts <- as.matrix(read.delim(gzfile("data/processed/counts_main.tsv.gz"), row.names = 1, check.names = FALSE))
counts <- counts[rownames(counts) %in% annot$GeneID[annot$GeneType == "protein-coding"], ]
meta <- samples[match(colnames(counts), samples$gsm), ]
meta$condition <- factor(meta$paper_group, levels = c("OH", "OT2D"))
sym <- annot$Symbol[match(as.integer(rownames(counts)), annot$GeneID)]

xl <- "paper/41598_2025_29141_MOESM2_ESM.xlsx"
paper <- rbind(as.data.frame(read_excel(xl, "Table S1", skip = 2)), as.data.frame(read_excel(xl, "Table S2", skip = 2)))
paper_dir <- setNames(sign(paper$log2FoldChange), paper[[1]])

# 검증 값: 라벨 없이 31명 전체로 크기 인자만 구한 log2 정규화 count
sf <- estimateSizeFactorsForMatrix(counts)
logn <- log2(t(t(counts) / sf) + 1)
rownames(logn) <- sym

auc_dir <- function(v, yy, d) {  # 방향 d(+1: OT2D 가 높음)로 고정한 AUC = Mann-Whitney
  v <- d * v; r <- rank(v); n1 <- sum(yy == "OT2D"); n0 <- sum(yy == "OH")
  (sum(r[yy == "OT2D"]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

splits <- read.delim("data/processed/04_splits_100.tsv")
one_seed <- function(s) {
  sp <- splits[splits$seed == s, ]
  tr <- match(sp$gsm[sp$set == "train"], colnames(counts)); te <- match(sp$gsm[sp$set == "test"], colnames(counts))
  bt <- table(meta$series[tr])
  if (any(bt < 2)) return(data.frame(seed = s, note = "학습에 1명뿐인 배치", n_deg = NA, overlap = NA, aim2 = NA, fhod3 = NA,
                                     auc_leaky = NA, auc_clean = NA))
  adj <- suppressMessages(ComBat_seq(counts[, tr], batch = meta$series[tr], group = meta$condition[tr]))
  dds <- DESeq(DESeqDataSetFromMatrix(adj, meta[tr, ], ~condition), quiet = TRUE)
  res <- results(dds, contrast = c("condition", "OT2D", "OH"))
  deg <- which(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 2)
  yy <- meta$condition[te]
  a_clean <- sapply(deg, function(i) auc_dir(logn[i, te], yy, sign(res$log2FoldChange[i])))
  a_leaky <- sapply(names(paper_dir), function(g) auc_dir(logn[g, te], yy, paper_dir[[g]]))
  data.frame(seed = s, note = "", n_deg = length(deg), overlap = sum(sym[deg] %in% names(paper_dir)),
             aim2 = "AIM2" %in% sym[deg], fhod3 = "FHOD3" %in% sym[deg],
             auc_leaky = median(a_leaky), auc_clean = if (length(deg)) median(a_clean) else NA)
}
out <- do.call(rbind, mclapply(1:100, one_seed, mc.cores = 8))
write.table(out, "results/tables/05_leakage_100seeds.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
ok <- out[out$note == "", ]
cat(sprintf("실행된 시드 %d / 100 (건너뜀 %d: 학습에 1명뿐인 배치)\n", nrow(ok), sum(out$note != "")))
cat(sprintf("학습 20명만으로 고른 DEG 수: 중앙 %g (범위 %g–%g), 논문 210개와 겹침 중앙 %g\n",
            median(ok$n_deg), min(ok$n_deg), max(ok$n_deg), median(ok$overlap)))
cat(sprintf("AIM2 포함 %d번, FHOD3 포함 %d번\n", sum(ok$aim2), sum(ok$fhod3)))
cat(sprintf("검증 AUC 중앙값의 중앙: 누수 있음(논문 210) %.3f / 누수 없음 %.3f\n",
            median(ok$auc_leaky), median(ok$auc_clean, na.rm = TRUE)))
