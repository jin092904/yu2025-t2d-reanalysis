# 유전자별로 "분할 100번 중 검증 AUC = 1 인 횟수"를 센다(04_roc.R 와 같은 입력·분할·방향, 새 분석 아님).
# 04_test_auc_100seeds_R.tsv 는 분할별 요약(AUC 1 인 DEG 비율)만 있어서, AIM2 82번 · FHOD3 100번과 같은 단위로 비교하려고 만든다.
# 출력: results/tables/07_gene_auc1_splits.tsv (유전자, AUC 1 분할 수, 31명 전체 AUC)
# 실행: Rscript src/r/07_gene_auc_100splits.R  (프로젝트 폴더에서)
suppressPackageStartupMessages(library(pROC))
d <- read.delim(gzfile("data/processed/ml_input_vst.tsv.gz"), check.names = FALSE)
y <- factor(d$condition, levels = c("OH", "OT2D")); x <- as.matrix(d[, -(1:2)])
lfc <- read.delim("results/tables/02_paper_vs_R_group_pc.tsv", check.names = FALSE)
up <- setNames(lfc$log2FoldChange > 0, lfc$Symbol)[colnames(x)]
sp <- read.delim("data/processed/04_splits_100.tsv")
stopifnot(all(sp$gsm[sp$seed == 1] == d$gsm))
a <- function(g, idx) as.numeric(auc(roc(y[idx], x[idx, g], levels = c("OH", "OT2D"), direction = if (up[g]) "<" else ">", quiet = TRUE)))
n1 <- setNames(numeric(ncol(x)), colnames(x))
for (s in 1:100) {
  te <- which(sp$set[sp$seed == s] == "test")
  n1 <- n1 + sapply(colnames(x), function(g) a(g, te) == 1)
}
out <- data.frame(gene = colnames(x), splits_auc1 = n1, auc_all31 = sapply(colnames(x), a, idx = seq_along(y)))
out <- out[order(-out$splits_auc1, out$gene), ]
seeds <- read.delim("results/tables/04_test_auc_100seeds_R.tsv")
stopifnot(out["AIM2", "splits_auc1"] == sum(seeds$AIM2 == 1), out["FHOD3", "splits_auc1"] == sum(seeds$FHOD3 == 1))  # 04 결과와 같은지
write.table(out, "results/tables/07_gene_auc1_splits.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("100번 모두 AUC 1:", sum(out$splits_auc1 == 100), "개 —", paste(out$gene[out$splits_auc1 == 100], collapse = ", "), "\n")
cat("82번 이상:", sum(out$splits_auc1 >= 82), "개 / 유전자별 중앙값:", median(out$splits_auc1), "번 / 31명 전체 AUC 1:", sum(out$auc_all31 == 1), "개\n")
