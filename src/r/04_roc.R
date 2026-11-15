# Stage 4 (R): 논문 허브 유전자(AIM2·FHOD3)의 학습/검증 ROC AUC 와 95% CI(pROC, DeLong).
# 사전 명세(2026-11-09 고정): Stage 3 사전 명세와 같은 입력(vst)·같은 분할(set.seed(100); sample(31, 20)).
#   단독 AUC = 유전자 값 그대로 roc(), 방향은 논문 결과(AIM2 상향, FHOD3 하향)로 고정.
#   두 유전자 모형 = 학습 자료 glm(이항, AIM2 + FHOD3) 의 예측확률.
# 이어서 같은 절차를 Stage 3 와 같은 시드 1~100 으로 반복해 검증 AUC 분포를 보고, DEG 210개 전체의 검증 AUC 와 비교한다.
# 실행: Rscript src/r/04_roc.R  (프로젝트 폴더에서)
suppressPackageStartupMessages(library(pROC))

d <- read.delim(gzfile("data/processed/ml_input_vst.tsv.gz"), check.names = FALSE)
y <- factor(d$condition, levels = c("OH", "OT2D")); x <- as.matrix(d[, -(1:2)])
lfc <- read.delim("results/tables/02_paper_vs_R_group_pc.tsv", check.names = FALSE)
up <- setNames(lfc$log2FoldChange > 0, lfc$Symbol)[colnames(x)]  # DEG 방향(논문 표 S1·S2)
stopifnot(!anyNA(up), up["AIM2"], !up["FHOD3"])

auc_ci <- function(score, yy, case_high = TRUE) {
  r <- roc(yy, score, levels = c("OH", "OT2D"), direction = if (case_high) "<" else ">", quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  c(auc = ci[2], lo = max(ci[1], 0), hi = min(ci[3], 1))
}

one_split <- function(tr) {
  te <- setdiff(seq_along(y), tr)
  fit <- suppressWarnings(glm(y[tr] ~ AIM2 + FHOD3, family = binomial, data = as.data.frame(x[tr, ])))
  p <- predict(fit, newdata = as.data.frame(x), type = "response")
  rbind(cbind(set = "train", model = c("AIM2", "FHOD3", "AIM2+FHOD3"),
              rbind(auc_ci(x[tr, "AIM2"], y[tr]), auc_ci(x[tr, "FHOD3"], y[tr], FALSE), auc_ci(p[tr], y[tr]))),
        cbind(set = "test", model = c("AIM2", "FHOD3", "AIM2+FHOD3"),
              rbind(auc_ci(x[te, "AIM2"], y[te]), auc_ci(x[te, "FHOD3"], y[te], FALSE), auc_ci(p[te], y[te]))))
}

# 1) 사전 명세 분할
set.seed(100)
tr <- sort(sample(length(y), 20))
main <- as.data.frame(one_split(tr))
main$n_OH <- ifelse(main$set == "train", sum(y[tr] == "OH"), sum(y[-tr] == "OH"))
main$n_OT2D <- ifelse(main$set == "train", sum(y[tr] == "OT2D"), sum(y[-tr] == "OT2D"))
main$paper_auc <- c(1, 1, 1, 0.933, 1, 1)
write.table(main, "results/tables/04_roc_R_main.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(main)

# 2) 시드 1~100: 허브 유전자와 DEG 210개 전체의 검증 AUC (방향 고정)
splits <- list(); rows <- list()
for (s in 1:100) {
  set.seed(s)
  tr <- sort(sample(length(y), 20)); te <- setdiff(seq_along(y), tr)
  splits[[s]] <- data.frame(seed = s, gsm = d$gsm, set = ifelse(seq_along(y) %in% tr, "train", "test"))
  auc_all <- sapply(colnames(x), function(g) auc(roc(y[te], x[te, g], levels = c("OH", "OT2D"),
                                                        direction = if (up[g]) "<" else ">", quiet = TRUE)))
  rows[[s]] <- data.frame(seed = s, n_OH_test = sum(y[te] == "OH"), n_OT2D_test = sum(y[te] == "OT2D"),
                          AIM2 = auc_all["AIM2"], FHOD3 = auc_all["FHOD3"], deg_median = median(auc_all),
                          deg_ge09 = mean(auc_all >= 0.9), deg_eq1 = mean(auc_all == 1))
}
write.table(do.call(rbind, splits), "data/processed/04_splits_100.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
seeds <- do.call(rbind, rows)
write.table(seeds, "results/tables/04_test_auc_100seeds_R.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n시드 100개 검증 AUC 중앙값: AIM2", median(seeds$AIM2), "/ FHOD3", median(seeds$FHOD3),
    "/ DEG 210개 중앙값의 중앙값", median(seeds$deg_median), "\n")
cat("검증 AUC >= 0.9 인 DEG 비율(시드별 중앙값):", median(seeds$deg_ge09), "\n")
