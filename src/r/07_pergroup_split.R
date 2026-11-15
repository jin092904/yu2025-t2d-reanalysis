# Stage 7: 논문 Results 의 분할(군마다 65%)로 세 방법 전체를 돌린 결과 요약(Methods 의 단순 무작위 분할과 비교용).
# 입력: 03_ml.R split=pergroup 의 결과(시드 100 한 번 + 시드 1~100 반복). 새 계산은 시드 100 분할의 검증 AUC 뿐(04_roc.R 와 같은 방식 · 방향).
#   results/tables/03_ml_R_vst_table1_pergroup_no.tsv, 03_seed_stability_vst_table1_pergroup_no.tsv,
#   data/processed/03_split_vst_table1_pergroup_no.tsv (시드 100 학습 시료)
# 출력: results/tables/07_pergroup_summary.tsv (항목 · 값), 07_pergroup_hub_frequency.tsv (유전자별 세 방법 공통 비율)
# 실행: Rscript src/r/07_pergroup_split.R  (프로젝트 폴더에서)
suppressPackageStartupMessages(library(pROC))

rd <- function(f) read.delim(f, check.names = FALSE, colClasses = "character")
d <- read.delim(gzfile("data/processed/ml_input_vst.tsv.gz"), check.names = FALSE)
y <- factor(d$condition, levels = c("OH", "OT2D")); x <- as.matrix(d[, -(1:2)])

# 1) 시드 100 한 번: 세 방법이 고른 유전자, 검증 구성, 허브 유전자의 검증 AUC
one <- rd("results/tables/03_ml_R_vst_table1_pergroup_no.tsv")
tr_gsm <- rd("data/processed/03_split_vst_table1_pergroup_no.tsv")$gsm
tr <- which(d$gsm %in% tr_gsm); te <- setdiff(seq_along(y), tr)
auc1 <- function(g, high) as.numeric(auc(roc(y[te], x[te, g], levels = c("OH", "OT2D"), direction = if (high) "<" else ">", quiet = TRUE)))
hub <- one$gene[one$method == "hub"]
sm <- data.frame(item = c("seed100_n_LASSO", "seed100_n_SVM-RFE", "seed100_n_RF", "seed100_hub", "seed100_test_OH", "seed100_test_OT2D",
                          "seed100_test_AUC_AIM2", "seed100_test_AUC_FHOD3"),
                 value = c(sum(one$method == "LASSO"), sum(one$method == "SVM-RFE"), sum(one$method == "RF"),
                           paste(hub, collapse = ","), sum(y[te] == "OH"), sum(y[te] == "OT2D"),
                           sprintf("%.3f", auc1("AIM2", TRUE)), sprintf("%.3f", auc1("FHOD3", FALSE))))

# 2) 시드 1~100: 세 방법 공통 유전자의 빈도, 두 유전자가 함께 · 정확히 둘만 고른 분할 수
st <- rd("results/tables/03_seed_stability_vst_table1_pergroup_no.tsv")
hubs <- lapply(split(st$genes[st$method == "hub"], st$seed[st$method == "hub"]), function(g) if (g == "") character(0) else strsplit(g, ",")[[1]])
stopifnot(length(hubs) == 100)
freq <- sort(table(unlist(hubs)) / 100, decreasing = TRUE)
both <- sum(sapply(hubs, function(h) all(c("AIM2", "FHOD3") %in% h)))
exact <- sum(sapply(hubs, function(h) setequal(h, c("AIM2", "FHOD3"))))
empty <- sum(sapply(hubs, length) == 0)
sm <- rbind(sm, data.frame(item = c("seeds_both_AIM2_FHOD3", "seeds_exactly_AIM2_FHOD3", "seeds_no_common_gene",
                                    "top_gene", "top_gene_pct", "FHOD3_pct", "AIM2_pct"),
                           value = c(both, exact, empty, names(freq)[1], round(100 * freq[[1]]),
                                     round(100 * ifelse(is.na(freq["FHOD3"]), 0, freq[["FHOD3"]])),
                                     round(100 * ifelse(is.na(freq["AIM2"]), 0, freq[["AIM2"]])))))
write.table(sm, "results/tables/07_pergroup_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = names(freq), hub = as.numeric(freq)), "results/tables/07_pergroup_hub_frequency.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
print(sm, row.names = FALSE)
print(head(freq, 10))
