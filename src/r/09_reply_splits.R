# 원저자 답장의 저자 코드 · 저자 분할 방식으로 분할 안정성을 본다(그림 2b · 2c 의 입력).
# 저자 분할 = 09_reply_author_split_detail.R 에서 찾은 설정: 입력 vst, 시료 순서 GSE175495 → GSE166047 → GSE141432
#   (시리즈 안 GSM 오름차순), set.seed(시드) → group_by(Group) %>% sample_frac(0.65). 시드 100 이 논문 분할이다.
# 시드 1~100 마다 분할 시드만 바꾸고 나머지 단계는 저자 코드처럼 set.seed(100) 뒤에 돌린다(분할 효과만 남긴다).
# 출력
#   data/processed/09_splits_authors_100.tsv        시드별 학습/검증 GSM (05_leakage.R splits= 로 넘긴다)
#   results/tables/09_stability_authors_100.tsv     시드별 LASSO · SVM-RFE · RF 상위 10 · 세 방법 공통
#   results/tables/09_auc1_authors_100.tsv          논문 DEG 210개 각각이 검증 AUC = 1 인 분할 수(방향 = 논문 log2FC 부호)
#   results/tables/09_auc1_per_split_authors_100.tsv 분할마다 검증 AUC = 1 인 DEG 수
# 실행: Rscript src/r/09_reply_splits.R  (프로젝트 폴더에서, 16코어, 약 10분)
suppressPackageStartupMessages({library(dplyr); library(glmnet); library(caret); library(kernlab); library(randomForest); library(pROC)})

d <- read.delim(gzfile("data/processed/ml_input_vst.tsv.gz"), check.names = FALSE)
s <- ifelse(grepl("^GSM53354", d$gsm), 1, ifelse(grepl("^GSM50607", d$gsm), 2, 3))  # GSE175495, GSE166047, GSE141432
d <- d[order(s, as.integer(sub("GSM", "", d$gsm))), ]
genes <- colnames(d)[-(1:2)]
sig_gene_matrix <- data.frame(Group = factor(d$condition), row = d$gsm, d[, genes], check.names = FALSE)
lfc <- read.delim("results/tables/02_paper_vs_R_group_pc.tsv", check.names = FALSE)
up <- setNames(lfc$log2FoldChange > 0, lfc$Symbol)[genes]

author_split <- function(seed) {
  set.seed(seed)
  training <- sig_gene_matrix %>% group_by(Group) %>% sample_frac(size = 0.65, replace = F)
  list(training = training, testing = anti_join(sig_gene_matrix, training, by = "row"))
}

author_select <- function(training) {  # 답장 코드 그대로(분할 뒤 단계)
  training_x <- as.matrix(training[, genes])
  y01 <- ifelse(training$Group == "OH", 0, 1)
  set.seed(100)
  cvfit <- cv.glmnet(training_x, y01, alpha = 1, family = "binomial", nfolds = 10)
  b <- coef(cvfit, s = cvfit$lambda.min)
  lasso <- setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
  training_y <- droplevels(training$Group)
  set.seed(100)
  svm_model <- train(training_x, training_y, method = "svmRadial", trControl = trainControl(method = "cv", number = 10), tuneLength = 10)
  set.seed(100)
  control <- rfeControl(functions = caretFuncs, method = "cv", number = 10)
  set.seed(100)
  rfe.result <- rfe(training_x, training_y, sizes = c(3 * (1:10), 23, 10 * (4:20)), rfeControl = control, method = "svmRadial",
                    tuneGrid = expand.grid(sigma = svm_model$bestTune$sigma, C = svm_model$bestTune$C))
  set.seed(100)
  rf_model <- randomForest(x = training_x, y = as.factor(training_y), importance = TRUE, proximity = TRUE)
  imp <- as.data.frame(importance(rf_model))
  rf10 <- rownames(imp)[order(-imp$MeanDecreaseGini)][1:10]
  list(LASSO = lasso, `SVM-RFE` = predictors(rfe.result), RF = rf10, hub = Reduce(intersect, list(lasso, predictors(rfe.result), rf10)))
}

auc1 <- function(testing) sapply(genes, function(g)
  as.numeric(auc(roc(testing$Group, testing[[g]], levels = c("OH", "OT2D"), direction = if (up[[g]]) "<" else ">", quiet = TRUE))) == 1)

runs <- parallel::mclapply(1:100, function(seed) {
  sp <- author_split(seed)
  list(split = rbind(data.frame(seed = seed, gsm = sp$training$row, set = "train"), data.frame(seed = seed, gsm = sp$testing$row, set = "test")),
       sel = author_select(sp$training), auc1 = auc1(sp$testing))
}, mc.cores = 16)

splits <- do.call(rbind, lapply(runs, `[[`, "split"))
write.table(splits, "data/processed/09_splits_authors_100.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
stab <- do.call(rbind, lapply(1:100, function(s) do.call(rbind, lapply(names(runs[[s]]$sel), function(m)
  data.frame(seed = s, method = m, n = length(runs[[s]]$sel[[m]]), genes = paste(runs[[s]]$sel[[m]], collapse = ","))))))
write.table(stab, "results/tables/09_stability_authors_100.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
a1 <- Reduce(`+`, lapply(runs, `[[`, "auc1"))
per_split <- sapply(runs, function(r) sum(r$auc1))
write.table(data.frame(seed = 1:100, n_auc1 = per_split, AIM2 = sapply(runs, function(r) r$auc1[["AIM2"]]),
                       FHOD3 = sapply(runs, function(r) r$auc1[["FHOD3"]])), "results/tables/09_auc1_per_split_authors_100.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = genes, splits_auc1 = a1)[order(-a1, genes), ], "results/tables/09_auc1_authors_100.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

hub <- strsplit(stab$genes[stab$method == "hub"], ",")
freq <- sort(table(unlist(hub)), decreasing = TRUE)
both <- sum(sapply(hub, function(h) all(c("AIM2", "FHOD3") %in% h)))
cat(sprintf("시드 100(논문 분할) 공통: %s\n", stab$genes[stab$method == "hub" & stab$seed == 100]))
cat(sprintf("세 방법 공통 빈도 상위: %s\n", paste(sprintf("%s %d", names(freq)[1:8], freq[1:8]), collapse = ", ")))
cat(sprintf("AIM2 %d번, FHOD3 %d번, 둘 다 %d번, 정확히 AIM2 · FHOD3 둘만 %d번, 공통 없음 %d번\n",
            sum(sapply(hub, function(h) "AIM2" %in% h)), sum(sapply(hub, function(h) "FHOD3" %in% h)), both,
            sum(sapply(hub, function(h) setequal(h, c("AIM2", "FHOD3")))), sum(lengths(hub) == 0 | sapply(hub, function(h) identical(h, "")))))
cat(sprintf("검증 AUC = 1 인 DEG 수(분할별): 중앙 %g (%.0f%%), 범위 %g–%g / AIM2 %d번, FHOD3 %d번, 100번 모두 1: %s\n",
            median(per_split), 100 * median(per_split) / length(genes), min(per_split), max(per_split), a1[["AIM2"]], a1[["FHOD3"]],
            paste(names(a1)[a1 == 100], collapse = ", ")))
