# 원저자 답장 점검 ② 세부: 09_reply_author_code.R 에서 논문 SVM sigma(0.004400041)가 그대로 나온 설정
#   (입력 vst, 시료 순서 GSE175495 → GSE166047 → GSE141432, 시리즈 안은 GSM 오름차순)을 다시 돌려
#   논문 표 2 · 그림 값과 하나씩 맞춰 본다: LASSO λ 경로와 논문 λmin(0.01307997)에서의 유전자, SVM-RFE 크기별 정확도,
#   RF 상위 10, 검증 시료 구성과 AIM2 · FHOD3 단일 유전자 검증 AUC(논문 0.933 · 1.000).
# 실행: Rscript src/r/09_reply_author_split_detail.R  (프로젝트 폴더에서)
suppressPackageStartupMessages({library(dplyr); library(glmnet); library(caret); library(kernlab); library(randomForest); library(pROC)})

pub <- list(
  LASSO = c("CD52", "PROK1", "AIM2", "CHI3L1", "UGT1A1", "PNPLA1", "STMN2", "CLEC12A", "PPP4R4", "CCL18", "FHOD3"),
  SVM = c("AIM2", "AKR1B15", "ALAS2", "ALPK3", "AQP9", "BCL2A1", "CA1", "CD52", "CHI3L1", "CR1L", "FADS1", "FASN", "FHOD3",
          "GLT1D1", "GPAT3", "IFIT1B", "LRRC3C", "MCEMP1", "MEIKIN", "NFE2", "PHLDA2", "PLEK2", "PNPLA1"),
  RF = c("ALPK3", "MCEMP1", "SLC4A1", "ALAS2", "PHLDA2", "FADS1", "FHOD3", "AIM2", "RUNDC3A", "GPAT3"))

d <- read.delim(gzfile("data/processed/ml_input_vst.tsv.gz"), check.names = FALSE)
s <- ifelse(grepl("^GSM53354", d$gsm), 1, ifelse(grepl("^GSM50607", d$gsm), 2, 3))  # GSE175495, GSE166047, GSE141432
d <- d[order(s, as.integer(sub("GSM", "", d$gsm))), ]
genes <- colnames(d)[-(1:2)]
sig_gene_matrix <- data.frame(Group = factor(d$condition), row = d$gsm, d[, genes], check.names = FALSE)

set.seed(100)
training <- sig_gene_matrix %>% group_by(Group) %>% sample_frac(size = 0.65, replace = F)
testing <- anti_join(sig_gene_matrix, training, by = "row")
training_x <- as.matrix(training[, genes])
y01 <- ifelse(training$Group == "OH", 0, 1)

set.seed(100)
cvfit <- cv.glmnet(training_x, y01, alpha = 1, family = "binomial", nfolds = 10)
lam <- cvfit$lambda
k_min <- which(lam == cvfit$lambda.min); k_pub <- which.min(abs(lam - 0.01307997))
at <- function(k) { b <- as.matrix(coef(cvfit$glmnet.fit, s = lam[k]))[-1, 1]; names(b)[b != 0] }
lasso_pub <- at(k_pub)

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
imp <- as.data.frame(importance(rf_model)); rf10 <- rownames(imp)[order(-imp$MeanDecreaseGini)][1:10]

auc1 <- function(g, dir) as.numeric(auc(roc(testing$Group, testing[[g]], levels = c("OH", "OT2D"), direction = dir, quiet = TRUE)))  # AIM2 상향 "<", FHOD3 하향 ">"
res <- data.frame(item = c("train_OH", "train_OT2D", "test_OH", "test_OT2D", "test_gsm",
                           "lambda_min", "lambda_pub_on_path", "lambda_pub_rel_diff", "lambda_steps_min_to_pub",
                           "lasso_at_min", "lasso_at_pub", "lasso_at_pub_vs_published",
                           "svm_sigma", "svm_C", "rfe_best_size", "rfe_acc_best", "rfe_acc_23", "rfe_genes_vs_published",
                           "rf10", "rf10_vs_published", "hub_at_min", "hub_with_lasso_at_pub",
                           "test_auc_AIM2", "test_auc_FHOD3"),
                  value = c(sum(training$Group == "OH"), sum(training$Group == "OT2D"), sum(testing$Group == "OH"), sum(testing$Group == "OT2D"),
                            paste(testing$row, collapse = ","),
                            signif(cvfit$lambda.min, 8), signif(lam[k_pub], 8), signif(lam[k_pub] / 0.01307997 - 1, 3), k_pub - k_min,
                            paste(at(k_min), collapse = ","), paste(lasso_pub, collapse = ","),
                            sprintf("%d of %d published; extra %s", length(intersect(lasso_pub, pub$LASSO)), length(pub$LASSO),
                                    paste(setdiff(lasso_pub, pub$LASSO), collapse = ",")),
                            signif(svm_model$bestTune$sigma, 10), svm_model$bestTune$C, rfe.result$optsize,
                            round(max(rfe.result$results$Accuracy), 4), round(rfe.result$results$Accuracy[rfe.result$results$Variables == 23], 4),
                            sprintf("%d of %d published in our %d", length(intersect(predictors(rfe.result), pub$SVM)), length(pub$SVM), length(predictors(rfe.result))),
                            paste(rf10, collapse = ","), sprintf("%d of 10", length(intersect(rf10, pub$RF))),
                            paste(Reduce(intersect, list(at(k_min), predictors(rfe.result), rf10)), collapse = ","),
                            paste(Reduce(intersect, list(lasso_pub, predictors(rfe.result), rf10)), collapse = ","),
                            round(auc1("AIM2", "<"), 3), round(auc1("FHOD3", ">"), 3)))
print(res, right = FALSE)
cat("\nRFE 크기별 정확도\n"); print(rfe.result$results[, c("Variables", "Accuracy", "Kappa")], digits = 4)
res <- rbind(res, data.frame(item = "generated", value = format(Sys.time(), "%Y-%m-%d %H:%M %Z")))
write.table(res, "results/tables/09_reply_author_split_detail.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
