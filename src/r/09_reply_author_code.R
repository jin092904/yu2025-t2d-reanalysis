# 원저자 답장(2026-11-17) 점검 ②: 저자가 보내 준 특징 선택 코드를 우리 자료에 글자 그대로 돌린다.
# 답장 코드: set.seed(100) → group_by(Group) %>% sample_frac(0.65) → 단계마다 set.seed(100) 후
#   cv.glmnet(이항, 10겹) λmin / caret train(svmRadial, 10겹 CV, tuneLength 10) → rfe(caretFuncs, 10겹,
#   sizes c(3*(1:10), 23, 10*(4:20)), 고른 sigma · C) / randomForest(importance = TRUE, proximity = TRUE) Gini 상위 10.
# 답장에 없는 것: 입력값 종류("batch-corrected and normalized")와 시료 행 순서. 그래서
#   입력 4종(vst · vst_nb · norm · lognorm, 03_ml_inputs.R) × GEO 시리즈 순서 6가지(시리즈 안은 GSM 오름차순)를 모두 돌린다.
# 저자 결과와 같은 분할 · 입력인지는 논문이 밝힌 값 두 개로 본다: LASSO λmin = 0.01307997, SVM sigma = 0.004400041(C = 0.5).
# 실행: Rscript src/r/09_reply_author_code.R  (프로젝트 폴더에서, 8코어)
suppressPackageStartupMessages({library(dplyr); library(glmnet); library(caret); library(kernlab); library(randomForest)})

series_of <- function(gsm) ifelse(grepl("^GSM53354", gsm), "GSE175495", ifelse(grepl("^GSM42027", gsm), "GSE141432", "GSE166047"))
perms <- list(c("GSE175495", "GSE141432", "GSE166047"), c("GSE175495", "GSE166047", "GSE141432"),
              c("GSE141432", "GSE175495", "GSE166047"), c("GSE141432", "GSE166047", "GSE175495"),
              c("GSE166047", "GSE175495", "GSE141432"), c("GSE166047", "GSE141432", "GSE175495"))
grid <- expand.grid(input = c("vst", "vst_nb", "norm", "lognorm"), order = seq_along(perms), stringsAsFactors = FALSE)

author_code <- function(input, order) {
  d <- read.delim(gzfile(sprintf("data/processed/ml_input_%s.tsv.gz", input)), check.names = FALSE)
  s <- factor(series_of(d$gsm), levels = perms[[order]])
  d <- d[order(s, as.integer(sub("GSM", "", d$gsm))), ]
  genes <- colnames(d)[-(1:2)]
  sig_gene_matrix <- data.frame(Group = factor(d$condition), row = d$gsm, d[, genes], check.names = FALSE)

  # divide Training and test set (답장 코드 그대로)
  set.seed(100)
  training <- sig_gene_matrix %>% group_by(Group) %>% sample_frac(size = 0.65, replace = F)
  testing <- anti_join(sig_gene_matrix, training, by = "row")
  training_x <- as.matrix(training[, genes])
  training_y <- training$Group

  # LASSO reg
  training_y <- ifelse(training_y == "OH", 0, 1)
  set.seed(100)
  fit <- glmnet(training_x, training_y, alpha = 1)
  set.seed(100)
  cvfit <- cv.glmnet(training_x, training_y, alpha = 1, family = "binomial", nfolds = 10)
  coef <- coef(cvfit, s = cvfit$lambda.min)
  lasso <- setdiff(rownames(coef)[coef[, 1] != 0], "(Intercept)")

  # svm-rfe
  training_y <- droplevels(training$Group)
  set.seed(100)
  train_ctrl <- trainControl(method = "cv", number = 10)
  set.seed(100)
  svm_model <- train(training_x, training_y, method = "svmRadial", trControl = train_ctrl, tuneLength = 10)
  best.gamma <- svm_model$bestTune$sigma
  best.cost <- svm_model$bestTune$C
  set.seed(100)
  control <- rfeControl(functions = caretFuncs, method = "cv", number = 10)
  tune.grid <- expand.grid(sigma = best.gamma, C = best.cost)
  set.seed(100)
  rfe.result <- rfe(training_x, training_y, sizes = c(3 * (1:10), 23, 10 * (4:20)), rfeControl = control,
                    method = "svmRadial", tuneGrid = tune.grid)
  rfe_genes <- predictors(rfe.result)

  # randomForest
  set.seed(100)
  rf_model <- randomForest(x = training_x, y = as.factor(training_y), importance = TRUE, proximity = TRUE)
  imp <- as.data.frame(importance(rf_model))
  rf10 <- rownames(imp)[order(-imp$MeanDecreaseGini)][1:10]

  hub <- Reduce(intersect, list(lasso, rfe_genes, rf10))
  data.frame(input = input, order = paste(perms[[order]], collapse = ">"),
             test_OH = sum(testing$Group == "OH"), test_OT2D = sum(testing$Group == "OT2D"),
             lambda_min = cvfit$lambda.min, lambda_match = abs(cvfit$lambda.min - 0.01307997) < 5e-9,
             sigma = best.gamma, C = best.cost, sigma_match = abs(best.gamma - 0.004400041) < 5e-10,
             n_lasso = length(lasso), n_rfe = length(rfe_genes), hub = paste(hub, collapse = ","),
             AIM2_in = paste0(c("L", "S", "R")[c("AIM2" %in% lasso, "AIM2" %in% rfe_genes, "AIM2" %in% rf10)], collapse = ""),
             FHOD3_in = paste0(c("L", "S", "R")[c("FHOD3" %in% lasso, "FHOD3" %in% rfe_genes, "FHOD3" %in% rf10)], collapse = ""),
             lasso = paste(lasso, collapse = ","), rf10 = paste(rf10, collapse = ","))
}

out <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(k) author_code(grid$input[k], grid$order[k]), mc.cores = 8))
out$generated <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
write.table(out, "results/tables/09_reply_author_code.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out[, c("input", "order", "test_OH", "test_OT2D", "lambda_min", "lambda_match", "sigma", "C", "sigma_match",
              "n_lasso", "n_rfe", "hub", "AIM2_in", "FHOD3_in")], digits = 8, row.names = FALSE)
