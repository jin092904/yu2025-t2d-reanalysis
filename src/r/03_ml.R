# Stage 3 (R): 학습/검증 분할 → LASSO · SVM-RFE · RF 로 유전자를 고르고 세 방법의 교집합을 구한다.
# 사전 명세(2026-11-09 고정, Methods 문장 그대로):
#   입력 vst(blind = TRUE), 논문 DEG 210개, 시료 순서 = 표 1 순서(GSE175495, GSE141432, GSE166047),
#   set.seed(100) 한 번 → sample(31, 20) 으로 학습 20 / 검증 11 → cv.glmnet(이항, 10겹) λmin
#   → caret::rfe(svmRadial, sigma 0.004400041, C 0.5 고정, 10겹 CV, 크기 1~210) → randomForest(ntree 500) Gini 상위 10.
# 탐색 인자: input=vst|vst_nb|norm|lognorm, order=table1|gsm, split=sample|pergroup, reseed=no|yes(cv.glmnet 직전 set.seed(100))
# "search" 를 주면 LASSO 만 모든 조합으로 돌려 λmin 을 표로 낸다(논문 λmin = 0.01307997).
# "stability" 를 주면 시드 1~100 으로 전체 절차를 반복한다(분할 운 확인, 8코어).
# 실행: Rscript src/r/03_ml.R [search|stability] [input=..] [order=..] [split=..] [reseed=..]  (프로젝트 폴더에서)
suppressPackageStartupMessages({library(glmnet); library(caret); library(randomForest)})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(input = "vst", order = "table1", split = "sample", reseed = "no")
for (a in grep("=", args, value = TRUE)) opt[[sub("=.*", "", a)]] <- sub(".*=", "", a)

read_input <- function(input, order) {
  d <- read.delim(gzfile(sprintf("data/processed/ml_input_%s.tsv.gz", input)), check.names = FALSE)
  if (order == "gsm") d <- d[order(as.integer(sub("GSM", "", d$gsm))), ]
  list(gsm = d$gsm, y = factor(d$condition, levels = c("OH", "OT2D")), x = as.matrix(d[, -(1:2)]))
}

split_train <- function(y, split) {
  if (split == "sample") return(sort(sample(length(y), 20)))  # Methods: 무작위 65% (n = 20)
  sort(unlist(lapply(levels(y), function(g) {                  # Results: 군별 65%
    i <- which(y == g); i[sample(length(i), round(0.65 * length(i)))]
  })))
}

run_lasso <- function(x, y, reseed) {
  if (reseed == "yes") set.seed(100)
  cv <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10, keep = TRUE)
  b <- as.matrix(coef(cv, s = "lambda.min"))[-1, 1]
  list(cv = cv, genes = names(b)[b != 0], coef = b[b != 0])
}

if ("search" %in% args) {
  grid <- expand.grid(input = c("vst", "vst_nb", "norm", "lognorm"), order = c("table1", "gsm"),
                      split = c("sample", "pergroup"), reseed = c("no", "yes"), stringsAsFactors = FALSE)
  grid$lambda_min <- NA; grid$n_genes <- NA; grid$n_train <- NA
  for (k in seq_len(nrow(grid))) {
    d <- read_input(grid$input[k], grid$order[k])
    set.seed(100)
    tr <- split_train(d$y, grid$split[k])
    l <- run_lasso(d$x[tr, ], d$y[tr], grid$reseed[k])
    grid$lambda_min[k] <- l$cv$lambda.min; grid$n_genes[k] <- length(l$genes); grid$n_train[k] <- length(tr)
  }
  grid$match <- abs(grid$lambda_min - 0.01307997) < 5e-9
  write.table(grid, "results/tables/03_lambda_search.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  print(grid[order(abs(grid$lambda_min - 0.01307997)), ], digits = 8, row.names = FALSE)
  quit(save = "no")
}

pipeline <- function(d, seed) {
  set.seed(seed)
  tr <- split_train(d$y, opt$split)
  x <- d$x[tr, ]; y <- d$y[tr]
  # 1) LASSO
  l <- run_lasso(x, y, opt$reseed)
  # 2) SVM-RFE: caret 의 rfe. 순위는 caretFuncs 가 varImp 로 매긴다(svmRadial 은 모형 고유 중요도가 없어 단변량 ROC 면적).
  svm <- rfe(x, y, sizes = 1:ncol(x), rfeControl = rfeControl(functions = caretFuncs, method = "cv", number = 10),
             method = "svmRadial", tuneGrid = data.frame(sigma = 0.004400041, C = 0.5), trControl = trainControl(method = "none"))
  # 3) RF: MeanDecreaseGini 상위 10
  rf <- randomForest(x, y, ntree = 500)
  gini <- sort(importance(rf)[, "MeanDecreaseGini"], decreasing = TRUE)
  genes <- list(LASSO = l$genes, `SVM-RFE` = predictors(svm), RF = names(gini)[1:10])
  genes$hub <- Reduce(intersect, genes)
  list(tr = tr, y = y, l = l, svm = svm, gini = gini, genes = genes)
}

tag <- paste(unlist(opt), collapse = "_")
d <- read_input(opt$input, opt$order)

if ("stability" %in% args) {  # 분할 운: 시드 1~100 에서 같은 절차를 반복
  runs <- parallel::mclapply(1:100, function(s) pipeline(d, s)$genes, mc.cores = 8)
  rows <- do.call(rbind, lapply(1:100, function(s) do.call(rbind, lapply(names(runs[[s]]), function(m)
    data.frame(seed = s, method = m, n = length(runs[[s]][[m]]), genes = paste(runs[[s]][[m]], collapse = ","))))))
  write.table(rows, sprintf("results/tables/03_seed_stability_%s.tsv", tag), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("시드 100개 완료\n")
  quit(save = "no")
}

r <- pipeline(d, 100)
tr <- r$tr; y <- r$y; l <- r$l; svm <- r$svm; gini <- r$gini
svm_genes <- r$genes$`SVM-RFE`; rf_genes <- r$genes$RF; hub <- r$genes$hub
out <- rbind(data.frame(method = "LASSO", rank = seq_along(l$genes), gene = l$genes, value = unname(l$coef)),
             data.frame(method = "SVM-RFE", rank = seq_along(svm_genes), gene = svm_genes, value = NA),
             data.frame(method = "RF", rank = 1:10, gene = rf_genes, value = unname(gini[1:10])),
             data.frame(method = "hub", rank = seq_along(hub), gene = hub, value = NA))
write.table(out, sprintf("results/tables/03_ml_R_%s.tsv", tag), sep = "\t", quote = FALSE, row.names = FALSE)

# Python 대조용: 같은 분할·같은 겹(LASSO foldid, rfe 겹), λ 경로, rfe 크기별 정확도, RF 전체 중요도
fold_rfe <- integer(length(tr))
for (f in seq_along(svm$control$index)) fold_rfe[setdiff(seq_along(tr), svm$control$index[[f]])] <- f
write.table(data.frame(gsm = d$gsm[tr], foldid_lasso = l$cv$foldid, fold_rfe = fold_rfe),
            sprintf("data/processed/03_split_%s.tsv", tag), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(lambda = l$cv$lambda, cvm = l$cv$cvm, nzero = l$cv$nzero),
            sprintf("results/tables/03_lasso_path_R_%s.tsv", tag), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(svm$results, sprintf("results/tables/03_rfe_sizes_R_%s.tsv", tag), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = names(gini), MeanDecreaseGini = gini), sprintf("results/tables/03_rf_gini_R_%s.tsv", tag),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("[%s] 학습 %d (OH %d, OT2D %d)\n", tag, length(tr), sum(y == "OH"), sum(y == "OT2D")))
cat(sprintf("LASSO λmin %.8f (논문 0.01307997), %d개: %s\n", l$cv$lambda.min, length(l$genes), paste(l$genes, collapse = ", ")))
cat(sprintf("SVM-RFE %d개 (논문 23): %s\n", length(svm_genes), paste(svm_genes, collapse = ", ")))
cat(sprintf("RF 상위 10: %s\n", paste(rf_genes, collapse = ", ")))
cat(sprintf("세 방법 공통: %s (논문 AIM2, FHOD3)\n", paste(hub, collapse = ", ")))
