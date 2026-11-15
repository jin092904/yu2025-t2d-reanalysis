# Stage 3 탐색: 논문 λmin(0.01307997)이 후보 학습 자료의 glmnet λ 경로 위에 정확히 놓이는지 본다.
# λ 경로는 자료로만 정해지고(난수 없음) λ_k = λmax × 0.01^((k−1)/99) 이므로, k 가 정수에 가까울수록 같은 학습 자료일 가능성이 크다.
# 후보: 입력 4종 × 시료 순서 4종(표 1, GSM, 군별 정렬 2종) × 분할 2종.
# 실행: Rscript src/r/03_lambda_path_check.R  (프로젝트 폴더에서)
suppressPackageStartupMessages(library(glmnet))
target <- 0.01307997
onpath <- function(x, y) { f <- glmnet(x, y, family = "binomial", alpha = 1); lm <- max(f$lambda)
  k <- 1 + 99 * log(target / lm) / log(0.01); c(lmax = lm, k = k) }
res <- list()
for (input in c("vst","vst_nb","norm","lognorm")) {
  d <- read.delim(gzfile(sprintf("data/processed/ml_input_%s.tsv.gz", input)), check.names = FALSE)
  y <- factor(d$condition, levels = c("OH","OT2D")); x <- as.matrix(d[, -(1:2)])
  ords <- list(table1 = seq_len(nrow(d)), gsm = order(as.integer(sub("GSM","",d$gsm))),
               grp_table1 = order(y, seq_len(nrow(d))), grp_gsm = order(y, as.integer(sub("GSM","",d$gsm))))
  for (o in names(ords)) for (sp in c("sample","pergroup")) {
    xo <- x[ords[[o]], ]; yo <- y[ords[[o]]]
    set.seed(100)
    tr <- if (sp == "sample") sample(31, 20) else unlist(lapply(levels(yo), function(g){i <- which(yo==g); i[sample(length(i), round(.65*length(i)))]}))
    r <- onpath(xo[tr, ], yo[tr]); res[[length(res)+1]] <- data.frame(input, order = o, split = sp, lmax = r["lmax"], k = r["k"])
  }
}
r <- do.call(rbind, res); r$frac <- abs(r$k - round(r$k))
write.table(r[order(r$frac), ], "results/tables/03_lambda_path_check.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(r[order(r$frac), ][1:6, ], digits = 6, row.names = FALSE)
