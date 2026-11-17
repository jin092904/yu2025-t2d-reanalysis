# 한 코호트(= 한 배치) 안에서 비만과 당뇨를 나눠 본다. 배치 보정이 필요 없다.
# 사전 명세(2026-11-09 고정):
#   GSE141432: Lean 9 · Obese NGT 8 · Obese T2D 8 (논문은 Lean 을 OH, Obese T2D 를 OT2D 로 썼고 Obese NGT 는 안 씀)
#   GSE166047 여성: NW 5 · OBF 5 · OBT2D 5 (논문은 OBT2D 만 씀. OBM 은 남성이라 뺌)
#   protein-coding 유전자, 코호트별 DESeq2(~group, 기본값, LFC 축소 없음).
#   비교 = 비만(비만·비당뇨 vs 정상체중), 당뇨(비만·당뇨 vs 비만·비당뇨), 논문형(비만·당뇨 vs 정상체중).
#   AIM2·FHOD3 는 정규화 count 로 정확 Wilcoxon 순위합 검정(명목 P) 도 낸다.
# 실행: Rscript src/r/05_within_cohort.R  (프로젝트 폴더에서)
suppressPackageStartupMessages({library(DESeq2); library(readxl)})

annot <- read.delim(gzfile("data/raw/Human.GRCh38.p13.annot.tsv.gz"), quote = "")[, c("GeneID", "Symbol", "GeneType")]
pc <- annot$GeneID[annot$GeneType == "protein-coding"]
samples <- read.delim("data/processed/samples.tsv")
xl <- "paper/41598_2025_29141_MOESM2_ESM.xlsx"
paper <- rbind(as.data.frame(read_excel(xl, "Table S1", skip = 2)), as.data.frame(read_excel(xl, "Table S2", skip = 2)))
names(paper)[1] <- "Symbol"

cohorts <- list(
  GSE141432 = c(lean = "Lean", obese = "Obese NGT", t2d = "Obese T2D"),
  GSE166047 = c(lean = "NW", obese = "OBF", t2d = "OBT2D"))
CMPS <- list(c("obesity", "obese", "lean"), c("diabetes", "t2d", "obese"), c("paper_like", "t2d", "lean"))   # 비교 이름, 비교군, 기준군
de_rows <- list(); expr_rows <- list(); test_rows <- list()
for (g in names(cohorts)) {
  lab <- cohorts[[g]]
  raw <- read.delim(gzfile(sprintf("data/raw/%s_raw_counts_GRCh38.p13_NCBI.tsv.gz", g)), row.names = 1, check.names = FALSE)
  s <- samples[samples$series == g & samples$raw_group %in% lab, ]
  cnt <- as.matrix(raw[rownames(raw) %in% pc, s$gsm])
  s$group <- factor(names(lab)[match(s$raw_group, lab)], levels = c("lean", "obese", "t2d"))
  dds <- DESeq(DESeqDataSetFromMatrix(cnt, s, ~group), quiet = TRUE)
  sym <- annot$Symbol[match(as.integer(rownames(dds)), annot$GeneID)]
  for (cmp in CMPS) {
    r <- as.data.frame(results(dds, contrast = c("group", cmp[2], cmp[3])))
    r$Symbol <- sym
    r <- r[r$Symbol %in% paper$Symbol, c("Symbol", "baseMean", "log2FoldChange", "pvalue", "padj")]
    de_rows[[length(de_rows) + 1]] <- cbind(cohort = g, comparison = cmp[1], r)
  }
  nc <- counts(dds, normalized = TRUE)
  for (gene in c("AIM2", "FHOD3")) {
    v <- log2(nc[which(sym == gene), ] + 1)
    expr_rows[[length(expr_rows) + 1]] <- data.frame(cohort = g, gsm = s$gsm, group = s$group, gene = gene, log2_norm = v)
    for (cmp in CMPS) {
      a <- v[s$group == cmp[2]]; b <- v[s$group == cmp[3]]
      w <- wilcox.test(a, b, exact = TRUE)
      test_rows[[length(test_rows) + 1]] <- data.frame(cohort = g, gene = gene, comparison = cmp[1],
        n1 = length(a), n2 = length(b), median_diff = median(a) - median(b), auc = unname(w$statistic) / (length(a) * length(b)),
        p_wilcox = w$p.value)
    }
  }
}
de <- do.call(rbind, de_rows)
de <- merge(de, setNames(paper[, c("Symbol", "log2FoldChange")], c("Symbol", "paper_lfc")), by = "Symbol")
de$same_direction <- sign(de$log2FoldChange) == sign(de$paper_lfc)
write.table(de, "results/tables/05_within_cohort_deg210.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, expr_rows), "results/tables/05_within_cohort_aim2_fhod3_expr.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
tests <- do.call(rbind, test_rows)
write.table(tests, "results/tables/05_within_cohort_aim2_fhod3_tests_R.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(tests, digits = 3)
# 210개 요약과 그림은 src/py/05_figures.py 에서 낸다.
