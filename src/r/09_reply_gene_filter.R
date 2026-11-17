# 원저자 답장(2026-11-17) 점검 ①: "Genes with zero counts across all samples were removed before DESeq2 analysis."
# 이 문장을 여러 가지로 읽어 ComBat-Seq(batch = GEO 시리즈, group = 연구군) + DESeq2 를 다시 돌리고,
# 논문 DEG 210개(보충 표 S1–S2)가 몇 개 나오는지와 논문 대비 padj 비(우리 / 논문)를 본다.
#   none          거르기 없음(= 02 의 R_group)
#   allzero_pre   ComBat-Seq 전에 31개 시료 모두 0 인 유전자 제거
#   allzero_post  ComBat-Seq 뒤, DESeq2 전에 모두 0 인 유전자 제거(답장 문장 그대로)
#   anyzero_pre   ComBat-Seq 전에 한 시료라도 0 인 유전자 제거
#   anyzero_post  ComBat-Seq 뒤, 한 시료라도 0 인 유전자 제거
#   pc            NCBI 주석표의 protein-coding 만(= 02 의 R_group_pc, 기준)
#   pc_allzero_post  protein-coding 만 + ComBat-Seq 뒤 모두 0 인 유전자 제거(답장의 거르기를 protein-coding 위에 더함)
# padj 비가 1 이면 검정한 유전자 수(독립 거르기 뒤)가 논문과 같다는 뜻이다.
# 실행: Rscript src/r/09_reply_gene_filter.R  (프로젝트 폴더에서)
suppressPackageStartupMessages({library(sva); library(DESeq2); library(readxl)})

annot <- read.delim(gzfile("data/raw/Human.GRCh38.p13.annot.tsv.gz"), quote = "")[, c("GeneID", "Symbol", "GeneType")]
samples <- read.delim("data/processed/samples.tsv")
counts <- as.matrix(read.delim(gzfile("data/processed/counts_main.tsv.gz"), row.names = 1, check.names = FALSE))
meta <- samples[match(colnames(counts), samples$gsm), ]
meta$condition <- factor(meta$paper_group, levels = c("OH", "OT2D"))
stopifnot(identical(meta$gsm, colnames(counts)))

xl <- "paper/41598_2025_29141_MOESM2_ESM.xlsx"
paper <- do.call(rbind, lapply(c("Table S1", "Table S2"), function(s) {
  d <- as.data.frame(read_excel(xl, sheet = s, skip = 2)); names(d)[1] <- "Symbol"; d[, c("Symbol", "baseMean", "log2FoldChange", "pvalue", "padj")]
}))
stopifnot(nrow(paper) == 210)

run <- function(name) {
  m <- counts
  if (name == "allzero_pre") m <- m[rowSums(m) > 0, ]
  if (name == "anyzero_pre") m <- m[rowSums(m == 0) == 0, ]
  if (startsWith(name, "pc")) m <- m[rownames(m) %in% annot$GeneID[annot$GeneType == "protein-coding"], ]
  adj <- ComBat_seq(m, batch = meta$series, group = meta$condition)
  if (name %in% c("allzero_post", "pc_allzero_post")) adj <- adj[rowSums(adj) > 0, ]
  if (name == "anyzero_post") adj <- adj[rowSums(adj == 0) == 0, ]
  dds <- DESeq(DESeqDataSetFromMatrix(adj, colData = meta, design = ~condition), quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", "OT2D", "OH")))
  res$Symbol <- annot$Symbol[match(as.integer(rownames(res)), annot$GeneID)]
  deg <- !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 2
  hit <- match(paper$Symbol, res$Symbol)
  ratio <- res$padj[hit] / paper$padj
  rel <- function(k) max(abs(res[[k]][hit] / paper[[k]] - 1), na.rm = TRUE)  # 논문 대비 최대 상대차
  data.frame(filter = name, genes = nrow(adj), tested = sum(!is.na(res$padj)), deg = sum(deg),
             paper_recovered = sum(paper$Symbol %in% res$Symbol[deg]), extra = sum(!res$Symbol[deg] %in% paper$Symbol),
             padj_ratio_median = median(ratio, na.rm = TRUE), padj_ratio_min = min(ratio, na.rm = TRUE),
             padj_ratio_max = max(ratio, na.rm = TRUE), basemean_max_rel_diff = rel("baseMean"),
             lfc_max_rel_diff = rel("log2FoldChange"), pvalue_max_rel_diff = rel("pvalue"))
}

out <- do.call(rbind, lapply(c("none", "allzero_pre", "allzero_post", "anyzero_pre", "anyzero_post", "pc", "pc_allzero_post"), function(n) {
  r <- run(n); print(r); r
}))
out$generated <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
write.table(out, "results/tables/09_reply_gene_filter.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
