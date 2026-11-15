# Stage 7: Matters Arising 원고용 영문 그림(ggplot2). 저장된 결과 파일만 읽는다(새 수치 없음).
# 양식은 원 논문 그림 6·7(ggpubr 상자그림, 정상 파랑 · 당뇨 빨강, 비교 괄호와 P 값)에 맞춘다.
#   Figure 1 코호트 안 비교: (a, b) GSE141432 FHOD3 · AIM2, (c, d) GSE166047 FHOD3 · AIM2
#   Figure 2 재현 · 분할 · 누수: (a) ComBat-seq 설정별 DEG 수, (b) 분할 방식 두 가지에서 세 방법 공통 비율, (c) 누수 비교
#   Supplementary Figure 1: 논문 DEG 가 코호트 안 비교마다 논문과 같은 방향인 비율
# 출력: manuscript/figures/Figure1_v7 · Figure2_v7 · SupplementaryFigure1_v1 .{pdf,png} (cairo 벡터 PDF + 600 dpi PNG)
# 실행: Rscript src/r/07_manuscript_figures.R  (프로젝트 폴더에서)
suppressPackageStartupMessages({library(ggplot2); library(ggpubr); library(patchwork)})

rd <- function(name, ...) read.delim(file.path("results/tables", name), check.names = FALSE, ...)
OUT <- "manuscript/figures"
FONT <- "Liberation Sans"                         # Arial 과 글자 폭이 같은 글꼴
BLUE <- "#3C5BD9"; RED <- "#E02B2B"; GREY <- "#9E9E9E"; LIGHT <- "#D9D9D9"; TEAL <- "#2A9D8F"   # 원 논문 그림 6·7 의 파랑·빨강
theme_ms <- theme_classic(base_size = 8, base_family = FONT) +
  theme(plot.tag = element_text(face = "bold", size = 12), axis.text = element_text(colour = "black"),
        plot.title = element_text(size = 8, hjust = 0.5), legend.key.size = unit(3, "mm"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "in", dpi = 600, type = "cairo")
  cat(file.path(OUT, paste0(name, ".pdf")), "\n")
}
sup <- function(x) chartr("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079", x)
pfmt <- function(p) {  # P < 0.001 은 "8.2 × 10⁻⁵" 꼴
  e <- floor(log10(p)); m <- p / 10^e
  ifelse(p < 0.001, sprintf("P = %.1f \u00d7 10%s", m, sup(as.character(e))),
         ifelse(p > 0.995, "P = 1.0", sprintf("P = %.2g", p)))   # 본문 표기 "P = 1.0" 과 맞춤
}

## Figure 2: 재현 · 분할 · 누수
# (a) R 네 설정. 논문 목록에 든 DEG(빨강)와 나머지(회색)를 쌓은 막대. Python 결과는 Methods 에 적는다.
spec <- rd("02_deg_by_spec.tsv")
COND <- c(R = "Batch only\n(as written)", R_pc = "Protein-coding\nonly", R_group = "Group\ncovariate", R_group_pc = "Group covariate\n+ protein-coding")
sr <- spec[spec$run %in% names(COND), ]
sr$cond <- factor(COND[sr$run], levels = rev(COND))
sa <- rbind(data.frame(cond = sr$cond, part = "In published list", n = sr$paper_recovered),
            data.frame(cond = sr$cond, part = "Other DEGs", n = sr$deg - sr$paper_recovered))
sa$part <- factor(sa$part, levels = c("Other DEGs", "In published list"))   # 빨강이 0 쪽(왼쪽)
pa <- ggplot(sa, aes(x = n, y = cond, fill = part)) +
  geom_col(width = 0.62) +
  geom_text(data = sr, aes(x = deg, y = cond, label = sprintf("%d (%d)", deg, paper_recovered)), inherit.aes = FALSE,
            hjust = -0.1, size = 2.5, family = FONT) +
  geom_vline(xintercept = 210, linetype = "dashed", colour = RED, linewidth = 0.4) +
  annotate("text", x = 218, y = 4.45, label = "published 210", hjust = 0, size = 2.4, colour = RED, family = FONT) +
  scale_fill_manual(values = c(`In published list` = RED, `Other DEGs` = LIGHT), breaks = c("In published list", "Other DEGs"), name = NULL) +
  scale_x_continuous(limits = c(0, 830), breaks = seq(0, 800, 200), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = expression("DEGs (padj < 0.05, |log"[2] * "FC| ">= 2 * ")"), y = NULL) + theme_ms +
  theme(legend.position = c(0.72, 0.84), legend.background = element_blank(), legend.text = element_text(size = 6.5),
        legend.key.size = unit(2.6, "mm"))

# (b) 분할 방식 두 가지: Methods 의 단순 무작위(사전 정의)와 Results 의 군마다 65%(07_pergroup_split.R)
hubs_of <- function(f) {
  st <- rd(f, colClasses = "character"); st <- st[st$method == "hub", ]
  lapply(split(st$genes, st$seed), function(g) if (g == "") character(0) else strsplit(g, ",")[[1]])
}
freq_of <- function(h) { tb <- table(unlist(h)) / length(h); setNames(as.numeric(tb), names(tb)) }
hm <- hubs_of("03_seed_stability_vst_table1_sample_no.tsv"); hr <- hubs_of("03_seed_stability_vst_table1_pergroup_no.tsv")
stopifnot(length(hm) == 100, length(hr) == 100)
fm <- freq_of(hm); fr <- freq_of(hr)
both <- sapply(list(hm, hr), function(H) sum(sapply(H, function(h) all(c("AIM2", "FHOD3") %in% h))))
exact <- sapply(list(hm, hr), function(H) sum(sapply(H, function(h) setequal(h, c("AIM2", "FHOD3")))))
genes <- unique(c(names(sort(fm, decreasing = TRUE))[1:8], names(sort(fr, decreasing = TRUE))[1:8], "AIM2", "FHOD3"))
val <- function(f, g) ifelse(is.na(f[g]), 0, f[g])
f <- rbind(data.frame(gene = genes, split = "Simple random (Methods)", pct = 100 * val(fm, genes)),
           data.frame(gene = genes, split = "Stratified (Results)", pct = 100 * val(fr, genes)))
ord <- genes[order(-val(fm, genes), -val(fr, genes))]              # 단순 무작위 기준 내림차순, 같으면 군별 분할 기준
f$gene <- factor(f$gene, levels = rev(ord))
f$split <- factor(f$split, levels = c("Stratified (Results)", "Simple random (Methods)"))   # 묶음 안에서 단순 무작위가 위
hubcol <- ifelse(rev(ord) %in% c("AIM2", "FHOD3"), RED, "black")   # 논문 허브 유전자 = 빨간 이름(막대는 회색 두 톤)
pb <- ggplot(f, aes(x = pct, y = gene, fill = split)) +
  geom_col(position = position_dodge(0.8), width = 0.76) +
  geom_text(aes(label = round(pct)), position = position_dodge(0.8), hjust = -0.25, size = 2.0, family = FONT, colour = "grey25") +
  annotate("text", x = 55.5, y = 2.5, hjust = 1, vjust = 0.5, size = 2.2, family = FONT, colour = RED, lineheight = 0.95,
           label = sprintf("AIM2 & FHOD3 together:\nsimple random %d/100\nstratified %d/100", both[1], both[2])) +
  scale_fill_manual(values = c(`Simple random (Methods)` = "grey40", `Stratified (Results)` = "grey78"),
                    breaks = c("Simple random (Methods)", "Stratified (Results)"), name = NULL) +
  scale_x_continuous(limits = c(0, 56), breaks = seq(0, 50, 10), expand = c(0, 0)) +
  labs(x = "Selected by all three methods\n(% of 100 random splits)", y = NULL) + theme_ms +
  theme(axis.text.y = suppressWarnings(element_text(colour = hubcol, face = ifelse(hubcol == RED, "bold", "plain"))),
        legend.position = "top", legend.direction = "vertical", legend.justification = "left",
        legend.margin = margin(0, 0, -4, 0), legend.text = element_text(size = 6.5), legend.key.size = unit(2.6, "mm"))

# (c) 누수: 같은 분할끼리 선으로 잇고, 두 중앙값을 굵은 선으로 잇는다
lk <- rd("05_leakage_100seeds.tsv")
lk <- lk[!is.na(lk$auc_clean), ]
long <- rbind(data.frame(seed = lk$seed, sel = "All 31\n(published)", auc = lk$auc_leaky),   # 왼쪽 = 논문 방식, 오른쪽 = 학습만(본문 "0.733 → 0.664" 순서)
              data.frame(seed = lk$seed, sel = "Training\nonly", auc = lk$auc_clean))
long$sel <- factor(long$sel, levels = c("All 31\n(published)", "Training\nonly"))
med <- aggregate(auc ~ sel, long, median)
n_lower <- sum(lk$auc_clean < lk$auc_leaky)
pc <- ggplot(long, aes(x = sel, y = auc)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = GREY, linewidth = 0.4) +
  annotate("text", x = 2.55, y = 0.5, label = "0.5", hjust = 0, vjust = -0.3, size = 2.2, colour = "grey30", family = FONT) +
  geom_line(aes(group = seed), colour = LIGHT, linewidth = 0.2, alpha = 0.7) +
  geom_point(aes(colour = sel), size = 0.7, alpha = 0.55) +
  geom_line(data = med, aes(group = 1), colour = "black", linewidth = 0.8) +
  geom_crossbar(data = med, aes(ymin = auc, ymax = auc), width = 0.4, linewidth = 0.45, colour = "black") +
  geom_text(data = med, aes(label = sprintf("%.3f", auc)), nudge_x = c(-0.36, 0.36), hjust = c(1, 0), size = 2.6, family = FONT) +
  annotate("text", x = 1.5, y = 0.875, label = sprintf("Lower in %d of %d splits", n_lower, nrow(lk)), size = 2.4, colour = "grey25", family = FONT) +
  scale_x_discrete(expand = expansion(add = 0.95)) +
  scale_y_continuous(limits = c(0.45, 0.885), breaks = seq(0.45, 0.85, 0.1)) +
  scale_colour_manual(values = c(RED, "grey35"), guide = "none") +
  labs(x = sprintf("DEG selection (%d splits)", nrow(lk)), y = "Test-set AUC\n(median over DEGs)") + theme_ms

fig_repro <- (pa | pb | pc) + plot_layout(widths = c(1.15, 1.1, 0.85)) + plot_annotation(tag_levels = "a")
save_fig(fig_repro, "Figure2_v7", 7.1, 3.0)

## Figure 1: 코호트 안 비교, Supplementary Figure 1: 논문 DEG 방향 비율
ex <- rd("05_within_cohort_aim2_fhod3_expr.tsv")
tests <- rd("05_within_cohort_aim2_fhod3_tests_R.tsv")
GROUPS <- c(lean = "Lean", obese = "Obese\nnon-T2D", t2d = "Obese\nT2D")
ex$grp <- factor(GROUPS[ex$group], levels = GROUPS)
COH <- c(GSE141432 = "GSE141432 (men and women)", GSE166047 = "GSE166047 (women)")
panel <- function(gene, cohort) {   # y축과 괄호 높이는 유전자별로 두 코호트가 같다(위·아래 줄을 눈으로 비교할 수 있게)
  g <- ex[ex$gene == gene, ]
  d <- g[g$cohort == cohort, ]
  t <- tests[tests$gene == gene & tests$cohort == cohort, ]   # 비만 · 당뇨 · 논문형(정상체중 대 비만 당뇨) 세 비교
  top <- max(g$log2_norm); rng <- diff(range(g$log2_norm))
  br <- data.frame(group1 = GROUPS[c("lean", "obese", "lean")], group2 = GROUPS[c("obese", "t2d", "t2d")],
                   p = pfmt(t$p_wilcox[match(c("obesity", "diabetes", "paper_like"), t$comparison)]),
                   y.position = top + rng * c(0.08, 0.20, 0.32))
  ggplot(d, aes(x = grp, y = log2_norm)) +
    geom_boxplot(aes(fill = grp), width = 0.55, outlier.shape = NA, linewidth = 0.3, alpha = 0.9) +
    geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 0.7) +   # 시드 고정: PDF · PNG · 다시 돌려도 점 위치가 같다
    stat_pvalue_manual(br, label = "p", size = 2.3, tip.length = 0.01, family = FONT, bracket.size = 0.3) +
    scale_fill_manual(values = c(BLUE, GREY, RED), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    coord_cartesian(ylim = c(min(g$log2_norm), top + rng * 0.36)) +
    labs(x = NULL, y = bquote(.(gene) * ", log"[2] * "(normalized count + 1)"), title = COH[[cohort]]) + theme_ms
}
pe_dat <- rd("05_within_cohort_deg210_summary.tsv")
pe_dat$pct <- 100 * pe_dat$same_direction / pe_dat$genes
pe_dat$comp <- factor(c(obesity = "Obesity", diabetes = "T2D", paper_like = "Published\ncontrast")[pe_dat$comparison],
                      levels = c("Obesity", "T2D", "Published\ncontrast"))
pe <- ggplot(pe_dat, aes(x = comp, y = pct, fill = cohort)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%s", pct, ifelse(cohort == "GSE166047" & comparison == "obesity", "\u2020", ""))),
            position = position_dodge(0.75), vjust = -0.3, size = 2.3, family = FONT) +   # † = 원 분석 시료가 없는 유일한 비교
  geom_hline(yintercept = 50, linetype = "dashed", colour = GREY, linewidth = 0.4) +
  scale_fill_manual(values = c(GSE141432 = "#A9D8D1", GSE166047 = "#F2A33A"), name = NULL,   # a–d 의 파랑(정상체중)과 겹치지 않게
                    labels = c(GSE141432 = "GSE141432 (used for selection)", GSE166047 = "GSE166047")) +   # GSE141432 = 순환 비교라 옅게
  scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = NULL, y = "Published DEGs changing in the published direction (%)",
       caption = "\u2020 Only comparison without samples\nfrom the original analysis") + theme_ms +
  theme(legend.position = c(0.47, 0.97), legend.background = element_blank(),
        plot.caption = element_text(colour = "grey30", hjust = 0.5, size = 6.5, lineheight = 0.95))

fig_cohort <- ((panel("FHOD3", "GSE141432") | panel("AIM2", "GSE141432")) /      # 본문 순서: GSE141432 FHOD3(a) → AIM2(b) → GSE166047(c, d)
               (panel("FHOD3", "GSE166047") | panel("AIM2", "GSE166047"))) + plot_annotation(tag_levels = "a")
save_fig(fig_cohort, "Figure1_v7", 5.6, 4.6)
save_fig(pe, "SupplementaryFigure1_v1", 3.6, 3.8)
