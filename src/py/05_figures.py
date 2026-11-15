"""Stage 5 정리: ① 코호트 안 비교(Python 으로 Wilcoxon·AUC 재계산해 R 과 대조) ② 누수 비교, 그림 3장.

입력: results/tables/05_within_cohort_*.tsv (src/r/05_within_cohort.R), 05_leakage_100seeds.tsv (src/r/05_leakage.R)
출력: results/tables/05_within_cohort_aim2_fhod3_tests_py.tsv, 05_within_cohort_deg210_summary.tsv,
      results/figures/05a_aim2_fhod3_within_cohort, 05b_deg210_direction, 05c_leakage (.png, .pdf)
실행: python src/py/05_figures.py  (프로젝트 폴더에서)
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Noto Sans CJK JP"  # 한글
ROOT = Path(__file__).resolve().parents[2]
TAB, FIG = ROOT / "results/tables", ROOT / "results/figures"
GROUP_LABEL = {"lean": "정상체중", "obese": "비만·비당뇨", "t2d": "비만·당뇨"}
COHORT_LABEL = {"GSE141432": "GSE141432 (남녀, 같은 배치)", "GSE166047": "GSE166047 (여성, 같은 배치)"}
COMPS = [("obesity", "obese", "lean"), ("diabetes", "t2d", "obese"), ("paper_like", "t2d", "lean")]


def save(fig, name):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300)


# ① AIM2·FHOD3: Python 으로 정확 Wilcoxon·AUC 를 다시 계산해 R 과 대조
ex = pd.read_csv(TAB / "05_within_cohort_aim2_fhod3_expr.tsv", sep="\t")
rows = []
for (cohort, gene), d in ex.groupby(["cohort", "gene"]):
    for comp, g1, g2 in COMPS:
        a, b = d.log2_norm[d.group == g1].values, d.log2_norm[d.group == g2].values
        u = mannwhitneyu(a, b, alternative="two-sided", method="exact")
        rows.append(dict(cohort=cohort, gene=gene, comparison=comp, auc=u.statistic / (len(a) * len(b)), p_wilcox=u.pvalue))
py = pd.DataFrame(rows)
py.to_csv(TAB / "05_within_cohort_aim2_fhod3_tests_py.tsv", sep="\t", index=False)
r = pd.read_csv(TAB / "05_within_cohort_aim2_fhod3_tests_R.tsv", sep="\t")
m = r.merge(py, on=["cohort", "gene", "comparison"], suffixes=("_R", "_py"))
print("① AIM2·FHOD3 R vs Python 최대 차: AUC %.1e, P %.1e" % ((m.auc_R - m.auc_py).abs().max(), (m.p_wilcox_R - m.p_wilcox_py).abs().max()))
print(m[["cohort", "gene", "comparison", "n1", "n2", "median_diff", "auc_R", "p_wilcox_R"]].round(4).to_string(index=False))

# ① DEG 210개: 코호트 안 비교에서 논문과 같은 방향인 비율, padj < 0.05 수
de = pd.read_csv(TAB / "05_within_cohort_deg210.tsv", sep="\t")
de = de[de.log2FoldChange.notna()].copy()  # GSE166047(여성)에서 Y 염색체 유전자 HSFY1 은 전부 0 이라 빠진다 → 209개
de["same_sig"] = (de.padj < 0.05) & de.same_direction
de["lfc_paper_dir"] = de.log2FoldChange * np.sign(de.paper_lfc)
summ = (de.groupby(["cohort", "comparison"])
        .agg(genes=("Symbol", "size"), same_direction=("same_direction", "sum"), same_dir_padj05=("same_sig", "sum"),
             median_lfc_in_paper_direction=("lfc_paper_dir", "median"))
        .reset_index())
summ.to_csv(TAB / "05_within_cohort_deg210_summary.tsv", sep="\t", index=False)
print("\n① 논문 DEG 210개, 코호트 안 비교\n", summ.round(2).to_string(index=False))

# 그림 05a: AIM2·FHOD3, 코호트 × 유전자
fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4), sharey="row")
rng = np.random.default_rng(0)
colors = {"lean": "#64748b", "obese": "#0f766e", "t2d": "#b45309"}
for i, gene in enumerate(["AIM2", "FHOD3"]):
    for j, cohort in enumerate(["GSE141432", "GSE166047"]):
        ax, d = axes[i, j], ex[(ex.gene == gene) & (ex.cohort == cohort)]
        for k, g in enumerate(["lean", "obese", "t2d"]):
            v = d.log2_norm[d.group == g].values
            ax.scatter(k + rng.uniform(-0.12, 0.12, len(v)), v, s=26, color=colors[g], alpha=0.85)
            ax.hlines(np.median(v), k - 0.25, k + 0.25, color="black", lw=1.8)
        p = r.set_index(["cohort", "gene", "comparison"]).p_wilcox
        ax.set_title(f"{gene} · {cohort}\n당뇨 비교 P = {p[(cohort, gene, 'diabetes')]:.3g} · 비만 비교 P = {p[(cohort, gene, 'obesity')]:.3g}",
                     fontsize=9.5)
        ax.set_xticks(range(3), [GROUP_LABEL[g] for g in ["lean", "obese", "t2d"]], fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        if j == 0:
            ax.set_ylabel(f"{gene} log2(정규화 count + 1)", fontsize=9.5)
fig.suptitle("같은 코호트 안: FHOD3 는 비만에서 이미 내려가고, AIM2 의 당뇨 차이는 두 코호트에서 방향이 다르다", fontsize=10.5)
save(fig, "05a_aim2_fhod3_within_cohort")

# 그림 05b: DEG 210개가 논문과 같은 방향인 비율
fig, ax = plt.subplots(figsize=(7, 4.4))
comp_label = {"obesity": "비만\n(비만·비당뇨 vs 정상체중)", "diabetes": "당뇨\n(비만·당뇨 vs 비만·비당뇨)", "paper_like": "논문형\n(비만·당뇨 vs 정상체중)"}
w = 0.38
for j, (cohort, col) in enumerate([("GSE141432", "#005AAB"), ("GSE166047", "#B45309")]):
    s = summ[summ.cohort == cohort].set_index("comparison").loc[["obesity", "diabetes", "paper_like"]]
    pct = s.same_direction / s.genes * 100
    x = np.arange(3) + (j - 0.5) * w
    ax.bar(x, pct, width=w, color=col, label=COHORT_LABEL[cohort])
    for xi, (v, k) in enumerate(zip(pct, s.same_dir_padj05)):
        ax.text(x[xi], v + 1.5, f"{v:.0f}%\n(유의 {k})", ha="center", fontsize=8.5)
ax.axhline(50, ls="--", color="#94a3b8", lw=1)
ax.text(2.62, 51, "우연 50%", fontsize=8.5, color="#64748b", ha="right")
ax.set_xticks(range(3), [comp_label[c] for c in ["obesity", "diabetes", "paper_like"]], fontsize=9)
ax.set_ylabel("논문 DEG 중 같은 방향 (%)")
ax.set_ylim(0, 118)
ax.legend(fontsize=8.5, frameon=False, loc="upper left")
ax.set_title("논문 DEG: GSE166047 에서는 당뇨 비교의 방향이 오히려 반대 (210개, GSE166047 은 209개)", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
save(fig, "05b_deg210_direction")

# ② 누수: 같은 검증 11명, 고를 때 검증을 썼나만 다름
lk = pd.read_csv(TAB / "05_leakage_100seeds.tsv", sep="\t", keep_default_na=True)
ok = lk[lk.auc_clean.notna()].copy()
diff = ok.auc_leaky - ok.auc_clean
print(f"\n② 누수: 실행 시드 {len(ok)}개(건너뜀 {lk.auc_leaky.isna().sum()}), 검증 AUC 중앙값의 중앙 — 누수 {ok.auc_leaky.median():.3f}, "
      f"누수 없음 {ok.auc_clean.median():.3f}, 차의 중앙 {diff.median():.3f}, 누수 쪽이 높은 시드 {(diff > 0).sum()}/{len(ok)}")
print(f"   학습만으로 고른 DEG 수 중앙 {ok.n_deg.median():.0f} (범위 {ok.n_deg.min():.0f}–{ok.n_deg.max():.0f}), 논문 210개와 겹침 중앙 {ok.overlap.median():.0f}, "
      f"AIM2 포함 {int(ok.aim2.sum())}, FHOD3 포함 {int(ok.fhod3.sum())}")
fig, ax = plt.subplots(figsize=(5.6, 4.8))
for _, rw in ok.iterrows():
    ax.plot([0, 1], [rw.auc_clean, rw.auc_leaky], color="#94a3b8", lw=0.7, alpha=0.6)
ax.scatter(np.zeros(len(ok)), ok.auc_clean, s=14, color="#005AAB", zorder=3)
ax.scatter(np.ones(len(ok)), ok.auc_leaky, s=14, color="#B45309", zorder=3)
for x, v in ((0, ok.auc_clean.median()), (1, ok.auc_leaky.median())):
    ax.hlines(v, x - 0.18, x + 0.18, color="black", lw=2)
    ax.text(x + 0.22, v, f"{v:.2f}", va="center", fontsize=10)
ax.set_xticks([0, 1], ["학습 20명으로만\n유전자 선택", "31명 전체로 선택\n(논문 DEG 210개)"])
ax.set_xlim(-0.4, 1.5)
ax.set_ylabel("검증 11명 AUC (유전자별 중앙값)")
ax.set_title(f"검증 자료를 고를 때 쓰면 AUC 가 {diff.median():.2f} 오른다 (시드 {len(ok)}개)", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)
save(fig, "05c_leakage")
