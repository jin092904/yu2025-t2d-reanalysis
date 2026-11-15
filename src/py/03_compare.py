"""Stage 3 대조: 표 2(논문 유전자 목록) vs R 사전 명세 vs Python, 그리고 시드 100개 반복(분할 운).

입력: results/tables/03_ml_{R,py}_vst_table1_sample_no.tsv, 03_seed_stability_vst_table1_sample_no.tsv
출력: results/tables/03_table2_comparison.tsv, 03_seed_hub_frequency.tsv, results/figures/03_hub_frequency.{png,pdf}
실행: python src/py/03_compare.py  (프로젝트 폴더에서)
"""
from collections import Counter
from pathlib import Path

import matplotlib
import pandas as pd
from scipy.stats import rankdata

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Noto Sans CJK JP"  # 한글
ROOT = Path(__file__).resolve().parents[2]
TAB, FIG = ROOT / "results/tables", ROOT / "results/figures"
TAG = "vst_table1_sample_no"

# 논문 표 2 (paper/PMC12748723.xml 표 2 에서 옮김)
PAPER = {
    "LASSO": "CD52 PROK1 AIM2 CHI3L1 UGT1A1 PNPLA1 STMN2 CLEC12A PPP4R4 CCL18 FHOD3".split(),
    "SVM-RFE": ("AIM2 AKR1B15 ALAS2 ALPK3 AQP9 BCL2A1 CA1 CD52 CHI3L1 CR1L FADS1 FASN FHOD3 GLT1D1 GPAT3 "
                "IFIT1B LRRC3C MCEMP1 MEIKIN NFE2 PHLDA2 PLEK2 PNPLA1").split(),
    "RF": "ALPK3 MCEMP1 SLC4A1 ALAS2 PHLDA2 FADS1 FHOD3 AIM2 RUNDC3A GPAT3".split(),
    "hub": ["AIM2", "FHOD3"],
}
assert [len(PAPER[k]) for k in ("LASSO", "SVM-RFE", "RF")] == [11, 23, 10]


def genes_by_method(path):
    t = pd.read_csv(path, sep="\t")
    return {m: t.gene[t.method == m].tolist() for m in PAPER}


R, PY = genes_by_method(TAB / f"03_ml_R_{TAG}.tsv"), genes_by_method(TAB / f"03_ml_py_{TAG}.tsv")
rows = []
for m in PAPER:
    p, r, q = set(PAPER[m]), set(R[m]), set(PY[m])
    rows.append(dict(method=m, paper_n=len(p), R_n=len(r), R_in_paper=len(r & p), py_n=len(q), py_in_paper=len(q & p),
                     R_py_shared=len(r & q), R_genes=", ".join(R[m]), py_genes=", ".join(PY[m])))
cmp = pd.DataFrame(rows)
cmp.to_csv(TAB / "03_table2_comparison.tsv", sep="\t", index=False)
print(cmp.drop(columns=["R_genes", "py_genes"]).to_string(index=False), "\n")

# 사전 명세 학습 20명에서 유전자 하나로 두 군이 완전히 갈리는(max(AUC, 1−AUC) = 1) 유전자
x = pd.read_csv(ROOT / "data/processed/ml_input_vst.tsv.gz", sep="\t").set_index("gsm")
x = x.loc[pd.read_csv(ROOT / f"data/processed/03_split_{TAG}.tsv", sep="\t").gsm]
yy = (x.pop("condition") == "OT2D").values
r = x.apply(rankdata)
auc = (r[yy].sum() - yy.sum() * (yy.sum() + 1) / 2) / (yy.sum() * (~yy).sum())
perfect = sorted(auc.index[(auc == 1) | (auc == 0)])
pd.Series(perfect, name="gene").to_csv(TAB / "03_train_perfect_auc_genes.tsv", sep="\t", index=False)
print(f"학습 20명에서 AUC = 1 인 유전자 {len(perfect)}개: {', '.join(perfect)}\n")

# 시드 100개: 방법별로 각 유전자가 뽑힌 비율
st = pd.read_csv(TAB / f"03_seed_stability_{TAG}.tsv", sep="\t", keep_default_na=False)
nseed = st.seed.nunique()
freq = {m: Counter(g for s in st.genes[st.method == m] for g in s.split(",") if g) for m in PAPER}
genes = sorted(set().union(*freq.values()))
tab = pd.DataFrame({m: [freq[m][g] / nseed for g in genes] for m in PAPER}, index=genes).sort_values("hub", ascending=False)
tab.to_csv(TAB / "03_seed_hub_frequency.tsv", sep="\t")
hub_sets = st.genes[st.method == "hub"]
print(f"시드 {nseed}개 — 세 방법 공통 유전자 수 분포:", Counter(st.n[st.method == "hub"]).most_common())
print("공통이 정확히 {AIM2, FHOD3}:", (hub_sets.map(lambda s: set(s.split(",")) == {"AIM2", "FHOD3"})).sum())
print("AIM2·FHOD3 둘 다 공통에 포함:", (hub_sets.map(lambda s: {"AIM2", "FHOD3"} <= set(s.split(",")))).sum())
print("방법별 선택 비율(AIM2, FHOD3):\n", tab.loc[["AIM2", "FHOD3"]].round(2))
print("공통 비율 상위 10:\n", tab.head(10).round(2))

# 그림: 세 방법이 모두 고른 비율 상위 10개 + 논문 허브 유전자(순위 밖이어도 넣음, 주황)
sel = list(tab.head(10).index) + [g for g in PAPER["hub"] if g not in tab.head(10).index]
top = tab.loc[sel].sort_values("hub").copy()
fig, ax = plt.subplots(figsize=(6.5, 5))
colors = ["#B45309" if g in PAPER["hub"] else "#94a3b8" for g in top.index]
ax.barh(top.index, top.hub * 100, color=colors)
for i, v in enumerate(top.hub * 100):
    ax.text(v + 1, i, f"{v:.0f}%", va="center", fontsize=10)
ax.set_xlabel(f"세 방법이 모두 고른 비율 (분할 시드 {nseed}개)")
ax.set_xlim(0, 105)
a, f = (round(tab.hub[g] * nseed) for g in ("AIM2", "FHOD3"))
ax.set_title(f"분할 {nseed}번 중 논문 허브 유전자(주황)는 AIM2 {a}번 · FHOD3 {f}번만 다시 나온다", fontsize=12)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"03_hub_frequency.{ext}", dpi=300)
