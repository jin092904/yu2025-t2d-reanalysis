"""Stage 4 (Python): src/r/04_roc.R 와 같은 AUC·DeLong 95% CI 를 scikit-learn + 직접 구현으로 계산하고 R 과 대조한다.

분할은 R 이 낸 data/processed/04_splits_100.tsv(시드 1~100)를 그대로 쓴다(Python 은 R 의 난수열을 만들 수 없다). 사전 명세 분할 = 시드 100.
DeLong 분산은 위치값(placement) 공식으로 직접 계산한다(pROC 의 method = "delong" 과 같은 식).
출력: results/tables/04_roc_py_main.tsv, 04_test_auc_100seeds_py.tsv, 04_roc_R_vs_py.tsv, results/figures/04_test_auc_100seeds.{png,pdf}
실행: python src/py/04_roc.py  (프로젝트 폴더에서)
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Noto Sans CJK JP"  # 한글
ROOT = Path(__file__).resolve().parents[2]
TAB, FIG = ROOT / "results/tables", ROOT / "results/figures"

d = pd.read_csv(ROOT / "data/processed/ml_input_vst.tsv.gz", sep="\t").set_index("gsm")
y = (d.pop("condition") == "OT2D").astype(int)
lfc = pd.read_csv(TAB / "02_paper_vs_R_group_pc.tsv", sep="\t").set_index("Symbol").log2FoldChange
sign = np.where(lfc.loc[d.columns] > 0, 1, -1)  # DEG 방향: 하향 유전자는 값을 뒤집어 AUC 를 잰다
splits = pd.read_csv(ROOT / "data/processed/04_splits_100.tsv", sep="\t")


def delong_ci(score, yy):
    """AUC 와 DeLong 95% CI(위치값 공식). AUC = 1 이면 분산 0 → CI 1–1 (pROC 와 같음)."""
    pos, neg = score[yy == 1], score[yy == 0]
    psi = (pos[:, None] > neg[None, :]) + 0.5 * (pos[:, None] == neg[None, :])
    auc = psi.mean()
    v10, v01 = psi.mean(1), psi.mean(0)
    se = np.sqrt(v10.var(ddof=1) / len(pos) + v01.var(ddof=1) / len(neg))
    return auc, max(auc - norm.ppf(0.975) * se, 0), min(auc + norm.ppf(0.975) * se, 1)


# 1) 사전 명세 분할(시드 100)
s100 = splits[splits.seed == 100].set_index("gsm").loc[d.index, "set"]
tr, te = (s100 == "train").values, (s100 == "test").values
X2 = d[["AIM2", "FHOD3"]].values
lr = LogisticRegression(penalty=None, max_iter=100000).fit(X2[tr], y.values[tr])  # 완전 분리면 계수가 커질 뿐 순위는 정해진다
prob = lr.predict_proba(X2)[:, 1]
rows = []
for set_name, m in (("train", tr), ("test", te)):
    yy = y.values[m]
    for model, score in (("AIM2", d.AIM2.values), ("FHOD3", -d.FHOD3.values), ("AIM2+FHOD3", prob)):
        auc, lo, hi = delong_ci(score[m], yy)
        rows.append(dict(set=set_name, model=model, auc=auc, lo=lo, hi=hi))
main = pd.DataFrame(rows)
main.to_csv(TAB / "04_roc_py_main.tsv", sep="\t", index=False)
r_main = pd.read_csv(TAB / "04_roc_R_main.tsv", sep="\t")
cmp = r_main.merge(main, on=["set", "model"], suffixes=("_R", "_py"))
cmp.to_csv(TAB / "04_roc_R_vs_py.tsv", sep="\t", index=False)
print(cmp[["set", "model", "n_OH", "n_OT2D", "auc_R", "lo_R", "auc_py", "lo_py", "paper_auc"]].round(4).to_string(index=False))

# 2) 시드 1~100: 검증 AUC (허브 유전자, DEG 210개)
rows = []
for s, g in splits.groupby("seed"):
    m = (g.set_index("gsm").loc[d.index, "set"] == "test").values
    auc_all = np.array([roc_auc_score(y.values[m], sign[j] * d.values[m, j]) for j in range(d.shape[1])])
    a = dict(zip(d.columns, auc_all))
    rows.append(dict(seed=s, AIM2=a["AIM2"], FHOD3=a["FHOD3"], deg_median=np.median(auc_all),
                     deg_ge09=(auc_all >= 0.9).mean(), deg_eq1=(auc_all == 1).mean()))
seeds = pd.DataFrame(rows)
seeds.to_csv(TAB / "04_test_auc_100seeds_py.tsv", sep="\t", index=False)
r_seeds = pd.read_csv(TAB / "04_test_auc_100seeds_R.tsv", sep="\t")
diff = (r_seeds[["AIM2", "FHOD3", "deg_median", "deg_ge09"]] - seeds[["AIM2", "FHOD3", "deg_median", "deg_ge09"]]).abs().max()
print("\n시드 100개 R vs Python 최대 차:", diff.round(12).to_dict())
print("검증 AUC = 1 인 시드: AIM2", (seeds.AIM2 == 1).sum(), "/ FHOD3", (seeds.FHOD3 == 1).sum())
print("검증 AUC 최솟값: AIM2 %.3f, FHOD3 %.3f" % (seeds.AIM2.min(), seeds.FHOD3.min()))
print("DEG 210개 검증 AUC 중앙값(시드별)의 범위: %.3f–%.3f, 중앙 %.3f" % (seeds.deg_median.min(), seeds.deg_median.max(), seeds.deg_median.median()))
print("검증 AUC = 1 인 DEG 비율 중앙: %.3f, ≥ 0.9 비율 중앙: %.3f" % (seeds.deg_eq1.median(), seeds.deg_ge09.median()))
comp = r_seeds.groupby(["n_OH_test", "n_OT2D_test"]).size()
print("검증 11명 구성(OH, OT2D) 분포:", comp.to_dict())

# 3) 그림: 시드 100개 검증 AUC — AIM2, FHOD3, DEG 210개 중앙값
fig, ax = plt.subplots(figsize=(6.5, 4.6))
cols = [("AIM2", "AIM2", "#B45309"), ("FHOD3", "FHOD3", "#B45309"), ("deg_median", "DEG 210개\n중앙값", "#64748b")]
rng = np.random.default_rng(0)
for i, (c, label, color) in enumerate(cols):
    v = seeds[c].values
    ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=14, alpha=0.6, color=color)
    ax.hlines(np.median(v), i - 0.3, i + 0.3, color="black", lw=2)
    ax.text(i + 0.33, np.median(v), f"{np.median(v):.2f}", va="center", fontsize=10)
ax.set_xticks(range(3), [c[1] for c in cols])
ax.set_ylabel("검증 11명 AUC (분할 시드 100개)")
ax.set_ylim(0.7, 1.02)
ax.set_title(f"검증 AUC: 허브 유전자는 대부분 1, DEG 210개 전체도 중앙값 {seeds.deg_median.median():.2f}", fontsize=12)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"04_test_auc_100seeds.{ext}", dpi=300)
