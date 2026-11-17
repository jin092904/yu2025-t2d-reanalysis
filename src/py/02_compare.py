"""논문 표 S1·S2(DEG 210개)와 R·Python 실행 결과를 유전자 단위로 비교한다.

입력: results/tables/02_deseq2_{R,R_pc,R_group,R_group_pc,py,py_group_pc}.tsv.gz (02_combat_deseq2.R/.py 출력)
출력: results/tables/02_deg_by_spec.tsv, 02_paper_vs_R_group_pc.tsv, results/figures/02_deg_by_spec.{png,pdf}
실행: python src/py/02_compare.py  (프로젝트 폴더에서)
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Noto Sans CJK JP"  # 한글

ROOT = Path(__file__).resolve().parents[2]
TAB, FIG = ROOT / "results/tables", ROOT / "results/figures"
FIG.mkdir(parents=True, exist_ok=True)

xl = ROOT / "paper/41598_2025_29141_MOESM2_ESM.xlsx"
paper = pd.concat([pd.read_excel(xl, sheet_name=s, header=2, index_col=0) for s in ("Table S1", "Table S2")])
paper.index.name = "Symbol"
assert len(paper) == 210 and (paper.log2FoldChange > 0).sum() == 148

# 1) 설정별 DEG 수와 논문 210개 중 몇 개를 찾았나
RUNS = [  # (파일 꼬리, 언어, ComBat-Seq 군 공변량, protein-coding 만, 구분)
    ("R", "R", "없음", "아니오", "사전 명세"),
    ("R_pc", "R", "없음", "예", "탐색"),
    ("R_group", "R", "있음", "아니오", "탐색"),
    ("R_group_pc", "R", "있음", "예", "탐색(재현)"),
    ("py", "Python", "없음", "아니오", "사전 명세"),
    ("py_group_pc", "Python", "있음", "예", "탐색(재현)"),
]
res, rows = {}, []
for tag, lang, grp, pc, kind in RUNS:
    r = pd.read_csv(TAB / f"02_deseq2_{tag}.tsv.gz", sep="\t")
    res[tag] = r
    deg = set(r.Symbol[r.deg.astype(bool)])
    rows.append(dict(run=tag, language=lang, combat_group=grp, protein_coding_only=pc, kind=kind,
                     genes=len(r), deg=len(deg), up=int((r.deg.astype(bool) & (r.log2FoldChange > 0)).sum()),
                     paper_recovered=len(deg & set(paper.index)), extra=len(deg - set(paper.index))))
by_spec = pd.DataFrame(rows)
by_spec.to_csv(TAB / "02_deg_by_spec.tsv", sep="\t", index=False)
print(by_spec.to_string(index=False), "\n")

# 2) 재현 설정(R_group_pc): 논문 210개 유전자별 수치 차
r = res["R_group_pc"].set_index("Symbol")
m = paper.join(r, rsuffix="_R")
cols = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
rel = pd.DataFrame({c: (m[c + "_R"] - m[c]).abs() / m[c].abs() for c in cols})
print("논문 vs R(재현 설정) 상대 차 최대:", rel.max().map("{:.1e}".format).to_dict())
ratio = (m.padj_R / m.padj)
print(f"padj 비(R/논문): 중앙 {ratio.median():.5f}, 범위 {ratio.min():.5f}–{ratio.max():.5f}")
tested = r.padj.notna().sum()
print(f"BH 대상 유전자: R {tested:,} → 논문 추정 {tested / ratio.median():,.0f}")
print("논문에만 있는 DEG:", m.index[~m.deg.astype(bool)].tolist(),
      m.loc[~m.deg.astype(bool), ["padj", "padj_R"]].round(5).values.tolist())
m[cols + [c + "_R" for c in cols]].to_csv(TAB / "02_paper_vs_R_group_pc.tsv", sep="\t")

# 3) 같은 설정의 R vs Python
rc = pd.read_csv(ROOT / "data/processed/combat_counts_R_group_pc.tsv.gz", sep="\t", index_col=0)
pc = pd.read_csv(ROOT / "data/processed/combat_counts_py_group_pc.tsv.gz", sep="\t", index_col=0).loc[rc.index, rc.columns]
same = (rc.values == pc.values)
print(f"\nComBat-Seq 보정 count R=Python: 칸 {same.mean():.4%}, 유전자(31칸 모두) {same.all(1).mean():.2%}")
a, b = res["R_group_pc"].set_index("GeneID"), res["py_group_pc"].set_index("GeneID").loc[res["R_group_pc"].GeneID]
both = a.deg.astype(bool) & b.deg.astype(bool)
print(f"DEG R {a.deg.sum()} / Python {b.deg.sum()} / 공통 {both.sum()}")
lfc = (a.log2FoldChange - b.log2FoldChange).abs()
print(f"log2FC 차(R−Python) 중앙 {lfc.median():.1e}, 99% {lfc.quantile(.99):.1e}")
only_r = a[a.deg.astype(bool) & ~b.deg.astype(bool)]
print("R 에만 있는 DEG 의 baseMean:", only_r.baseMean.round(2).tolist())

# 4) 그림: 설정별 DEG 수와 논문 일치 수
lab = {"R": "R 사전 명세\n(군 공변량 없음)", "R_pc": "R + protein-coding", "R_group": "R + 군 공변량",
       "R_group_pc": "R + 군 공변량\n+ protein-coding", "py": "Python 사전 명세", "py_group_pc": "Python + 군 공변량\n+ protein-coding"}
d = by_spec.iloc[::-1]
fig, ax = plt.subplots(figsize=(7, 4.6))
y = np.arange(len(d))
ax.barh(y, d.deg, color="#cbd5e1", label="찾은 DEG")
ax.barh(y, d.paper_recovered, color="#005AAB", label="그중 논문 210개와 같은 유전자")
for yi, (n, k) in enumerate(zip(d.deg, d.paper_recovered)):
    ax.text(n + 8, yi, f"{n} ({k})", va="center", fontsize=10)
ax.axvline(210, ls="--", color="#B45309", lw=1.2)
ax.text(212, len(d) - 0.45, "논문 210", color="#B45309", fontsize=10)
ax.set_yticks(y, [lab[t] for t in d.run], fontsize=10)
ax.set_xlabel("DEG 수 (padj < 0.05, |log2FC| ≥ 2)")
ax.set_xlim(0, 720)
ax.set_title("ComBat-Seq 에 군 공변량을 넣어야 논문 DEG 가 나온다", fontsize=12)
ax.legend(loc="lower right", fontsize=9, frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"02_deg_by_spec.{ext}", dpi=300)
print("\n그림:", FIG / "02_deg_by_spec.png")
