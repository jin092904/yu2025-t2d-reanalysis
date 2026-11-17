"""R 의 sva::ComBat_seq·DESeq2 를 Python 판(inmoose·pydeseq2)으로 따로 구현한다.

사전 명세는 src/r/02_combat_deseq2.R 와 같다. 탐색용 인자도 같다: group(ComBat-Seq 에 군 공변량), pc(protein-coding 만).
pydeseq2 의 alpha 기본값은 0.05 라서 DESeq2 results() 기본값(0.1)에 맞춘다(독립 필터링 문턱에 쓰인다).
실행: python src/py/02_combat_deseq2.py [group] [pc]  (프로젝트 폴더에서)
"""
import sys
from pathlib import Path

import pandas as pd
from inmoose.pycombat import pycombat_seq
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parents[2]
use_group, use_pc = "group" in sys.argv, "pc" in sys.argv
tag = "py" + ("_group" if use_group else "") + ("_pc" if use_pc else "")

samples = pd.read_csv(ROOT / "data/processed/samples.tsv", sep="\t").set_index("gsm")
counts = pd.read_csv(ROOT / "data/processed/counts_main.tsv.gz", sep="\t", index_col="GeneID")
meta = samples.loc[counts.columns]
annot = pd.read_csv(ROOT / "data/raw/Human.GRCh38.p13.annot.tsv.gz", sep="\t", usecols=["GeneID", "Symbol", "GeneType"])
if use_pc:
    counts = counts.loc[counts.index.isin(annot.GeneID[annot.GeneType == "protein-coding"])]

adj = pycombat_seq(counts, meta.series.tolist(), covar_mod=meta.paper_group.tolist() if use_group else None)
adj = adj.round().astype(int)  # R ComBat_seq 도 정수로 돌려준다

dds = DeseqDataSet(counts=adj.T, metadata=meta.rename(columns={"paper_group": "condition"})[["condition"]],
                   design="~condition", quiet=True, n_cpus=1)
dds.deseq2()
st = DeseqStats(dds, contrast=["condition", "OT2D", "OH"], alpha=0.1, quiet=True, n_cpus=1)
st.summary()
res = st.results_df.rename_axis("GeneID").reset_index()
res["GeneID"] = res.GeneID.astype(int)
res = res.merge(annot[["GeneID", "Symbol"]], on="GeneID", how="left")
res["deg"] = (res.padj < 0.05) & (res.log2FoldChange.abs() >= 2)

out = ROOT / "results/tables"
out.mkdir(parents=True, exist_ok=True)
res.to_csv(out / f"02_deseq2_{tag}.tsv.gz", sep="\t", index=False)
adj.to_csv(ROOT / f"data/processed/combat_counts_{tag}.tsv.gz", sep="\t")
up, down = (res.deg & (res.log2FoldChange > 0)).sum(), (res.deg & (res.log2FoldChange < 0)).sum()
print(f"[{tag}] 시료 {adj.shape[1]}, 유전자 {adj.shape[0]}, DEG {res.deg.sum()} (상향 {up}, 하향 {down})")
