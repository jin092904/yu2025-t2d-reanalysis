"""GEO 시료 정보로 표본표를 만들고, 논문 군(OH·OT2D)의 count 를 하나로 합친다.

출력: data/processed/samples.tsv (세 GEO 시리즈 전 시료), data/processed/counts_main.tsv.gz (GeneID × 분석 시료)
실행: python src/py/01_data.py  (프로젝트 폴더에서)
"""
import gzip
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW, OUT = ROOT / "data/raw", ROOT / "data/processed"
SERIES = ["GSE175495", "GSE141432", "GSE166047"]


def series_matrix(g):
    """series matrix 의 !Sample_ 줄을 시료 × 항목 표로 읽는다."""
    rows = {}
    for line in gzip.open(RAW / f"{g}_series_matrix.txt.gz", "rt", errors="ignore"):
        if line.startswith("!Sample_"):
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip('"') for v in vals]
            rows.setdefault(key, []).append(vals)
    d = pd.DataFrame({"gsm": rows["!Sample_geo_accession"][0], "title": rows["!Sample_title"][0]})
    for vals in rows["!Sample_characteristics_ch1"]:  # "key: value" 형식
        key = vals[0].split(":")[0].strip().lower()
        d[key] = [v.split(":", 1)[1].strip() for v in vals]
    return d


# 1) 시료 정보 → 원 연구의 군 이름(raw_group)
s1 = series_matrix("GSE175495")
s1["raw_group"] = s1.title.str.split("_").str[0] + "_" + s1.title.str.split("_").str[2]  # Old_AdiposeTissue 등
s1["sex"] = "M"  # 시리즈 설명: "healthy, active younger and older men"
s2 = series_matrix("GSE141432")
s2["raw_group"] = s2.condition
s3 = series_matrix("GSE166047")
s3["raw_group"] = s3.title.str.replace(r"\d+_S\d+$", "", regex=True)  # NW1_S1 → NW, OBT2D1_S11 → OBT2D
s3["sex"] = s3.gender.str[0].str.upper()
for d, g in ((s1, "GSE175495"), (s2, "GSE141432"), (s3, "GSE166047")):
    d["series"] = g
samples = pd.concat([d[["gsm", "series", "title", "raw_group", "sex"]] for d in (s1, s2, s3)], ignore_index=True)

# 2) 논문 군: 표 1 = GSE175495 고령 지방(OH), GSE141432 Lean(OH)·Obese T2D(OT2D), GSE166047 OBT2D(OT2D)
paper_group = {"Old_AdiposeTissue": "OH", "Lean": "OH", "Obese T2D": "OT2D", "OBT2D": "OT2D"}
samples["paper_group"] = samples.raw_group.map(paper_group)

# 3) NCBI 생성 count 행렬에 있는지
counts = {}
for g in SERIES:
    counts[g] = pd.read_csv(RAW / f"{g}_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="\t", index_col="GeneID")
assert counts["GSE175495"].index.equals(counts["GSE141432"].index) and counts["GSE175495"].index.equals(counts["GSE166047"].index)
in_matrix = set().union(*(c.columns for c in counts.values()))
samples["in_ncbi_matrix"] = samples.gsm.isin(in_matrix)
OUT.mkdir(parents=True, exist_ok=True)
samples.to_csv(OUT / "samples.tsv", sep="\t", index=False)

# 4) 표 1 대조: 데이터셋 × 군 별 GEO 시료 수, NCBI 행렬에 있는 수
paper_n = {("GSE175495", "OH"): 12, ("GSE141432", "OH"): 9, ("GSE141432", "OT2D"): 8, ("GSE166047", "OT2D"): 5}
tab = (samples.dropna(subset=["paper_group"])
       .groupby(["series", "paper_group"])
       .agg(geo=("gsm", "size"), ncbi_matrix=("in_ncbi_matrix", "sum"), male=("sex", lambda x: (x == "M").sum())))
tab["paper_table1"] = [paper_n[i] for i in tab.index]
print(tab, "\n")
print("GEO 에는 있으나 NCBI 행렬에 없는 분석 대상:",
      samples.query("paper_group.notna() and not in_ncbi_matrix")[["gsm", "title"]].values.tolist())

# 5) 분석 표본 = 논문 군 중 NCBI 행렬에 있는 시료
main = samples.query("paper_group.notna() and in_ncbi_matrix")
merged = pd.concat([counts[g][main.query("series == @g").gsm] for g in SERIES], axis=1)
merged.to_csv(OUT / "counts_main.tsv.gz", sep="\t")
print(f"\n분석 표본 {merged.shape[1]}명 (OH {sum(main.paper_group == 'OH')}, OT2D {sum(main.paper_group == 'OT2D')}),"
      f" 유전자 {merged.shape[0]:,}개 → {OUT / 'counts_main.tsv.gz'}")
