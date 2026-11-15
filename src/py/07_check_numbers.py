"""Stage 7: 원고 수치 점검. 원고(영문·한글)의 숫자가 결과 파일에서 나온 값과 같은지 본다.

1) 결과 파일에서 원고에 쓰는 값을 다시 계산해, 그 표기가 영문·한글 원고 본문에 둘 다 있는지 확인한다(없으면 실패).
2) 원고에 있는데 결과 파일 값도, 원 논문 값(아래 PAPER)도, 방법 설정값도 아닌 숫자를 모두 출력한다(사람이 확인).
3) 영문과 한글 원고에 나오는 숫자의 종류가 같은지 비교한다(번역·윤문 중 숫자가 바뀌거나 빠졌는지).
실행: python src/py/07_check_numbers.py  (프로젝트 폴더에서)  → 실패가 있으면 exit 1
"""
import csv
import gzip
import re
import statistics as st
import sys
from collections import Counter
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"
MS = ROOT / "manuscript"
EN = MS / "Manuscript_MattersArising_v1.7.md"
KO = MS / "Manuscript_MattersArising_v1.7_KO.md"


def rows(name):
    op = gzip.open if name.endswith(".gz") else open
    with op(T / name, "rt", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sci(p):  # 8.2 × 10⁻⁵ 꼴(그림·원고와 같은 표기)
    e = int(f"{p:e}".split("e")[1])
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{p / 10**e:.1f} × 10{str(e).translate(sup)}"


def pf(p):
    return sci(p) if p < 0.001 else f"{p:.2g}"


reg = {}  # 이름 → 원고 표기

# Stage 2: ComBat-seq·DESeq2 설정별 DEG 수
spec = {r["run"]: r for r in rows("02_deg_by_spec.tsv")}
reg["DEG 원문대로(R)"] = spec["R"]["deg"]
FIG = {"DEG 군 공변량(R)": spec["R_group"]["deg"]}   # 그림 1a 에만
reg["DEG 군 공변량 + protein-coding(R)"] = spec["R_group_pc"]["deg"]
reg["DEG Python 같은 설정"] = spec["py_group_pc"]["deg"]
reg["DEG R·Python 공통"] = spec["py_group_pc"]["paper_recovered"]
cmp_ = rows("02_paper_vs_R_group_pc.tsv")
rel = max(abs(float(r[k + "_R"]) - float(r[k])) / abs(float(r[k]))
          for r in cmp_ if r["padj_R"] not in ("", "NA") for k in ("baseMean", "log2FoldChange", "pvalue"))
reg["표 S1·S2 최대 상대 차"] = sci(rel)
rip = next(r for r in cmp_ if r["Symbol"] == "RIPPLY2")
reg["RIPPLY2 padj 논문"] = f"{float(rip['padj']):.4f}"
reg["RIPPLY2 padj 우리"] = f"{float(rip['padj_R']):.4f}"
reg["padj 배율"] = f"{st.median(float(r['padj_R']) / float(r['padj']) for r in cmp_ if r['padj_R'] not in ('', 'NA')):.3f}"

# Stage 3: 분할·선택
reg["λ 탐색 조합 수"] = str(len(rows("03_lambda_search.tsv")))
assert not any(r["match"] == "TRUE" for r in rows("03_lambda_search.tsv"))
by_seed = {}
for r in rows("03_seed_stability_vst_table1_sample_no.tsv"):
    by_seed.setdefault(r["seed"], []).append(set(r["genes"].split(",")) if r["genes"] else set())
hubs = [set.intersection(*v) for v in by_seed.values()]
reg["정확히 AIM2·FHOD3 분할 수"] = str(sum(h == {"AIM2", "FHOD3"} for h in hubs))
reg["둘 다 포함 분할 수"] = str(sum({"AIM2", "FHOD3"} <= h for h in hubs))
freq = {r[""]: float(r["hub"]) for r in rows("03_seed_hub_frequency.tsv")}
for g in ("STOX1", "FHOD3", "AIM2"):
    reg[f"공통 선택 비율 {g}"] = f"{round(freq[g] * 100)}%"
reg["학습 완전 분리 유전자 수"] = str(len(rows("03_train_perfect_auc_genes.tsv")))
pg = {r["item"]: r["value"] for r in rows("07_pergroup_summary.tsv")}   # Results 방식(군별) 분할, 07_pergroup_split.R
reg["군별 분할 FHOD3 비율"], reg["군별 분할 AIM2 비율"] = pg["FHOD3_pct"] + "%", pg["AIM2_pct"] + "%"
reg["군별 분할 둘 다"] = pg["seeds_both_AIM2_FHOD3"]
assert pg["top_gene"] == "FHOD3" and pg["seed100_hub"] == "CD52,FHOD3" and pg["seeds_exactly_AIM2_FHOD3"] == "0"
g1 = rows("07_gene_auc1_splits.tsv")
assert sum(int(r["splits_auc1"]) == 100 for r in g1) == 5                           # 100번 모두 AUC 1 = 31명 전체 완전 분리 DEG 5개('five DEGs including FHOD3')
assert sum(float(r["auc_all31"]) == 1 for r in g1) == 5
reg["RF 상위 10 R·Python 공통"] = next(r["R_py_shared"] for r in rows("03_table2_comparison.tsv") if r["method"] == "RF")

# Stage 4: 검증 AUC 100분할(배치 보정값)
a = rows("04_test_auc_100seeds_R.tsv")
reg["AIM2 AUC 1 분할 수"] = str(sum(float(r["AIM2"]) == 1 for r in a))
reg["FHOD3 AUC 1 분할 수"] = str(sum(float(r["FHOD3"]) == 1 for r in a))
FIG["DEG 검증 AUC 중앙값"] = f"{st.median(float(r['deg_median']) for r in a):.3f}"   # 원고 본문에는 없는 보조 값
eq1 = st.median(float(r["deg_eq1"]) for r in a)
reg["분할당 AUC 1 DEG 수(중앙값)"] = str(round(eq1 * 210))
reg["분할당 AUC 1 DEG 비율"] = f"{round(eq1 * 100)}%"

# Stage 5: 누수
lk = [r for r in rows("05_leakage_100seeds.tsv") if r["auc_clean"] not in ("", "NA")]
reg["누수 평가 분할 수"] = str(len(lk))
reg["누수 제외 분할 수"] = str(100 - len(lk))
reg["누수: 31명 선택 AUC"] = f"{st.median(float(r['auc_leaky']) for r in lk):.3f}"
reg["누수: 학습만 선택 AUC"] = f"{st.median(float(r['auc_clean']) for r in lk):.3f}"
reg["누수: 짝 차이 중앙값"] = f"{st.median(float(r['auc_leaky']) - float(r['auc_clean']) for r in lk):.3f}"
reg["누수: 낮아진 분할 수"] = str(sum(float(r["auc_clean"]) < float(r["auc_leaky"]) for r in lk))
tf = lambda v: str(v).strip().lower() == "true"
reg["누수 없이 고른 목록에 FHOD3"] = str(sum(tf(r["fhod3"]) for r in lk))
reg["누수 없이 고른 목록에 AIM2"] = str(sum(tf(r["aim2"]) for r in lk))
nd = [int(r["n_deg"]) for r in lk]
reg["학습만 선택 DEG 수 중앙값"], reg["최소"], reg["최대"] = str(round(st.median(nd))), str(min(nd)), str(max(nd))
comp = [r for r in rows("07_hub_by_split_composition.tsv") if r["n"] != "100"][0]
reg["검증 6·5 분할 수"] = comp["n"]
FIG["그중 AIM2"], FIG["그중 FHOD3"] = comp["AIM2"], comp["FHOD3"]   # 원고 본문에는 없는 보조 값
assert comp["both"] == "0"   # 원고는 낱말로 씀: "never selected together" / "한 번도 함께 선택되지 않았다"(아래에서 확인)

# Stage 5: 코호트 안 비교
OPT = {}  # 원고에 안 써도 되는 값(그림에만 있거나 범위로만 씀)
for r in rows("05_within_cohort_aim2_fhod3_tests_R.tsv"):
    k = f"{r['cohort']} {r['gene']} {r['comparison']}"
    p = float(r["p_wilcox"])
    in_text = (r["cohort"] == "GSE141432") != (r["comparison"] == "paper_like")          # GSE141432 비만·당뇨, GSE166047 논문형
    in_text = in_text or (r["cohort"] == "GSE166047" and r["comparison"] == "obesity")   # GSE166047 비만 · 당뇨 P 는 범위(0.15–0.55)로만
    (reg if in_text else OPT)[k + " P"] = "1.0" if p == 1 else pf(p)
    d = float(r["median_diff"])
    diff_in = in_text and not (r["cohort"] == "GSE166047" and r["comparison"] == "obesity")
    diff_in = diff_in or (r["cohort"] == "GSE166047" and r["gene"] == "AIM2" and r["comparison"] == "diabetes")   # AIM2 당뇨 비교 −1.07
    (reg if diff_in else OPT)[k + " 중앙값 차"] = f"{abs(d):.2f}"
for r in rows("05_within_cohort_deg210_summary.tsv"):
    k = f"{r['cohort']} DEG {r['comparison']}"
    tgt = reg if r["cohort"] == "GSE166047" or r["comparison"] != "paper_like" else OPT   # GSE141432 비만·당뇨 90%·94% 도 본문에
    tgt[k + " 같은 방향 %"] = f"{round(100 * int(r['same_direction']) / int(r['genes']))}%"
    OPT[k + " 유의"] = r["same_dir_padj05"]
reg["GSE166047 잴 수 있는 DEG"] = next(r["genes"] for r in rows("05_within_cohort_deg210_summary.tsv") if r["cohort"] == "GSE166047")
# 비만·당뇨 비교는 비만·비당뇨 5명을 공유하는 거울상: 비만 방향 유전자 중 당뇨에서 뒤집힌 수
d2 = {(r["Symbol"], r["comparison"]): tf(r["same_direction"]) for r in rows("05_within_cohort_deg210.tsv") if r["cohort"] == "GSE166047"}
ob = [g for (g, c), v in d2.items() if c == "obesity" and v]
reg["GSE166047 비만 방향 유전자"] = str(len(ob))
reg["그중 당뇨에서 뒤집힘"] = str(sum(not d2[(g, "diabetes")] for g in ob))
reg["5 대 5 최소 P"] = f"{2 / comb(10, 5):.2g}"
sx = [r for r in csv.DictReader(open(ROOT / "data/processed/samples.tsv", encoding="utf-8"), delimiter="\t") if r["raw_group"] == "Lean"]
reg["GSE141432 정상체중 여성"], reg["GSE141432 정상체중 남성"] = str(sum(r["sex"] == "F" for r in sx)), str(sum(r["sex"] == "M" for r in sx))

# GSE166047 군별 나이·BMI(GEO series matrix 의 군 요약값)
with gzip.open(ROOT / "data/raw/GSE166047_series_matrix.txt.gz", "rt", encoding="utf-8") as f:
    design = next(l for l in f if l.startswith("!Series_overall_design"))
for grp, lab in (("CTRL", "정상체중"), ("OBF", "비만"), ("OBT2D", "비만 당뇨")):
    m = re.search(rf"\({grp}[,;] age ([\d.]+) ± [\d.]+ years, BMI ([\d.]+) ±", design)
    reg[f"GSE166047 {lab} 나이"] = m.group(1)
    (OPT if lab == "정상체중" else reg)[f"GSE166047 {lab} BMI"] = m.group(2)

# 원 논문에 적힌 값(표 1·본문·그림), 방법 설정값, 참고문헌 서지 숫자는 결과 파일 밖이라 따로 둔다
PAPER = {"34", "21", "13", "12", "66.0", "5.0", "66", "5", "11", "9", "31", "18", "20", "210", "0.933", "1.000", "65", "35", "65%", "4",
         "0.01307997", "24.8", "23.75", "33.75", "38.1", "28/30", "30", "28", "6", "8", "2", "0.05", "1", "3",
         "2025", "15", "44904", "8"}
SETTINGS = {"100", "10", "500", "0.004400041", "0.5", "50%", "16", "1.021", "32", "1–210", "1–100", "95%", "99.97%", "209", "202", "200", "1.5", "1–12",
            "4.6.1", "3.60.0", "1.52.0", "5.0", "7.0.1", "0.9.33", "4.7.1.2", "1.19.1", "3.10", "0.5.2", "0.8.1",
            "1.7.2", "4.4.0", "38", "2026", "13.p13", "50612"}   # 50612 = 소속 우편번호


def body(p):  # 참고문헌 목록(서지 숫자)은 빼고 본다. 참고문헌 절만(다음 절 앞까지) 잘라 낸다
    return re.sub(r"\n## (?:References|참고문헌)\n.*?(?=\n## )", "\n", p.read_text(encoding="utf-8"), flags=re.S)


NUM = re.compile(r"\d+(?:\.\d+)?(?: × 10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)?%?")


def tokens(text):
    text = re.sub(r"`[^`]*`|^>.*$|<!--.*?-->", " ", text, flags=re.M)                # 코드 표기 · 안내 인용문 · 주석
    text = re.sub(r"GSE\d+|GRCh38\.p13|log₂|type 2 diabetes|제2형|sample\(31, 20\)|set\.seed\(100\)", " ", text)
    text = re.sub(r"\d+\.\d+\.\d+(?:\.\d+)?", " ", text)                                  # 소프트웨어 버전
    text = re.sub(r"https?://\S+", " ", text)                                            # 주소(DOI · GEO 링크)
    text = re.sub(r"\[\d+(?:[,–]\d+)*\]", " ", text)                                      # 참고문헌 번호 [1] · [7,8] · [2–4]
    text = re.sub(r"(?<![⁻⁰¹²³⁴⁵⁶⁷⁸⁹])[⁰¹²³⁴⁵⁶⁷⁸⁹]+", " ", text)                       # 위첨자(소속 번호 · m²; 지수 10⁻⁵ 는 남김)
    return Counter(m.group(0) for m in re.finditer(r"(?<![A-Za-z가-힣\d.])" + NUM.pattern + r"(?![A-Za-z\d])", text))


fail = 0
texts = {"EN": body(EN), "KO": body(KO)}
for lang, phrase in (("EN", "never selected together"), ("KO", "한 번도 함께 선택되지 않았")):
    if phrase not in texts[lang]:
        fail = 1
        print(f"[{lang}] 원고에 없음: 6·5 분할에서 둘 다 0번 ({phrase})")
for lang, txt in texts.items():
    toks = tokens(txt)                                  # 낱말 경계로 센 숫자(부분 문자열 일치로 통과하지 않게)
    for k, v in reg.items():
        if v not in toks:
            fail = 1
            print(f"[{lang}] 원고에 없음: {k} = {v}")
allowed = set(reg.values()) | set(OPT.values()) | set(FIG.values()) | PAPER | SETTINGS
for lang, txt in texts.items():
    stray = sorted(t for t in tokens(txt) if t not in allowed and t.rstrip("%") not in allowed)
    if stray:
        print(f"[{lang}] 결과·논문·설정 어디에도 없는 숫자(확인 필요): {stray}")
ce, ck = set(tokens(texts["EN"])), set(tokens(texts["KO"]))   # 같은 숫자를 한 언어에서만 되풀이하는 것은 허용(낱말 'six' · '11월' 등)
WORDS = {"6", "11", "1"}                                          # 영문에서 낱말로 쓴 수(six splits, November, a single)
diff = sorted((ce ^ ck) - WORDS)
if diff:
    fail = 1
    print("영문·한글 한쪽에만 있는 숫자:", diff)
print("결과 파일 값", len(reg), "개 확인 (본문 밖 보조 값", len(OPT) + len(FIG), "개)" + (" — 실패" if fail else " — 영문·한글 숫자 일치, 결과 파일 값 모두 원고에 있음"))
sys.exit(fail)
