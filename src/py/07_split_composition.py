"""논문 Results 의 분할(군마다 65%, 검증 OH 6 · OT2D 5)과 같은 구성의 분할만 골라 허브 유전자 빈도를 다시 센다.

새 계산은 없다. 03_ml.R · 04_roc.R 의 시드 1~100 결과(03_seed_stability, 04_test_auc_100seeds_R)를 검증 구성으로 나눠 요약한다.
단순 무작위 분할을 검증 구성으로 조건부로 묶으면, 그 안에서는 군별 무작위 분할과 분포가 같다.
출력: results/tables/07_hub_by_split_composition.tsv
실행: python src/py/07_split_composition.py  (프로젝트 폴더에서)
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"


def rows(name):
    with open(T / name, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


hub = {}
for r in rows("03_seed_stability_vst_table1_sample_no.tsv"):
    hub.setdefault(r["seed"], []).append(set(r["genes"].split(",")) if r["genes"] else set())
hub = {k: set.intersection(*v) for k, v in hub.items()}
comp = {r["seed"]: (r["n_OH_test"], r["n_OT2D_test"]) for r in rows("04_test_auc_100seeds_R.tsv")}

out = []
for label, keep in [("전체 100분할", lambda c: True), ("검증 OH 6 · OT2D 5 (논문 Results 의 군별 분할 구성)", lambda c: c == ("6", "5"))]:
    seeds = [s for s in hub if keep(comp[s])]
    n = len(seeds)
    cnt = {}
    for s in seeds:
        for g in hub[s]:
            cnt[g] = cnt.get(g, 0) + 1
    top = sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[:3]
    out.append({"splits": label, "n": n,
                "AIM2": cnt.get("AIM2", 0), "FHOD3": cnt.get("FHOD3", 0),
                "both": sum({"AIM2", "FHOD3"} <= hub[s] for s in seeds),
                "exact": sum(hub[s] == {"AIM2", "FHOD3"} for s in seeds),
                "top3": ", ".join(f"{g} {c}" for g, c in top)})

with open(T / "07_hub_by_split_composition.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
    w.writeheader(); w.writerows(out)
for o in out:
    print(o)
