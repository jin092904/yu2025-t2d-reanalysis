"""Stage 0: 원자료와 논문 파일을 내려받는다.

논문 Methods: "Raw count matrices and the human genome annotation table provided by
the National Center for Biotechnology Information" → GEO 의 NCBI 생성 count 를 쓴다.
실행: python src/py/00_download.py  (프로젝트 폴더에서)
"""
import hashlib
import time
from datetime import datetime
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
RAW, PAPER = ROOT / "data/raw", ROOT / "paper"
GEO = "https://www.ncbi.nlm.nih.gov/geo/download/"
SERIES = ["GSE175495", "GSE141432", "GSE166047"]

files = []  # (저장 위치, URL)
for g in SERIES:
    # GEO 다운로드 링크(?type=rnaseq_counts&acc=...)가 넘겨주는 FTP 위치
    files.append((RAW / f"{g}_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                  f"https://ftp.ncbi.nlm.nih.gov/geo/rnaseq_counts/{g[:6]}nnn/{g}_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
    files.append((RAW / f"{g}_series_matrix.txt.gz",
                  f"https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:6]}nnn/{g}/matrix/{g}_series_matrix.txt.gz"))
files.append((RAW / "Human.GRCh38.p13.annot.tsv.gz",
              f"{GEO}?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"))
# 원 연구자가 올린 지방조직 count (표본 수 12 vs 9 확인용)
files.append((RAW / "GSE175495_Adipose_simple_counts.txt.gz",
              "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175495/suppl/GSE175495_Adipose_simple_counts.txt.gz"))
ESM = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-025-29141-9/MediaObjects/"
files += [
    (PAPER / "s41598-025-29141-9.pdf", "https://www.nature.com/articles/s41598-025-29141-9.pdf"),
    (PAPER / "PMC12748723.xml", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=12748723&rettype=xml"),
    (PAPER / "41598_2025_29141_MOESM1_ESM.docx", ESM + "41598_2025_29141_MOESM1_ESM.docx"),
    (PAPER / "41598_2025_29141_MOESM2_ESM.xlsx", ESM + "41598_2025_29141_MOESM2_ESM.xlsx"),
]

s = requests.Session()
# GEO 주석표는 다운로드 페이지 쿠키와 Referer 가 없으면 reCAPTCHA 페이지를 돌려준다.
# 반대로 GEO 는 브라우저 User-Agent 에 reCAPTCHA 를 걸고, nature.com 은 브라우저 User-Agent 가 필요하다.
ref = f"{GEO}?acc=GSE141432"
s.get(ref, timeout=60)
BROWSER = {"User-Agent": "Mozilla/5.0"}

rows = ["file\tbytes\tmd5\turl\tdownloaded"]
for path, url in files:
    path.parent.mkdir(parents=True, exist_ok=True)
    hdr = BROWSER if ("nature.com" in url or "springer.com" in url) else {"Referer": ref}
    for attempt in range(3):  # 서버가 가끔 500 을 돌려준다
        r = s.get(url, headers=hdr, timeout=300)
        if r.ok:
            break
        time.sleep(10)
    r.raise_for_status()
    body = r.content
    # gzip·PDF·docx/xlsx(zip)·XML 이 아니면 오류 페이지를 받은 것이다
    assert body[:2] == b"\x1f\x8b" or body[:4] in (b"%PDF", b"PK\x03\x04") or body[:5] == b"<?xml", f"{path.name}: 파일이 아님"
    path.write_bytes(body)
    rows.append(f"{path.relative_to(ROOT)}\t{len(body)}\t{hashlib.md5(body).hexdigest()}\t{url}\t{datetime.now():%Y-%m-%d %H:%M}")
    print(f"{path.relative_to(ROOT)}  {len(body):,} bytes")

(RAW / "MANIFEST.tsv").write_text("\n".join(rows) + "\n")
print("MANIFEST:", RAW / "MANIFEST.tsv")
