"""Stage 7: 원고 그림 PNG(600 dpi)를 투고용 JPEG 로 바꾼다(Sci Rep 은 첫 투고 그림을 JPEG 로 권장). 흰 바탕에 얹어 저장한다.

입력: manuscript/figures/Figure1_v7.png, Figure2_v7.png, SupplementaryFigure1_v1.png (src/r/07_manuscript_figures.R)
출력: 같은 이름의 .jpg
실행: python src/py/07_figure_jpeg.py  (프로젝트 폴더에서)
"""
from pathlib import Path

from PIL import Image

FIG = Path(__file__).resolve().parents[2] / "manuscript/figures"
for name in ("Figure1_v7", "Figure2_v7", "SupplementaryFigure1_v1"):
    png = FIG / f"{name}.png"
    im = Image.open(png).convert("RGBA")
    bg = Image.new("RGB", im.size, "white")
    bg.paste(im, mask=im.split()[3])
    bg.save(png.with_suffix(".jpg"), "JPEG", quality=95, subsampling=0, dpi=(600, 600))
    print(png.with_suffix(".jpg"))
