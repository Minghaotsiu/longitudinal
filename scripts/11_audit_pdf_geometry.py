from __future__ import annotations

import json
import sys
from pathlib import Path

import pdfplumber


def overlap_area(a: dict, b: dict) -> float:
    left = max(a["x0"], b["x0"])
    right = min(a["x1"], b["x1"])
    top = max(a["top"], b["top"])
    bottom = min(a["bottom"], b["bottom"])
    if left >= right or top >= bottom:
        return 0.0
    return (right - left) * (bottom - top)


def audit_pdf(path: Path) -> dict:
    with pdfplumber.open(path) as pdf:
        page = pdf.pages[0]
        words = page.extract_words(
            keep_blank_chars=False,
            use_text_flow=True,
            extra_attrs=["fontname", "size"],
        )
        clipped = []
        for word in words:
            if (
                word["x0"] < -0.1
                or word["top"] < -0.1
                or word["x1"] > page.width + 0.1
                or word["bottom"] > page.height + 0.1
            ):
                clipped.append(word["text"])

        text_overlaps = []
        for index, first in enumerate(words):
            first_area = max((first["x1"] - first["x0"]) * (first["bottom"] - first["top"]), 1e-9)
            for second in words[index + 1 :]:
                if second["x0"] > first["x1"]:
                    continue
                second_area = max(
                    (second["x1"] - second["x0"]) * (second["bottom"] - second["top"]),
                    1e-9,
                )
                shared = overlap_area(first, second)
                if shared / min(first_area, second_area) > 0.35:
                    text_overlaps.append(
                        {
                            "a": first["text"],
                            "b": second["text"],
                            "fraction_of_smaller_box": shared / min(first_area, second_area),
                        }
                    )

        return {
            "pdf": str(path),
            "page_width_pt": page.width,
            "page_height_pt": page.height,
            "word_count": len(words),
            "minimum_font_size_pt": min(
                (float(word.get("size", 0)) for word in words), default=None
            ),
            "clipped_words": clipped,
            "text_overlap_count": len(text_overlaps),
            "text_overlaps": text_overlaps[:30],
            "status": "PASS" if not clipped and not text_overlaps else "REVIEW REQUIRED",
        }


def main() -> int:
    if len(sys.argv) < 3:
        print("usage: audit_pdf_geometry.py OUTPUT.json FILE.pdf [FILE.pdf ...]")
        return 2
    output = Path(sys.argv[1])
    reports = [audit_pdf(Path(path)) for path in sys.argv[2:]]
    output.write_text(json.dumps(reports, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(reports, ensure_ascii=False, indent=2))
    return 0 if all(report["status"] == "PASS" for report in reports) else 1


if __name__ == "__main__":
    raise SystemExit(main())
