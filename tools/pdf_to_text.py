"""Extract text from PDFs (first N pages) into .txt files next to them, using pypdf.

Usage: python tools/pdf_to_text.py <pdf_or_dir> [max_pages]
"""
import sys
from pathlib import Path

from pypdf import PdfReader


def convert(pdf: Path, max_pages: int) -> Path:
    reader = PdfReader(str(pdf))
    n = min(len(reader.pages), max_pages)
    parts = []
    for i in range(n):
        try:
            parts.append(f"\n===== page {i + 1} =====\n" + (reader.pages[i].extract_text() or ""))
        except Exception as e:  # noqa: BLE001
            parts.append(f"\n===== page {i + 1} (extraction error: {e}) =====\n")
    out = pdf.with_suffix(".txt")
    out.write_text("".join(parts), encoding="utf-8")
    return out


if __name__ == "__main__":
    target = Path(sys.argv[1])
    max_pages = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    pdfs = [target] if target.is_file() else sorted(target.glob("*.pdf"))
    for p in pdfs:
        out = convert(p, max_pages)
        print(f"{p.name}: {len(PdfReader(str(p)).pages)} pages -> {out.name} ({out.stat().st_size} bytes)")
