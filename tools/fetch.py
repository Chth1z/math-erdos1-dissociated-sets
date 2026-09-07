"""Fetch a URL and dump readable plain text (for literature checks).

Usage: python tools/fetch.py URL [OUTFILE]
"""
import html
import re
import sys
import urllib.request


def fetch_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
    raw = re.sub(r"(?is)<(script|style|noscript).*?</\1>", " ", raw)
    raw = re.sub(r"(?i)<br\s*/?>|</p>|</div>|</li>|</h\d>|</tr>", "\n", raw)
    text = re.sub(r"<[^>]+>", " ", raw)
    text = html.unescape(text)
    text = re.sub(r"[ \t\r\f\v]+", " ", text)
    text = re.sub(r"\n\s*\n+", "\n\n", text)
    return text.strip()


if __name__ == "__main__":
    url = sys.argv[1]
    text = fetch_text(url)
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w", encoding="utf-8") as f:
            f.write(text)
        print(f"saved {len(text)} chars to {sys.argv[2]}")
    else:
        sys.stdout.reconfigure(encoding="utf-8")
        print(text)
