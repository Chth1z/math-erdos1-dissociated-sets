"""Summarise an erdosproblems.com list page (as dumped by fetch.py).

Prints, for each problem: status, problem number, tags, and the statement.
Usage: python tools/summarize_list.py literature/raw/prizes500.txt [--open-only]
"""
import re
import sys

STATUS_RE = re.compile(r"^(OPEN|SOLVED|DISPROVED|DECIDABLE|PROVED|PARTIALLY SOLVED|UNKNOWN)[A-Z (),/]*$")


def main(path: str, open_only: bool) -> None:
    lines = [l.strip() for l in open(path, encoding="utf-8").read().split("\n")]
    i = 0
    items = []
    while i < len(lines):
        l = lines[i]
        if STATUS_RE.match(l):
            status = l
            # statement is the first non-empty line after the status/explanation/prize lines
            j = i + 1
            stmt = None
            num = None
            tags = []
            seen_prize = False
            while j < len(lines) and j < i + 40:
                s = lines[j]
                if re.match(r"^-\s*\$\d+", s):
                    seen_prize = True
                    j += 1
                    continue
                m = re.match(r"^#(\d+)\s*:", s)
                if m:
                    num = int(m.group(1))
                    # tags follow on the next few non-empty lines until a long paragraph
                    k = j + 1
                    while k < len(lines) and k < j + 8:
                        t = lines[k].strip(" |")
                        if t and len(t) < 40 and not t.startswith("#"):
                            tags.append(t)
                        elif t:
                            break
                        k += 1
                    break
                if s and stmt is None and seen_prize:
                    stmt = s
                j += 1
            if num is not None:
                items.append((num, status, ", ".join(tags), stmt or ""))
            i = j
        i += 1
    for num, status, tags, stmt in items:
        if open_only and not status.startswith("OPEN"):
            continue
        stmt_short = stmt if len(stmt) <= 600 else stmt[:600] + " ..."
        print(f"#{num} [{status}] ({tags})\n    {stmt_short}\n")
    print(f"total listed: {len(items)}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main(sys.argv[1], "--open-only" in sys.argv)
