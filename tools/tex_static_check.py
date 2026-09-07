"""Static checks for a single-file LaTeX manuscript (no TeX engine required).

Checks: non-ASCII bytes, brace balance, environment balance, \\ref/\\eqref targets,
\\cite keys vs \\bibitem keys, duplicate labels, uncited bibitems.
Usage: python tools/tex_static_check.py path/to/main.tex
"""
import re
import sys
from pathlib import Path


def strip_comments(text: str) -> str:
    out = []
    for line in text.split("\n"):
        i = 0
        while True:
            j = line.find("%", i)
            if j == -1:
                out.append(line)
                break
            if j > 0 and line[j - 1] == "\\":
                i = j + 1
                continue
            out.append(line[:j])
            break
    return "\n".join(out)


def main(path: str) -> int:
    raw = Path(path).read_bytes()
    problems: list[str] = []
    for lineno, line in enumerate(raw.split(b"\n"), start=1):
        bad = [b for b in line if b > 127]
        if bad:
            problems.append(f"line {lineno}: non-ASCII bytes {bad[:5]}")
    text = strip_comments(raw.decode("utf-8", errors="replace"))

    depth = 0
    for lineno, line in enumerate(text.split("\n"), start=1):
        for k, ch in enumerate(line):
            if ch in "{}" and (k == 0 or line[k - 1] != "\\"):
                depth += 1 if ch == "{" else -1
                if depth < 0:
                    problems.append(f"line {lineno}: closing brace without opening")
                    depth = 0
    if depth != 0:
        problems.append(f"unbalanced braces: final depth {depth}")

    stack: list[tuple[str, int]] = []
    for m in re.finditer(r"\\(begin|end)\{([^}]*)\}", text):
        lineno = text.count("\n", 0, m.start()) + 1
        if m.group(1) == "begin":
            stack.append((m.group(2), lineno))
        else:
            if not stack or stack[-1][0] != m.group(2):
                problems.append(f"line {lineno}: \\end{{{m.group(2)}}} does not match {stack[-1] if stack else None}")
                if stack:
                    stack.pop()
            else:
                stack.pop()
    for env, lineno in stack:
        problems.append(f"line {lineno}: \\begin{{{env}}} never closed")

    labels = re.findall(r"\\label\{([^}]*)\}", text)
    dup = {l for l in labels if labels.count(l) > 1}
    if dup:
        problems.append(f"duplicate labels: {sorted(dup)}")
    label_set = set(labels)
    for m in re.finditer(r"\\(?:eq)?ref\{([^}]*)\}", text):
        if m.group(1) not in label_set:
            lineno = text.count("\n", 0, m.start()) + 1
            problems.append(f"line {lineno}: undefined reference {m.group(1)}")

    bibitems = re.findall(r"\\bibitem\{([^}]*)\}", text)
    cited: set[str] = set()
    for m in re.finditer(r"\\cite\{([^}]*)\}", text):
        for key in m.group(1).split(","):
            cited.add(key.strip())
    for key in sorted(cited - set(bibitems)):
        problems.append(f"cited key without bibitem: {key}")
    for key in bibitems:
        if key not in cited:
            problems.append(f"bibitem never cited: {key}")

    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", text, re.S)
    if abstract and "\\cite" in abstract.group(1):
        problems.append("abstract contains \\cite (E-JC requires a self-contained abstract)")

    sections = re.findall(r"^\\section\{", text, re.M)
    print(f"file: {path}")
    print(f"lines: {text.count(chr(10)) + 1}; numbered sections: {len(sections)}; labels: {len(label_set)}; bibitems: {len(bibitems)}; cited keys: {len(cited)}")
    if problems:
        print("PROBLEMS:")
        for p in problems:
            print("  -", p)
        return 1
    print("OK: no static problems found")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
