"""Query the arXiv API and print verifiable metadata (id, title, authors, year, journal-ref, doi).

Usage: python tools/arxiv_query.py "<search_query>" [max_results]
Examples of search_query syntax: all:"distinct subset sums"    ti:dissociated AND cat:math.NT
"""
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

NS = {"a": "http://www.w3.org/2005/Atom", "arxiv": "http://arxiv.org/schemas/atom"}


def query(q: str, max_results: int = 10):
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode(
        {"search_query": q, "start": 0, "max_results": max_results, "sortBy": "relevance"})
    req = urllib.request.Request(url, headers={"User-Agent": "paper-spine-research/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = resp.read()
    root = ET.fromstring(data)
    out = []
    for e in root.findall("a:entry", NS):
        aid = e.find("a:id", NS).text.strip().split("/abs/")[-1]
        title = " ".join(e.find("a:title", NS).text.split())
        authors = [x.find("a:name", NS).text for x in e.findall("a:author", NS)]
        published = e.find("a:published", NS).text[:10]
        jref = e.find("arxiv:journal_ref", NS)
        doi = e.find("arxiv:doi", NS)
        cats = [c.attrib.get("term") for c in e.findall("a:category", NS)]
        summary = " ".join(e.find("a:summary", NS).text.split())
        out.append({
            "id": aid, "title": title, "authors": authors, "published": published,
            "journal_ref": jref.text if jref is not None else "", "doi": doi.text if doi is not None else "",
            "categories": cats, "summary": summary,
        })
    return out


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    q = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    for r in query(q, n):
        print(f"[{r['id']}] {r['title']}")
        print(f"    authors: {', '.join(r['authors'])} | published: {r['published']} | cats: {','.join(r['categories'])}")
        if r["journal_ref"] or r["doi"]:
            print(f"    journal-ref: {r['journal_ref']} | doi: {r['doi']}")
        print(f"    abstract: {r['summary'][:400]}{'...' if len(r['summary']) > 400 else ''}")
