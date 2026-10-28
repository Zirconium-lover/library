"""Слияние biblio/parts/*.csv в biblio/catalog.csv (дедупликация по DOI и названию, стабильные id)."""
import csv, glob, os, re, sys
from difflib import SequenceMatcher

B = os.path.dirname(os.path.abspath(__file__))
COLS = ["id", "direction", "subtopics", "year", "title", "authors", "journal", "type", "doi",
        "doi_check", "access", "summary_ru", "useful_for", "values", "values_source", "priority",
        "overlap_reorient", "pdf_in_library", "url", "notes"]
MULTI = {"direction", "subtopics", "pdf_in_library"}
RANK = {"PDF": 4, "V2": 3, "CAT": 3, "V1": 2, "no-doi": 1, "": 0}


def ndoi(d):
    d = (d or "").strip().lower()
    d = re.sub(r"^(https?://)?(dx\.)?doi\.org/", "", d)
    return re.sub(r"^doi:\s*", "", d).rstrip(".")


def ntitle(t):
    return re.sub(r"[^a-zа-я0-9]+", " ", (t or "").lower()).strip()


def read(p):
    with open(p, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def union(a, b):
    out = []
    for x in (a or "").split(";") + (b or "").split(";"):
        x = x.strip()
        if x and x not in out:
            out.append(x)
    return ";".join(sorted(out))


def merge(dst, src):
    for c in COLS:
        if c in ("id", "overlap_reorient"):
            continue
        s = (src.get(c) or "").strip()
        if not s:
            continue
        if c in MULTI:
            dst[c] = union(dst.get(c), s)
        elif c == "doi_check":
            if RANK.get(s, 0) > RANK.get(dst.get(c, ""), 0):
                dst[c] = s
        elif c == "priority":
            dst[c] = min(dst.get(c) or "C", s)
        elif c in ("values", "notes", "useful_for") and dst.get(c) and s not in dst[c]:
            dst[c] = dst[c] + " | " + s
        elif not (dst.get(c) or "").strip():
            dst[c] = s


def main():
    cat_p = os.path.join(B, "catalog.csv")
    cat = read(cat_p) if os.path.exists(cat_p) else []
    reo = {ndoi(r["doi"]): r["n"] for r in read(os.path.join(B, "existing_reorientation.csv")) if r["doi"]}
    reo_t = {ntitle(r["title"]): r["n"] for r in read(os.path.join(B, "existing_reorientation.csv"))}
    by_doi = {ndoi(r["doi"]): r for r in cat if r["doi"]}
    by_t = {ntitle(r["title"]): r for r in cat}
    nxt = max([int(r["id"]) for r in cat] or [0]) + 1
    added = 0
    for p in sorted(glob.glob(os.path.join(B, "parts", "*.csv"))):
        if p.endswith((".needs.csv", ".unverified.csv")):
            continue
        for r in read(p):
            r["doi"] = ndoi(r.get("doi"))
            if not r.get("title"):
                continue
            t = ntitle(r["title"])
            hit = by_doi.get(r["doi"]) if r["doi"] else None
            if hit is None:
                hit = by_t.get(t)
            if hit is None and len(t) > 30:
                for k, v in by_t.items():
                    if abs(len(k) - len(t)) < 15 and SequenceMatcher(None, k, t).ratio() > 0.93 and \
                            (not r["doi"] or not v["doi"] or v["doi"] == r["doi"]):
                        hit = v
                        break
            if hit is not None:
                merge(hit, r)
                continue
            new = {c: "" for c in COLS}
            merge(new, r)
            new["doi"] = r["doi"]
            new["id"] = str(nxt)
            nxt += 1
            added += 1
            cat.append(new)
            if new["doi"]:
                by_doi[new["doi"]] = new
            by_t[t] = new
    for r in cat:
        r["overlap_reorient"] = reo.get(r["doi"], "") or reo_t.get(ntitle(r["title"]), "")
    cat.sort(key=lambda r: int(r["id"]))
    with open(cat_p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(cat)
    # needs_fulltext и unverified
    ids = {r["doi"]: r["id"] for r in cat if r["doi"]}
    ids_t = {ntitle(r["title"]): r["id"] for r in cat}
    needs, seen = [], set()
    for p in sorted(glob.glob(os.path.join(B, "parts", "*.needs.csv"))):
        for r in read(p):
            d = ndoi(r.get("doi"))
            key = (d or ntitle(r.get("title")), r.get("what_to_extract"))
            if key in seen:
                continue
            seen.add(key)
            needs.append({"id": ids.get(d) or ids_t.get(ntitle(r.get("title")), ""), "doi": d,
                          "title": r.get("title", ""), "what_to_extract": r.get("what_to_extract", "")})
    with open(os.path.join(B, "needs_fulltext.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, ["id", "doi", "title", "what_to_extract"])
        w.writeheader()
        w.writerows(needs)
    unv, seen = [], set()
    uc = ["title", "authors", "year", "journal", "why", "where_seen"]
    for p in sorted(glob.glob(os.path.join(B, "parts", "*.unverified.csv"))):
        for r in read(p):
            t = ntitle(r.get("title"))
            if t in seen or t in ids_t:
                continue
            seen.add(t)
            unv.append(r)
    with open(os.path.join(B, "candidates_unverified.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, uc, extrasaction="ignore")
        w.writeheader()
        w.writerows(unv)
    print(f"catalog: {len(cat)} (+{added}); needs: {len(needs)}; unverified: {len(unv)}")


if __name__ == "__main__":
    sys.exit(main())
