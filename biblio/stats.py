"""Обновляет раздел «Статистика» в biblio/README.md по catalog.csv."""
import csv, collections, os
B = os.path.dirname(os.path.abspath(__file__))
r = list(csv.DictReader(open(os.path.join(B, "catalog.csv"), encoding="utf-8")))
n = len(r)
c = collections.Counter(d for x in r for d in x["direction"].split(";"))
dc = collections.Counter(x["doi_check"] for x in r)
pr = collections.Counter(x["priority"] for x in r)
s = collections.Counter(t for x in r for t in x["subtopics"].split(";") if t)
cnt = lambda f: sum(1 for _ in csv.DictReader(open(os.path.join(B, f), encoding="utf-8")))
st = "\n".join(f"| {k} | {v} |" for k, v in sorted(s.items()))
sec = f"""## Статистика (на момент последней сборки)
Всего {n}. Направление 1 — {c['1']}, 2 — {c['2']}, 3 — {c['3']}.
doi_check: {', '.join(f'{k} {v} ({round(100*v/n)} %)' for k, v in dc.most_common())}.
Приоритет: A {pr['A']}, B {pr['B']}, C {pr['C']}.
Пересечение с каталогом пользователя — {sum(bool(x['overlap_reorient']) for x in r)}; PDF в library — {sum(bool(x['pdf_in_library']) for x in r)}.
needs_fulltext — {cnt('needs_fulltext.csv')} строк; candidates_unverified — {cnt('candidates_unverified.csv')}.

| подтема | записей |
|---|---|
{st}
"""
p = os.path.join(B, "README.md")
t = open(p, encoding="utf-8").read().split("## Статистика")[0]
open(p, "w", encoding="utf-8").write(t + sec)
print(sec.splitlines()[1:4])
