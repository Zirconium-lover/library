#!/usr/bin/env python3
"""Генерирует TikZ-код элементарных ячеек ионных кристаллов (NaCl, CsCl,
ZnS-сфалерит, CaF2) для проекции cubeproj (см. fizmat.sty).
Атомы выводятся от дальних к ближним (по убыванию координаты y),
чтобы ближние перекрывали дальние.
Запуск: python3 tools/gen_ionic_cells.py  ->  figures/cell-*.tex"""
import itertools, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'figures')

def fmt(v):
    s = f"{v:.3f}".rstrip('0').rstrip('.')
    return s if s not in ('-0', '') else '0'

def pt(p):
    return '(' + ','.join(fmt(c) for c in p) + ')'

def in_cell(p):
    return all(-1e-9 <= c <= 1 + 1e-9 for c in p)

def fcc():
    base = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    pts = set()
    for b in base:
        for t in itertools.product((0, 1), repeat=3):
            p = tuple(round(b[i] + t[i], 6) for i in range(3))
            if in_cell(p):
                pts.add(p)
    return sorted(pts)

def write(name, atoms, bonds, bondstyle):
    """atoms: список (точка, стиль); bonds: список пар точек"""
    lines = [f"% сгенерировано tools/gen_ionic_cells.py — не править вручную"]
    for a, b in bonds:
        lines.append(f"\\draw[{bondstyle}] {pt(a)} -- {pt(b)};")
    # от дальних (большой y) к ближним; при равном y — снизу вверх
    for p, st in sorted(atoms, key=lambda q: (-q[0][1], q[0][2], q[0][0])):
        lines.append(f"\\node[{st}] at {pt(p)} {{}};")
    with open(os.path.join(OUT, f"cell-{name}.tex"), 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')

def dist(a, b):
    return math.dist(a, b)

# --- NaCl: узлы сетки с шагом 1/2; чётная сумма индексов — Cl, нечётная — Na
atoms, pts = [], []
for i, j, k in itertools.product(range(3), repeat=3):
    p = (i / 2, j / 2, k / 2)
    pts.append(p)
    atoms.append((p, 'anion' if (i + j + k) % 2 == 0 else 'cation'))
bonds = [(a, b) for a, b in itertools.combinations(pts, 2) if abs(dist(a, b) - 0.5) < 1e-6]
write('nacl', atoms, bonds, 'bondline')

# --- CsCl: Cl в вершинах, Cs в центре
atoms = [(p, 'anion') for p in itertools.product((0, 1), repeat=3)]
atoms.append(((0.5, 0.5, 0.5), 'cation'))
bonds = [((0.5, 0.5, 0.5), p) for p in itertools.product((0, 1), repeat=3)]
write('cscl', atoms, bonds, 'bondline')

# --- ZnS (сфалерит): S — ГЦК, Zn — в половине тетраэдрических пустот
S = fcc()
Zn = [(0.25, 0.25, 0.25), (0.75, 0.75, 0.25), (0.75, 0.25, 0.75), (0.25, 0.75, 0.75)]
atoms = [(p, 'anion') for p in S] + [(p, 'cation') for p in Zn]
bonds = [(z, s) for z in Zn for s in S if abs(dist(z, s) - math.sqrt(3) / 4) < 1e-6]
write('zns', atoms, bonds, 'bondthick')

# --- CaF2 (флюорит): Ca — ГЦК, F — во всех восьми тетраэдрических пустотах
Ca = fcc()
F = list(itertools.product((0.25, 0.75), repeat=3))
atoms = [(p, 'cation') for p in Ca] + [(p, 'anion') for p in F]
bonds = [(f, c) for f in F for c in Ca if abs(dist(f, c) - math.sqrt(3) / 4) < 1e-6]
write('caf2', atoms, bonds, 'bondline')
print('ok')
