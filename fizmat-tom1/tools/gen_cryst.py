#!/usr/bin/env python3
"""Рисунки к главе 1, рассчитанные численно.

  figures/cr-young-polar.dat — модуль Юнга монокристаллов железа и
      вольфрама в зависимости от направления в плоскости (1-10):
      1/E(n) = S11 − 2(S11 − S12 − S44/2)(n1²n2² + n2²n3² + n3²n1²),
      n = (sinθ/√2, sinθ/√2, cosθ), θ — угол от [001] к [110];
      упругие постоянные (ГПа): Fe — c11 = 231, c12 = 135, c44 = 116;
      W — c11 = 522, c12 = 204, c44 = 161 (почти изотропен).
      Столбцы: x = E sinθ, y = E cosθ (полярная диаграмма).

Запуск: python3 tools/gen_cryst.py"""
import math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')


def compliances(c11, c12, c44):
    d = (c11 - c12) * (c11 + 2 * c12)
    return (c11 + c12) / d, -c12 / d, 1 / c44


def young(S, n):
    s11, s12, s44 = S
    n1, n2, n3 = n
    g = n1**2 * n2**2 + n2**2 * n3**2 + n3**2 * n1**2
    return 1 / (s11 - 2 * (s11 - s12 - s44 / 2) * g)


fe = compliances(231, 135, 116)
w = compliances(522, 204, 161)
th = np.radians(np.linspace(0, 360, 361))
rows = []
for t in th:
    n = (math.sin(t) / math.sqrt(2), math.sin(t) / math.sqrt(2), math.cos(t))
    ef, ew = young(fe, n), young(w, n)
    rows.append((ef * math.sin(t), ef * math.cos(t), ew * math.sin(t), ew * math.cos(t)))
with open(os.path.join(OUT, 'cr-young-polar.dat'), 'w') as f:
    f.write('xf yf xw yw\n')
    for r in rows:
        f.write(' '.join(f'{v:.2f}' for v in r) + '\n')
for name, S in (('Fe', fe), ('W', w)):
    e100 = young(S, (0, 0, 1)); e110 = young(S, (1 / math.sqrt(2), 1 / math.sqrt(2), 0))
    e111 = young(S, (1 / math.sqrt(3),) * 3)
    print(f'{name}: E<100> = {e100:.0f}, E<110> = {e110:.0f}, E<111> = {e111:.0f} ГПа')


# ---------------------------------------------------------------- 14 решёток Браве
def bravais():
    """figures/cr-bravais.tex — элементарные ячейки 14 решёток Браве в косоугольной
    проекции (x вправо, y «в глубину», z вверх, как стиль cubeproj)."""
    P = np.array([[1, 0.42, 0], [0, 0.30, 1]])        # проекция 3D -> 2D
    def pr(v):
        x, y = P @ np.asarray(v, float)
        return x, y
    s, c = 0.662, 0.749                                # ромбоэдр с α = 70°
    rot = lambda v, t: (v[0]*math.cos(t) - v[1]*math.sin(t), v[0]*math.sin(t) + v[1]*math.cos(t), v[2])
    rh = [rot((s*math.cos(math.radians(a)), s*math.sin(math.radians(a)), c), math.radians(-60)) for a in (0, 120, 240)]
    sc = 1.15
    cells = [
        ('триклинная $P$', [(1.1, 0, 0), (0.3, 1.0, 0), (0.25, 0.2, 1.1)], ''),
        ('моноклинная $P$', [(1.1, 0, 0), (0, 1.0, 0), (0.35, 0, 1.05)], ''),
        ('моноклинная $C$', [(1.1, 0, 0), (0, 1.0, 0), (0.35, 0, 1.05)], 'C'),
        ('тригональная $R$', [tuple(sc * x for x in v) for v in rh], ''),
        ('ромбическая $P$', [(1.25, 0, 0), (0, 0.9, 0), (0, 0, 1.0)], ''),
        ('ромбическая $C$', [(1.25, 0, 0), (0, 0.9, 0), (0, 0, 1.0)], 'C'),
        ('ромбическая $I$', [(1.25, 0, 0), (0, 0.9, 0), (0, 0, 1.0)], 'I'),
        ('ромбическая $F$', [(1.25, 0, 0), (0, 0.9, 0), (0, 0, 1.0)], 'F'),
        ('тетрагональная $P$', [(0.95, 0, 0), (0, 0.95, 0), (0, 0, 1.35)], ''),
        ('тетрагональная $I$', [(0.95, 0, 0), (0, 0.95, 0), (0, 0, 1.35)], 'I'),
        ('гексагональная $P$', [(1.1, 0, 0), (-0.55, 0.953, 0), (0, 0, 1.25)], ''),
        ('кубическая $P$', [(1.1, 0, 0), (0, 1.1, 0), (0, 0, 1.1)], ''),
        ('кубическая $I$', [(1.1, 0, 0), (0, 1.1, 0), (0, 0, 1.1)], 'I'),
        ('кубическая $F$', [(1.1, 0, 0), (0, 1.1, 0), (0, 0, 1.1)], 'F'),
    ]
    pos = [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (3, 1),
           (0, 2), (1, 2), (2, 2), (0, 3), (1, 3), (2, 3)]
    out = ['% создано tools/gen_cryst.py: 14 решёток Браве']
    for (name, (a, b, cc), cen), (col, row) in zip(cells, pos):
        a, b, cc = map(np.array, (a, b, cc))
        ox, oy = col * 3.4, -row * 3.05
        verts = {(i, j, k): i*a + j*b + k*cc for i in (0, 1) for j in (0, 1) for k in (0, 1)}
        # невидимая вершина — самая «дальняя»: максимум y − 0,3z − 0,4x
        hid = max(verts, key=lambda k: verts[k][1] - 0.3*verts[k][2] - 0.4*verts[k][0])
        def P2(v):
            x, y = pr(v); return f'({ox + x:.3f},{oy + y:.3f})'
        for k1 in verts:
            for d in range(3):
                if k1[d] == 0:
                    k2 = list(k1); k2[d] = 1; k2 = tuple(k2)
                    st = 'bravH' if hid in (k1, k2) else 'bravE'
                    out.append(f'\\draw[{st}] {P2(verts[k1])} -- {P2(verts[k2])};')
        for k1, v in verts.items():
            out.append(f'\\node[bravN] at {P2(v)} {{}};')
        extra = []
        if cen == 'I':
            extra = [(a + b + cc) / 2]
        elif cen == 'C':
            extra = [(a + b) / 2, (a + b) / 2 + cc]
        elif cen == 'F':
            extra = [(a + b) / 2, (a + b) / 2 + cc, (b + cc) / 2, (b + cc) / 2 + a, (a + cc) / 2, (a + cc) / 2 + b]
        for v in extra:
            out.append(f'\\node[bravC] at {P2(v)} {{}};')
        out.append(f'\\node[font=\\small, anchor=north] at ({ox + 0.75:.3f},{oy - 0.25:.3f}) {{{name}}};')
    with open(os.path.join(OUT, 'cr-bravais.tex'), 'w') as f:
        f.write('\n'.join(out) + '\n')


bravais()
print('готово: cr-bravais')


# ---------------------------------------------------------------- фуллерен C60
def c60():
    """figures/cr-c60.tex — молекула C60 (усечённый икосаэдр, ребро 2):
    вершины — циклические перестановки (0, ±1, ±3φ), (±1, ±(2+φ), ±2φ),
    (±φ, ±2, ±φ³); ребро — пара вершин на расстоянии 2; ортогональная
    проекция после поворота; передние пятиугольники закрашены."""
    phi = (1 + 5**0.5) / 2
    base = [(0, 1, 3*phi), (1, 2 + phi, 2*phi), (phi, 2, phi**3)]
    V = set()
    for b in base:
        for sx in (1, -1):
            for sy in (1, -1):
                for sz in (1, -1):
                    v = (b[0]*sx, b[1]*sy, b[2]*sz)
                    for k in range(3):
                        V.add(tuple(round(v[(i + k) % 3], 6) for i in range(3)))
    V = np.array(sorted(V))
    assert len(V) == 60
    a1, a2 = math.radians(18), math.radians(28)
    Rx = np.array([[1, 0, 0], [0, math.cos(a1), -math.sin(a1)], [0, math.sin(a1), math.cos(a1)]])
    Ry = np.array([[math.cos(a2), 0, math.sin(a2)], [0, 1, 0], [-math.sin(a2), 0, math.cos(a2)]])
    W = V @ (Ry @ Rx).T / 4.0
    edges = [(i, j) for i in range(60) for j in range(i + 1, 60) if abs(np.linalg.norm(V[i] - V[j]) - 2) < 1e-4]
    assert len(edges) == 90
    out = ['% создано tools/gen_cryst.py: молекула C60 (передние рёбра — сплошные)']
    # пятиугольники: вокруг 12 направлений (0, ±1, ±φ) и их перестановок
    dirs = set()
    for sy in (1, -1):
        for sz in (1, -1):
            d = (0, sy, sz * phi)
            for k in range(3):
                dirs.add(tuple(d[(i + k) % 3] for i in range(3)))
    for d in dirs:
        d = np.array(d) / np.linalg.norm(d)
        idx = np.argsort(-(V @ d))[:5]
        dw = (Ry @ Rx) @ d
        if dw[2] <= 0.05:
            continue
        pts = W[idx]
        cen = pts.mean(axis=0)
        ang = np.arctan2(pts[:, 1] - cen[1], pts[:, 0] - cen[0])
        pts = pts[np.argsort(ang)]
        out.append('\\fill[c60pent] ' + ' -- '.join(f'({p[0]:.3f},{p[1]:.3f})' for p in pts) + ' -- cycle;')
    for i, j in sorted(edges, key=lambda e: W[e[0]][2] + W[e[1]][2]):
        st = 'c60front' if W[i][2] + W[j][2] > 0 else 'c60back'
        out.append(f'\\draw[{st}] ({W[i][0]:.3f},{W[i][1]:.3f}) -- ({W[j][0]:.3f},{W[j][1]:.3f});')
    for p in W[np.argsort(W[:, 2])]:
        if p[2] > 0:
            out.append(f'\\fill[atomLine] ({p[0]:.3f},{p[1]:.3f}) circle (0.03);')
    with open(os.path.join(OUT, 'cr-c60.tex'), 'w') as f:
        f.write('\n'.join(out) + '\n')


# ---------------------------------------------------------------- мозаика Пенроуза
def penrose(n=5):
    """figures/cr-penrose.tex — мозаика Пенроуза из ромбов (подразбиение
    треугольников Робинсона, n итераций); тонкие ромбы 36°/144°, толстые
    72°/108°; выводится отношение числа толстых и тонких ромбов."""
    g = (1 + 5**0.5) / 2
    tri = []
    for i in range(10):
        B = complex(math.cos((2*i - 1)*math.pi/10), math.sin((2*i - 1)*math.pi/10))
        C = complex(math.cos((2*i + 1)*math.pi/10), math.sin((2*i + 1)*math.pi/10))
        if i % 2 == 0:
            B, C = C, B
        tri.append((0, 0j, B, C))
    for _ in range(n):
        new = []
        for col, A, B, C in tri:
            if col == 0:
                P = A + (B - A) / g
                new += [(0, C, P, B), (1, P, C, A)]
            else:
                Q = B + (A - B) / g; R = B + (C - B) / g
                new += [(1, R, C, A), (1, Q, R, B), (0, R, Q, A)]
        tri = new
    S = 3.4
    out = ['% создано tools/gen_cryst.py: мозаика Пенроуза']
    for col, A, B, C in tri:
        style = 'pthin' if col == 0 else 'pthick'
        out.append(f'\\fill[{style}] ({S*A.real:.3f},{S*A.imag:.3f}) -- ({S*B.real:.3f},{S*B.imag:.3f}) -- ({S*C.real:.3f},{S*C.imag:.3f}) -- cycle;')
    for col, A, B, C in tri:
        out.append(f'\\draw[pedge] ({S*C.real:.3f},{S*C.imag:.3f}) -- ({S*A.real:.3f},{S*A.imag:.3f}) -- ({S*B.real:.3f},{S*B.imag:.3f});')
    with open(os.path.join(OUT, 'cr-penrose.tex'), 'w') as f:
        f.write('\n'.join(out) + '\n')
    nthin = sum(1 for t in tri if t[0] == 0); nthick = len(tri) - nthin
    print(f'Пенроуз: треугольников {len(tri)}, толстых/тонких = {nthick/nthin:.4f} (τ = {g:.4f})')


c60(); penrose()
print('готово: cr-c60, cr-penrose')
