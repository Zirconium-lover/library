#!/usr/bin/env python3
"""Рисунки к п. 3.2 (зонная теория), рассчитанные численно.

  figures/kp-data.tex   — модель Кронига–Пенни (P = 3π/2): кривая f(αa),
                          разрешённые интервалы, закон дисперсии E(k)
  figures/bz-*.tex      — первые зоны Бриллюэна (ПК, ГЦК, ОЦК) в проекции,
                          невидимые рёбра — пунктиром
  figures/dos-sc.tex    — плотность состояний модели сильной связи для
                          простой кубической решётки
  figures/sq-contours.tex — изоэнергетические линии в первой зоне
                          квадратной решётки (та же модель, 2D)

Запуск: python3 tools/gen_bands.py"""
import itertools, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'figures')


def save(name, text):
    with open(os.path.join(OUT, name), 'w', encoding='utf-8') as f:
        f.write('% сгенерировано tools/gen_bands.py — не править вручную\n' + text)


def save_dat(name, pts, header='x y'):
    with open(os.path.join(OUT, name), 'w', encoding='utf-8') as f:
        f.write(header + '\n' + '\n'.join(' '.join(f'{v:.5f}' for v in p) for p in pts) + '\n')


def coords(pts, fmt='{:.4f}'):
    return ' '.join('(' + fmt.format(x) + ',' + fmt.format(y) + ')' for x, y in pts)


# ---------------------------------------------------------------- Крониг–Пенни
P = 1.5 * math.pi
XMAX = 4 * math.pi + 0.3


def f(x):
    return P * math.sin(x) / x + math.cos(x)


xs = np.linspace(1e-6, XMAX, 4001)
fx = np.array([f(x) for x in xs])
allowed = np.abs(fx) <= 1
# границы разрешённых интервалов (уточняем бисекцией)
def refine(a, b, target):
    for _ in range(60):
        m = 0.5 * (a + b)
        if (f(a) - target) * (f(m) - target) <= 0:
            b = m
        else:
            a = m
    return 0.5 * (a + b)

intervals = []
i = 0
while i < len(xs):
    if allowed[i]:
        j = i
        while j + 1 < len(xs) and allowed[j + 1]:
            j += 1
        lo = xs[i] if i == 0 else refine(xs[i - 1], xs[i], 1 if fx[i - 1] > 1 else -1)
        hi = xs[j] if j == len(xs) - 1 else refine(xs[j], xs[j + 1], 1 if fx[j + 1] > 1 else -1)
        intervals.append((lo, hi))
        i = j + 1
    else:
        i += 1
intervals = [iv for iv in intervals if iv[1] < XMAX - 0.05][:4]

out = []
# кривая f(αa), обрезанная сверху для графика
pts = [(x, f(x)) for x in np.linspace(0.3, XMAX, 600)]
save_dat('kp-f.dat', pts)
for n, (lo, hi) in enumerate(intervals, 1):
    out.append(f'\\def\\KPlo{"ABCD"[n-1]}{{{lo:.4f}}}\\def\\KPhi{"ABCD"[n-1]}{{{hi:.4f}}}')
# E(k): E в единицах ħ²/(2ma²) равна (αa)²; ka = arccos f(αa) + сдвиг на номер зоны
for n, (lo, hi) in enumerate(intervals, 1):
    pos, neg = [], []
    for x in np.linspace(lo, hi, 160):
        c = max(-1.0, min(1.0, f(x)))
        ka = math.acos(c)  # 0..π
        # расширенная схема зон: n-я зона занимает (n-1)π < |ka| < nπ
        kext = (n - 1) * math.pi + (ka if n % 2 == 1 else math.pi - ka)
        pos.append((kext / math.pi, x * x))
        neg.append((-kext / math.pi, x * x))
    save_dat(f'kp-band{n}.dat', sorted(pos))
    save_dat(f'kp-band{n}m.dat', sorted(neg))
save('kp-data.tex', '\n'.join(out) + '\n')
print('Крониг–Пенни: разрешённые интервалы αa:', [(round(a, 3), round(b, 3)) for a, b in intervals])


# ---------------------------------------------------------------- зоны Бриллюэна
def polyhedron(vertices, normals):
    """Грани: для каждой нормали — вершины, лежащие в её плоскости (max n·v)."""
    V = np.array(vertices, float)
    faces = []
    for n in normals:
        n = np.array(n, float)
        d = V @ n
        idx = [i for i in range(len(V)) if abs(d[i] - d.max()) < 1e-9]
        c = V[idx].mean(axis=0)
        # упорядочим по углу в плоскости грани
        u = V[idx[0]] - c
        u /= np.linalg.norm(u)
        w = np.cross(n / np.linalg.norm(n), u)
        idx.sort(key=lambda i: math.atan2((V[i] - c) @ w, (V[i] - c) @ u))
        faces.append((idx, n / np.linalg.norm(n)))
    return V, faces


def project_setup(az=28, el=18):
    a, e = math.radians(az), math.radians(el)
    # направление на наблюдателя
    view = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
    # экранные оси: x' — горизонталь, y' — вертикаль
    ex = np.array([-math.sin(a), math.cos(a), 0.0])
    ey = np.cross(view, ex)
    return view, ex, ey


def draw_bz(name, vertices, normals, points, scale=1.0):
    V, faces = polyhedron(vertices, normals)
    view, ex, ey = project_setup()
    P2 = [(scale * (v @ ex), scale * (v @ ey)) for v in V]
    edges = {}
    for fi, (idx, n) in enumerate(faces):
        vis = n @ view > 1e-9
        for a, b in zip(idx, idx[1:] + idx[:1]):
            key = tuple(sorted((a, b)))
            edges[key] = edges.get(key, False) or vis
    lines = []
    # заливка видимых граней (светлая, чтобы читалась объёмность)
    for idx, n in faces:
        if n @ view > 1e-9:
            shade = 8 + int(14 * (n @ np.array([0.3, -0.5, 0.8]) + 1) / 2)
            lines.append(f'\\fill[fmBlue!{shade}] ' + ' -- '.join(f'({P2[i][0]:.3f},{P2[i][1]:.3f})' for i in idx) + ' -- cycle;')
    for (a, b), vis in edges.items():
        st = 'bzfront' if vis else 'bzback'
        lines.append(f'\\draw[{st}] ({P2[a][0]:.3f},{P2[a][1]:.3f}) -- ({P2[b][0]:.3f},{P2[b][1]:.3f});')
    for lab, p in points.items():
        q = np.array(p, float)
        lines.append(f'\\coordinate (bz{lab}) at ({scale * (q @ ex):.3f},{scale * (q @ ey):.3f});')
    save(f'bz-{name}.tex', '\n'.join(lines) + '\n')


h = 0.5
sc_v = list(itertools.product((-h, h), repeat=3))
sc_n = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
draw_bz('sc', sc_v, sc_n, {'G': (0, 0, 0), 'X': (0.5, 0, 0), 'R': (0.5, 0.5, 0.5), 'M': (0.5, 0.5, 0)}, scale=1.5)

fcc_v = set()
for perm in itertools.permutations((0, 0.5, 1)):
    for s in itertools.product((-1, 1), repeat=3):
        fcc_v.add(tuple(p * q for p, q in zip(perm, s)))
fcc_v = sorted(fcc_v)
fcc_n = sc_n + list(itertools.product((-1, 1), repeat=3))
draw_bz('fcc', fcc_v, fcc_n, {'G': (0, 0, 0), 'X': (1, 0, 0), 'L': (0.5, 0.5, 0.5), 'W': (1, 0.5, 0)}, scale=1.0)

bcc_v = [v for v in set(itertools.permutations((1, 0, 0))) | set(itertools.permutations((-1, 0, 0)))]
bcc_v += list(itertools.product((-h, h), repeat=3))
bcc_n = []
for i, j in itertools.combinations(range(3), 2):
    for si, sj in itertools.product((-1, 1), repeat=2):
        n = [0, 0, 0]; n[i] = si; n[j] = sj
        bcc_n.append(tuple(n))
draw_bz('bcc', bcc_v, bcc_n, {'G': (0, 0, 0), 'H': (1, 0, 0), 'N': (0.5, 0.5, 0), 'P': (0.5, 0.5, 0.5)}, scale=1.0)
print('зоны Бриллюэна: ok')

# ---------------------------------------------------------------- плотность состояний (ПК, сильная связь)
M = 160
k = (np.arange(M) + 0.5) / M * 2 * np.pi - np.pi
cx = np.cos(k)
E = -2 * (cx[:, None, None] + cx[None, :, None] + cx[None, None, :])  # в единицах t
hist, edges = np.histogram(E.ravel(), bins=240, range=(-6, 6), density=True)
# лёгкое сглаживание
ker = np.array([1, 2, 3, 2, 1], float); ker /= ker.sum()
hs = np.convolve(hist, ker, mode='same')
centers = 0.5 * (edges[1:] + edges[:-1])
pts = [(-6.0, 0.0)] + [(c, h) for c, h in zip(centers, hs)] + [(6.0, 0.0)]
# свободные электроны с той же массой у дна зоны: E = -6 + k²  ->  g = sqrt(E+6)/(4π²)
free = [(e, math.sqrt(e + 6) / (4 * math.pi ** 2)) for e in np.linspace(-6, 1.0, 60)]
save_dat('dos-sc.dat', pts)
save_dat('dos-free.dat', free)
print('DOS: максимум', round(hs.max(), 4))

# ---------------------------------------------------------------- изоэнергетические линии (квадратная решётка)
# E = -2(cos kx a + cos ky a)  =>  cos ky = -E/2 - cos kx; точки сгущаются там,
# где касательная вертикальна (чтобы кривые были гладкими)
lines = []
def ky(x, Ev):
    return math.acos(max(-1.0, min(1.0, -Ev / 2 - math.cos(x))))
for Ev, col in [(-3.5, 'c1'), (-2.8, 'c2'), (-2.0, 'c3'), (-1.2, 'c4'), (-0.5, 'c5'), (0.0, 'c6'), (0.8, 'c7'), (2.0, 'c8')]:
    pieces = []
    if Ev < 0:      # замкнутая линия вокруг Γ: |kx| <= x0
        x0 = math.acos(-Ev / 2 - 1) if -Ev / 2 - 1 > -1 else math.pi
        xs_ = [x0 * math.sin(t) for t in np.linspace(-math.pi / 2, math.pi / 2, 90)]
        pieces.append(xs_)
    elif Ev == 0:   # квадрат, вписанный в зону
        pieces.append(list(np.linspace(-math.pi, math.pi, 3)))
    else:           # «карманы» в углах зоны: x1 <= |kx| <= π
        x1 = math.acos(1 - Ev / 2)
        xs_ = [x1 + (math.pi - x1) * (1 - math.cos(t)) for t in np.linspace(0, math.pi / 2, 60)]
        pieces.append(xs_)
        pieces.append([-x for x in xs_])
    for xs_ in pieces:
        up = [(x / math.pi, ky(x, Ev) / math.pi) for x in xs_]
        dn = [(x / math.pi, -ky(x, Ev) / math.pi) for x in xs_]
        mode = '' if Ev == 0 else '[smooth]'
        lines.append(f'\\draw[{col}] plot{mode} coordinates {{{coords(up, "{:.3f}")}}};')
        lines.append(f'\\draw[{col}] plot{mode} coordinates {{{coords(dn, "{:.3f}")}}};')
save('sq-contours.tex', '\n'.join(lines) + '\n')
print('контуры: ok')

# ---------------------------------------------------------------- перекрытие зон (плотность состояний двух зон)
# форма каждой зоны — та же плотность состояний ПК-модели, ширина зоны 1
def band(Eg, start):
    x = (Eg - start) * 12 - 6          # перевод в шкалу [-6, 6]
    return np.interp(x, centers, hs, left=0.0, right=0.0) * 12

grid = np.linspace(0, 2.3, 461)
for name, s2 in (('overlap', 0.75), ('gap', 1.25)):
    g1 = band(grid, 0.0)
    g2 = band(grid, s2)
    save_dat(f'dos-{name}.dat', list(zip(grid, g1, g2, g1 + g2)), header='E g1 g2 s')
print('перекрытие зон: ok')
