#!/usr/bin/env python3
"""Рисунки к п. 1.6 (кристаллография пластической деформации).

  figures/pl-slip-stereo.tex — проекция кубического кристалла вдоль [100]
        (центр — [100]): следы плоскостей симметрии {100} и {110},
        стандартный треугольник [100]–[110]–[111], ось растяжения F и её
        поворот к направлению скольжения [101]; полюса плоскостей и
        направления первичной, сопряжённой, критической и поперечной систем
  figures/pl-ipf.tex         — стандартный треугольник крупно (обратные
        полюсные фигуры): точка [321] и линия G ⊥ F для текстуры <110>
  figures/pl-pf123.tex       — прямая полюсная фигура (100) аксиальной
        текстуры F || [123]: (а) F в центре проекции, (б) F на основном круге
  figures/pl-hcp-planes.tex  — плоскости скольжения ГПУ-кристалла
        (сечения гексагональной призмы)

Кроме того, печатает таблицу факторов Шмида 12 систем {111}<110> для
оси растяжения F (по умолчанию [531]).

Проекция — стереографическая из точки, противоположной оси взгляда;
радиус основного круга r = 1 (масштаб задаётся в TikZ).

Запуск из каталога fizmat-tom1:  python3 tools/gen_plastic.py"""
import math, os
import numpy as np
from itertools import product, permutations

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')


def unit(v):
    v = np.array(v, float)
    return v / np.linalg.norm(v)


def proj100(v):
    """Стереографическая проекция при взгляде вдоль [100]: X ← y, Y ← z."""
    x, y, z = unit(v)
    return y / (1 + x), z / (1 + x)


def path(points):
    return ' -- '.join(f'({x:.4f},{y:.4f})' for x, y in points)


def idx(v):
    return ''.join(f'\\br{{{-c}}}' if c < 0 else str(c) for c in v)


def great_circle(n, keep, prj, npts=181):
    """Видимая часть большого круга, нормального к n (keep(v) — условие
    видимости). Возвращает список полилиний."""
    n = unit(n)
    a = np.cross(n, [0.3, 0.5, 0.8]); a /= np.linalg.norm(a)
    b = np.cross(n, a)
    segs, cur = [], []
    for s in np.radians(np.linspace(0, 360, npts)):
        v = math.cos(s) * a + math.sin(s) * b
        if keep(v):
            cur.append(prj(v))
        elif cur:
            segs.append(cur); cur = []
    if cur:
        segs.append(cur)
    # склеить разрыв на стыке 0/360°
    if len(segs) > 1 and keep(a):
        segs[0] = segs[-1] + segs[0]; segs.pop()
    return segs


# ------------------------------------------------------------ факторы Шмида
PLANES = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]
DIRS = [(1, -1, 0), (1, 0, -1), (0, 1, -1), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
SYSTEMS = [(p, d) for p in PLANES for d in DIRS if np.dot(p, d) == 0]


def schmid(F, p, d):
    F = unit(F)
    return float(np.dot(F, unit(p)) * np.dot(F, unit(d)))


def schmid_table(F=(5, 3, 1)):
    print(f'Факторы Шмида для F = {F}:')
    for p, d in sorted(SYSTEMS, key=lambda s: -abs(schmid(F, *s))):
        print(f'  ({idx(p)})[{idx(d)}]  {abs(schmid(F, p, d)):.3f}')


# ------------------------------------------------------------ системы ГЦК
def slip_stereo(F=(5, 3, 1), RA=2.4, ZB=10.0, OB=(3.55, -2.05)):
    """Панель (а): вся проекция вдоль [100] радиусом RA (см);
    панель (б): окрестность стандартного треугольника с увеличением ZB,
    начало координат проекции — в точке OB."""
    keep = lambda v: v[0] >= -1e-9
    A = lambda X, Y: (RA * X, RA * Y)
    B = lambda X, Y: (OB[0] + ZB * X, OB[1] + ZB * Y)
    pa = lambda pts, T: ' -- '.join(f'({T(x, y)[0]:.4f},{T(x, y)[1]:.4f})' for x, y in pts)
    L = ['% создано tools/gen_plastic.py: системы скольжения ГЦК, проекция вдоль [100]']

    def arc(a, b, n=40):
        a, b = unit(a), unit(b)
        om = math.acos(np.dot(a, b))
        return [proj100((math.sin((1 - t) * om) * a + math.sin(t * om) * b) / math.sin(om))
                for t in np.linspace(0, 1, n)]
    tri = arc((1, 0, 0), (1, 1, 0)) + arc((1, 1, 0), (1, 1, 1)) + arc((1, 1, 1), (1, 0, 0))
    tri2 = arc((1, 0, 0), (1, 1, 1)) + arc((1, 1, 1), (1, 0, 1)) + arc((1, 0, 1), (1, 0, 0))
    zones = [((0, 1, 0), 'zone100'), ((0, 0, 1), 'zone100'), ((0, 1, 1), 'zone100'),
             ((0, 1, -1), 'zone100'), ((1, 1, 0), 'zone110'), ((1, -1, 0), 'zone110'),
             ((1, 0, 1), 'zone110'), ((1, 0, -1), 'zone110')]
    # ---- (а)
    for T, clip in ((A, None), (B, (-0.03, -0.03, 0.445, 0.445))):
        L.append('\\begin{scope}')
        if clip:
            x0, y0 = T(clip[0], clip[1]); x1, y1 = T(clip[2], clip[3])
            L.append(f'\\clip ({x0:.4f},{y0:.4f}) rectangle ({x1:.4f},{y1:.4f});')
            L.append(f'\\fill[fmAmber!8] {pa(tri2, T)} -- cycle;')
        L.append(f'\\fill[fmAmber!25] {pa(tri, T)} -- cycle;')
        for n, st in zones:
            for seg in great_circle(n, keep, proj100):
                L.append(f'\\draw[{st}] {pa(seg, T)};')
        L.append(f'\\draw[fmAmber!80!black, line width=0.8pt] {pa(tri, T)} -- cycle;')
        L.append('\\end{scope}')
        if clip:
            L.append(f'\\draw[atomLine, line width=0.5pt] ({x0:.4f},{y0:.4f}) rectangle ({x1:.4f},{y1:.4f});')
            # рамка увеличиваемой области на панели (а)
            a0 = A(clip[0], clip[1]); a1 = A(clip[2], clip[3])
            L.append(f'\\draw[atomLine!70, line width=0.4pt] ({a0[0]:.4f},{a0[1]:.4f}) rectangle ({a1[0]:.4f},{a1[1]:.4f});')
            L.append(f'\\draw[atomLine!50, densely dotted] ({a1[0]:.4f},{a1[1]:.4f}) -- ({x0:.4f},{y1:.4f}) '
                     f'({a1[0]:.4f},{a0[1]:.4f}) -- ({x0:.4f},{y0:.4f});')
        else:
            L.append(f'\\draw[atomLine, line width=0.7pt] (0,0) circle ({RA});')
    # полюса {100}, {110}, {111} на (а)
    seen = set()
    for fam, st in [((1, 0, 0), 'pole100'), ((1, 1, 0), 'pole110'), ((1, 1, 1), 'pole111')]:
        for perm in set(permutations(fam)):
            for sg in product((1, -1), repeat=3):
                v = tuple(c * s for c, s in zip(perm, sg))
                if v in seen or v[0] < 0:
                    continue
                if v[0] == 0 and (v[1] < 0 or (v[1] == 0 and v[2] < 0)):
                    continue
                seen.add(v)
                X, Y = A(*proj100(v))
                L.append(f'\\node[{st}] at ({X:.4f},{Y:.4f}) {{}};')
    # полюса на (б)
    for v, st, name in [((1, 0, 0), 'pole100', 'B100'), ((1, 1, 0), 'pole110', 'B110'),
                        ((1, 0, 1), 'pole110', 'B101'), ((1, 1, 1), 'pole111', 'B111'),
                        ((2, 1, 1), 'poleH', 'B211')]:
        X, Y = B(*proj100(v))
        L.append(f'\\node[{st}] ({name}) at ({X:.4f},{Y:.4f}) {{}};')
    # ось F: поворот к [101] до стороны [100]-[111], перескок, двойное скольжение к [211]
    Fv = unit(F); m = unit((1, 0, 1))
    om = math.acos(np.dot(Fv, m))
    pts, cross = [], None
    for t in np.linspace(0, 1, 2000):
        v = (math.sin((1 - t) * om) * Fv + math.sin(t * om) * m) / math.sin(om)
        pts.append(proj100(v))
        if v[2] >= v[1]:
            cross = v; break
    over = [proj100(unit(cross + 0.05 * t * (m - Fv))) for t in np.linspace(0, 1, 30)]
    a, b = unit(cross), unit((2, 1, 1))
    om2 = math.acos(np.dot(a, b))
    dbl = [proj100((math.sin((1 - t) * om2) * a + math.sin(t * om2) * b) / math.sin(om2))
           for t in np.linspace(0, 0.93, 30)]
    L.append(f'\\draw[fmRed, line width=1pt, -{{Stealth[length=2mm]}}] {pa(pts, B)};')
    L.append(f'\\draw[fmRed, line width=0.8pt, dashed] {pa(over, B)};')
    L.append(f'\\draw[fmGreen!70!black, line width=1.2pt, -{{Stealth[length=1.8mm]}}] {pa(dbl, B)};')
    for name, v, T in [('F', F, B), ('FA', F, A), ('Fcross', cross, B)]:
        X, Y = T(*proj100(v))
        L.append(f'\\coordinate ({name}) at ({X:.4f},{Y:.4f});')
    for name, v in {'nPrim': (1, 1, -1), 'dPrim': (1, 0, 1), 'nConj': (1, -1, 1), 'dConj': (1, 1, 0),
                    'nCrit': (1, 1, 1), 'dCrit': (1, 0, -1), 'nCross': (1, -1, -1),
                    'A100': (1, 0, 0), 'A010': (0, 1, 0), 'A001': (0, 0, 1), 'A110': (1, 1, 0)}.items():
        X, Y = A(*proj100(v))
        L.append(f'\\coordinate ({name}) at ({X:.4f},{Y:.4f});')
    print('пересечение со стороной [100]-[111]:', np.round(cross / cross[1], 2))
    with open(os.path.join(OUT, 'pl-slip-stereo.tex'), 'w') as f:
        f.write('\n'.join(L) + '\n')


def ipf():
    """Стандартный треугольник крупно: точка [321] (текстура <123>) и
    множество направлений G ⊥ F для текстуры <110>."""
    L = ['% создано tools/gen_plastic.py: стандартный треугольник [100]-[110]-[111] (r = 1)']

    def arc(a, b, n=40):
        a, b = unit(a), unit(b)
        om = math.acos(np.dot(a, b))
        return [proj100((math.sin((1 - t) * om) * a + math.sin(t * om) * b) / math.sin(om))
                for t in np.linspace(0, 1, n)]
    tri = arc((1, 0, 0), (1, 1, 0)) + arc((1, 1, 0), (1, 1, 1)) + arc((1, 1, 1), (1, 0, 0))
    L.append(f'\\fill[fmAmber!15] {path(tri)} -- cycle;')
    L.append(f'\\draw[atomLine, line width=0.6pt] {path(tri)} -- cycle;')
    # G ⊥ [110]: проверка — приведение случайных G в треугольник
    rng = np.random.default_rng(1)
    red = []
    for s in rng.uniform(0, 2 * math.pi, 2000):
        G = math.cos(s) * unit((1, -1, 0)) + math.sin(s) * np.array([0, 0, 1.0])
        g = np.sort(np.abs(G))[::-1]
        red.append(g)
    red = np.array(red)
    on_edge = np.all((np.abs(red[:, 0] - red[:, 1]) < 1e-9) | (np.abs(red[:, 1] - red[:, 2]) < 1e-9))
    assert on_edge, 'G ⊥ <110> должно лежать на сторонах треугольника'
    L.append(f'\\draw[fmBlue, line width=2.2pt, line cap=round] '
             f'{path(arc((1, 1, 0), (1, 1, 1)) + arc((1, 1, 1), (1, 0, 0))[1:])};')
    for name, v in [('T100', (1, 0, 0)), ('T110', (1, 1, 0)), ('T111', (1, 1, 1)),
                    ('T321', (3, 2, 1)), ('T210', (2, 1, 0)), ('T211', (2, 1, 1)),
                    ('T221', (2, 2, 1)), ('T310', (3, 1, 0)), ('T311', (3, 1, 1))]:
        X, Y = proj100(v)
        L.append(f'\\coordinate ({name}) at ({X:.4f},{Y:.4f});')
    with open(os.path.join(OUT, 'pl-ipf.tex'), 'w') as f:
        f.write('\n'.join(L) + '\n')


# ------------------------------------------------------------ ППФ (100) для F || [123]
def pf123():
    def projz(v):
        x, y, z = unit(v)
        if z < 0:
            x, y, z = -x, -y, -z
        return x / (1 + z), y / (1 + z)
    c = [k / math.sqrt(14) for k in (1, 2, 3)]
    alphas = [math.acos(ci) for ci in c]
    print('углы конусов ППФ(100), F||[123]:', [round(math.degrees(a), 1) for a in alphas])
    cols = ['fmBlue', 'fmRed', 'fmGreen!70!black']
    L = ['% создано tools/gen_plastic.py: ППФ (100) аксиальной текстуры F || [123] (r = 1)']
    # (а) F в центре
    L.append('\\begin{scope}')
    L.append('\\draw[atomLine, line width=0.7pt] (0,0) circle (1);')
    for a, col in zip(alphas, cols):
        L.append(f'\\draw[{col}, line width=1.1pt] (0,0) circle ({math.tan(a / 2):.4f});')
    # полюса {100} одного кристалла (поворот вокруг F на угол 25°)
    Fz = unit((1, 2, 3))
    e1 = unit(np.cross(Fz, (0, 0, 1))); e2 = np.cross(Fz, e1)
    Rm = np.array([e1, e2, Fz])          # строки: оси образца в осях кристалла
    phi = math.radians(25)
    Rz = np.array([[math.cos(phi), -math.sin(phi), 0], [math.sin(phi), math.cos(phi), 0], [0, 0, 1]])
    for k, col in zip(range(3), cols):
        v = Rz @ (Rm @ np.eye(3)[k])
        X, Y = projz(v)
        L.append(f'\\fill[{col}] ({X:.4f},{Y:.4f}) circle (0.045);')
    L.append('\\fill (0,0) circle (0.03);')
    L.append('\\end{scope}')
    # (б) F на основном круге (вдоль оси Y рисунка)
    L.append('\\begin{scope}[shift={(2.7,0)}]')
    L.append('\\draw[atomLine, line width=0.7pt] (0,0) circle (1);')
    F = np.array([0, 1.0, 0]); A = np.array([0, 0, 1.0]); B = np.array([1.0, 0, 0])
    for a, col in zip(alphas, cols):
        for aa in (a, math.pi - a):
            pts = []
            for t in np.radians(np.linspace(0, 180, 91)):
                v = math.cos(aa) * F + math.sin(aa) * (math.sin(t) * A + math.cos(t) * B)
                pts.append(projz(v))
            L.append(f'\\draw[{col}, line width=1.1pt] {path(pts)};')
    L.append('\\end{scope}')
    with open(os.path.join(OUT, 'pl-pf123.tex'), 'w') as f:
        f.write('\n'.join(L) + '\n')


# ------------------------------------------------------------ плоскости ГПУ
def hcp_planes(ca=1.633):
    """Сечения гексагональной призмы (вершины основания ±a1, ±a2, ±a3,
    высота c) плоскостями скольжения. Плоскость (h k i l) в координатах
    r = u1 a1 + u2 a2 + w c: h u1 + k u2 + l w = n."""
    a1 = np.array([1.0, 0, 0]); a2 = np.array([-0.5, math.sqrt(3) / 2, 0])
    c = np.array([0, 0, ca])
    hexv = [np.array([math.cos(math.radians(60 * k)), math.sin(math.radians(60 * k)), 0])
            for k in range(6)]
    verts = hexv + [v + c for v in hexv]
    edges = [(k, (k + 1) % 6) for k in range(6)] + [(k + 6, (k + 1) % 6 + 6) for k in range(6)] \
        + [(k, k + 6) for k in range(6)]
    # обратная матрица: декартовы -> (u1, u2, w)
    Minv = np.linalg.inv(np.array([a1, a2, c]).T)

    def section(h, k, l, n):
        f = lambda r: float(np.dot([h, k, l], Minv @ r)) - n
        pts = []
        for i, j in edges:
            fi, fj = f(verts[i]), f(verts[j])
            if abs(fi) < 1e-9:
                pts.append(verts[i])
            if fi * fj < -1e-12:
                t = fi / (fi - fj)
                pts.append(verts[i] + t * (verts[j] - verts[i]))
        # уникальные, упорядоченные по углу вокруг центра
        uniq = []
        for p in pts:
            if all(np.linalg.norm(p - q) > 1e-6 for q in uniq):
                uniq.append(p)
        cen = sum(uniq) / len(uniq)
        nrm = np.array([h, k, l]) @ Minv          # нормаль (декартова)
        e1 = unit(uniq[0] - cen); e2 = np.cross(unit(nrm), e1)
        uniq.sort(key=lambda p: math.atan2(np.dot(p - cen, e2), np.dot(p - cen, e1)))
        return uniq

    # (4 индекса) -> (h, k, l) и сдвиг n, подобранный для наглядности
    # варианты выбраны так, чтобы плоскость не была видна «с ребра» в cubeproj
    planes = [('0001', (0, 0, 1), 0.5), ('01\\br10', (0, 1, 0), 0.0),
              ('0\\br111', (0, -1, 1), 1.0), ('0\\br112', (0, -1, 2), 1.0),
              ('11\\br21', (1, 1, 1), 1.0), ('11\\br22', (1, 1, 2), 1.0)]
    L = [f'% создано tools/gen_plastic.py: плоскости скольжения ГПУ (cubeproj; c/a = {ca:.3f})']
    for m, (lab, (h, k, l), n) in enumerate(planes):
        dx = (m % 3) * 3.5; dz = -(m // 3) * 3.0
        L.append(f'\\begin{{scope}}[shift={{({dx},0,{dz})}}]')
        P = lambda v: f'({v[0]:.4f},{v[1]:.4f},{v[2]:.4f})'
        hidden = [(1, 2), (2, 3), (2, 8)]          # невидимые рёбра (вершина 120° внизу)
        for i, j in hidden:
            L.append(f'\\draw[atomLine!45, dashed] {P(verts[i])} -- {P(verts[j])};')
        poly = section(h, k, l, n)
        L.append('\\fill[fmRed!35, opacity=0.85] ' +
                 ' -- '.join(f'({p[0]:.4f},{p[1]:.4f},{p[2]:.4f})' for p in poly) + ' -- cycle;')
        L.append('\\draw[fmRed!80!black, line width=0.7pt] ' +
                 ' -- '.join(f'({p[0]:.4f},{p[1]:.4f},{p[2]:.4f})' for p in poly) + ' -- cycle;')
        for i, j in edges:
            if (i, j) not in hidden:
                L.append(f'\\draw[atomLine, line width=0.5pt] {P(verts[i])} -- {P(verts[j])};')
        L.append(f'\\node[font=\\small, anchor=north] at (0,-1.05,0) {{${m + 1}$: $({lab})$}};')
        L.append('\\end{scope}')
    with open(os.path.join(OUT, 'pl-hcp-planes.tex'), 'w') as f:
        f.write('\n'.join(L) + '\n')


if __name__ == '__main__':
    schmid_table()
    slip_stereo(); ipf(); pf123(); hcp_planes()
    print('готово: pl-slip-stereo, pl-ipf, pl-pf123, pl-hcp-planes')
