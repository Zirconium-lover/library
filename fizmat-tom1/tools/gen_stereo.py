#!/usr/bin/env python3
"""Стереографические сетки и проекции к п. 1.2 (кристаллографические проекции).

Проекция из южного полюса (0, 0, −1) на экваториальную плоскость:
X = x/(1 + z), Y = y/(1 + z); радиус основного круга r = 1 (масштаб
задаётся в TikZ).

  figures/st-wulff.tex      — сетка Вульфа: меридианы (большие круги через
                              ось y) и параллели (малые круги вокруг оси y)
                              через 10° (каждые 30° — жирнее), 2°-штрихи
                              по основному кругу и экватору
  figures/st-boldyrev.tex   — сетка Болдырева: концентрические окружности
                              R = tg(θ/2) через 10° и радиусы через 10°
  figures/st-cubic001.tex   — стандартная проекция (001) кубического
                              кристалла: полюса {100}, {110}, {111},
                              {210}, {211}, {221} и следы плоскостей
                              {100} и {110} (большие круги)

Запуск: python3 tools/gen_stereo.py"""
import math, os
import numpy as np
from itertools import product, permutations

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')


def proj(p):
    x, y, z = p / np.linalg.norm(p)
    return x / (1 + z), y / (1 + z)


def path(points):
    return ' -- '.join(f'({x:.4f},{y:.4f})' for x, y in points)


def wulff():
    lines = ['% создано tools/gen_stereo.py: сетка Вульфа (r = 1)']
    t = np.radians(np.linspace(-90, 90, 91))
    for a in range(10, 180, 10):                      # меридианы
        al = math.radians(a)
        u = np.array([math.cos(al), 0, math.sin(al)]); v = np.array([0, 1.0, 0])
        pts = [proj(math.cos(s) * u + math.sin(s) * v) for s in t]
        w = 'wulffB' if a % 30 == 0 else 'wulffA'
        lines.append(f'\\draw[{w}] {path(pts)};')
    s = np.radians(np.linspace(0, 180, 91))
    for b in range(10, 90, 10):                       # параллели
        be = math.radians(b)
        for sign in (1, -1):
            pts = [proj(np.array([math.sin(be) * math.cos(q), sign * math.cos(be),
                                  math.sin(be) * math.sin(q)])) for q in s]
            w = 'wulffB' if b % 30 == 0 else 'wulffA'
            lines.append(f'\\draw[{w}] {path(pts)};')
    lines.append('\\draw[wulffB] (-1,0) -- (1,0) (0,-1) -- (0,1);')
    for d in range(0, 360, 2):                        # штрихи по основному кругу
        a = math.radians(d); L = 0.035 if d % 10 else 0.06
        lines.append(f'\\draw[wulffA] ({math.cos(a):.4f},{math.sin(a):.4f}) -- '
                     f'({(1 - L) * math.cos(a):.4f},{(1 - L) * math.sin(a):.4f});')
    lines.append('\\draw[wulffC] (0,0) circle (1);')
    with open(os.path.join(OUT, 'st-wulff.tex'), 'w') as f:
        f.write('\n'.join(lines) + '\n')


def boldyrev():
    lines = ['% создано tools/gen_stereo.py: сетка Болдырева (r = 1)']
    for th in range(10, 90, 10):
        R = math.tan(math.radians(th) / 2)
        w = 'wulffB' if th % 30 == 0 else 'wulffA'
        lines.append(f'\\draw[{w}] (0,0) circle ({R:.4f});')
    for ph in range(0, 180, 10):
        a = math.radians(ph)
        w = 'wulffB' if ph % 30 == 0 else 'wulffA'
        lines.append(f'\\draw[{w}] ({math.cos(a):.4f},{math.sin(a):.4f}) -- '
                     f'({-math.cos(a):.4f},{-math.sin(a):.4f});')
    lines.append('\\draw[wulffC] (0,0) circle (1);')
    with open(os.path.join(OUT, 'st-boldyrev.tex'), 'w') as f:
        f.write('\n'.join(lines) + '\n')


def idx(v):
    return ''.join(f'\\br{{{-c}}}' if c < 0 else str(c) for c in v)


def cubic001():
    lines = ['% создано tools/gen_stereo.py: стандартная проекция (001) куба (r = 1)']
    # следы плоскостей {100} и {110}: большие круги, нормальные к n
    t = np.radians(np.linspace(0, 360, 361))
    for n, style in [((1, 0, 0), 'zone100'), ((0, 1, 0), 'zone100'),
                     ((1, 1, 0), 'zone110'), ((1, -1, 0), 'zone110'),
                     ((1, 0, 1), 'zone110'), ((1, 0, -1), 'zone110'),
                     ((0, 1, 1), 'zone110'), ((0, 1, -1), 'zone110')]:
        n = np.array(n, float); n /= np.linalg.norm(n)
        a = np.cross(n, [0, 0, 1.0])
        if np.linalg.norm(a) < 1e-9:
            continue
        a /= np.linalg.norm(a); b = np.cross(n, a)
        if b[2] < 0:
            b = -b
        pts = [proj(math.cos(s) * a + math.sin(s) * b) for s in np.radians(np.linspace(0, 180, 121))]
        lines.append(f'\\draw[{style}] {path(pts)};')
    # полюса
    fams = {'100': 'pole100', '110': 'pole110', '111': 'pole111',
            '210': 'poleH', '211': 'poleH', '221': 'poleH'}
    seen = set()
    for fam, style in fams.items():
        base = [int(c) for c in fam]
        for perm in set(permutations(base)):
            for signs in product((1, -1), repeat=3):
                v = tuple(p * s for p, s in zip(perm, signs))
                if v in seen or v[2] < 0 or all(c == 0 for c in v):
                    continue
                seen.add(v)
                X, Y = proj(np.array(v, float))
                lab = idx(v)
                ang = math.degrees(math.atan2(Y, X)) if (X or Y) else -90
                anchor_ang = ang + 180
                big = fam in ('100', '110', '111')
                fs = '\\scriptsize' if big else '\\tiny'
                lines.append(f'\\node[{style}] at ({X:.4f},{Y:.4f}) {{}};')
                if big:
                    lines.append(f'\\node[font={fs}, inner sep=1pt, anchor={anchor_ang:.0f}] '
                                 f'at ({X:.4f},{Y:.4f}) {{${lab}$}};')
    with open(os.path.join(OUT, 'st-cubic001.tex'), 'w') as f:
        f.write('\n'.join(lines) + '\n')


wulff(); boldyrev(); cubic001()
print('готово: st-wulff, st-boldyrev, st-cubic001')


# ---------------------------------------------------------------- оси симметрии
def axes_row():
    """Стереографические проекции эквивалентных направлений для поворотных
    (1, 2, 3, 4, 6) и инверсионных (1̄, 2̄=m, 3̄, 4̄, 6̄) осей, перпендикулярных
    плоскости проекции; исходное направление — точка под углом 20° к оси x
    на расстоянии 0,62 r. Точка — верхняя полусфера, крестик — нижняя."""
    lines = ['% создано tools/gen_stereo.py: действие осей симметрии']
    rho, phi0 = 0.62, 20.0
    entries = [(n, False) for n in (1, 2, 3, 4, 6)] + [(n, True) for n in (1, 2, 3, 4, 6)]
    for k, (n, inv) in enumerate(entries):
        cx = (k % 5) * 2.6; cy = 0 if not inv else -2.75
        lines.append(f'\\draw[atomLine, line width=0.6pt] ({cx},{cy}) circle (1);')
        pts = set()
        ang, up = phi0, True
        for _ in range(2 * n if inv else n):
            pts.add((round(ang % 360, 3), up))
            if inv:
                ang += 360 / n + 180; up = not up
            else:
                ang += 360 / n
        for a, u in sorted(pts):
            x = cx + rho * math.cos(math.radians(a)); y = cy + rho * math.sin(math.radians(a))
            if u:
                lines.append(f'\\fill[fmBlue] ({x:.3f},{y:.3f}) circle (0.075);')
            else:
                lines.append(f'\\node[fmRed, font=\\small] at ({x:.3f},{y:.3f}) {{$\\times$}};')
        lab = (f'$\\br{{{n}}}$' if n != 2 else '$\\br2=m$') if inv else f'${n}$'
        lines.append(f'\\node[font=\\small] at ({cx},{cy - 1.32}) {{{lab}}};')
    with open(os.path.join(OUT, 'st-axes.tex'), 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------- группа m3m
def m3m():
    """Элементы симметрии группы m3̄m на проекции (001): плоскости симметрии
    {100} и {110} (большие круги), оси 4 (квадраты) вдоль <100>, оси 3
    (треугольники) вдоль <111>, оси 2 (овалы) вдоль <110>."""
    lines = ['% создано tools/gen_stereo.py: элементы симметрии m-3m на проекции (001)']
    for n in [(1, 0, 0), (0, 1, 0), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]:
        n = np.array(n, float); n /= np.linalg.norm(n)
        a = np.cross(n, [0, 0, 1.0])
        a /= np.linalg.norm(a); b = np.cross(n, a)
        if b[2] < 0:
            b = -b
        pts = [proj(math.cos(s) * a + math.sin(s) * b) for s in np.radians(np.linspace(0, 180, 121))]
        lines.append(f'\\draw[mirror] {path(pts)};')
    sym = {(1, 0, 0): 'ax4', (1, 1, 0): 'ax2', (1, 1, 1): 'ax3'}
    seen = set()
    for base, style in sym.items():
        for perm in set(permutations(base)):
            for signs in product((1, -1), repeat=3):
                v = tuple(p * s for p, s in zip(perm, signs))
                if v in seen or v[2] < 0:
                    continue
                seen.add(v)
                X, Y = proj(np.array(v, float))
                rot = math.degrees(math.atan2(Y, X)) if (X or Y) else 0
                extra = f', rotate={rot:.0f}' if style == 'ax2' else ''
                lines.append(f'\\node[{style}{extra}] at ({X:.4f},{Y:.4f}) {{}};')
    with open(os.path.join(OUT, 'st-m3m.tex'), 'w') as f:
        f.write('\n'.join(lines) + '\n')


axes_row(); m3m()
print('готово: st-axes, st-m3m')
