#!/usr/bin/env python3
"""Рисунки и численные данные к п. 1.7 (кристаллография границ зёрен).

  figures/gb-csl5.tex    — (а) «дихроматическая» картина: две квадратные
        решётки, повёрнутые на 36,87° вокруг [001], совпадающие узлы и
        ячейка РСУ Σ5; (б) симметричная граница наклона Σ5 (310)[001]
  figures/gb-energy.dat  — схематическая зависимость энергии границы
        наклона вокруг <110> ГЦК-металла от угла разворота (провалы при
        специальных разориентациях Σ3, Σ9, Σ11)
  figures/gb-mackenzie.dat — распределение углов разориентации (по
        минимальному углу) для случайно ориентированных кубических зёрен
        (распределение Маккензи; метод Монте-Карло)

Печатает: Σ для таблицы специальных разориентаций (проверка), предельные
углы разворота вокруг <100>, <110>, <111> и максимальный угол
разориентации кубических кристаллов.

Запуск из каталога fizmat-tom1:  python3 tools/gen_boundaries.py"""
import math, os
import numpy as np
from itertools import permutations, product

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')


def rot(axis, ang):
    u = np.array(axis, float); u /= np.linalg.norm(u)
    c, s = math.cos(ang), math.sin(ang)
    K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    return np.eye(3) + s * K + (1 - c) * K @ K


def cubic_ops():
    ops = []
    for p in permutations(range(3)):
        for sg in product((1, -1), repeat=3):
            M = np.zeros((3, 3))
            for i in range(3):
                M[i, p[i]] = sg[i]
            if np.linalg.det(M) > 0:
                ops.append(M)
    return ops                     # 24 поворота группы 432


OPS = cubic_ops()


def angle_of(R):
    return math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(R) - 1) / 2))))


def disorientation(R):
    return min(angle_of(S @ R) for S in OPS)


def sigma_of(R):
    for s in range(1, 200):
        if np.allclose(s * R, np.round(s * R), atol=1e-7 * s):
            return s
    return None


# ------------------------------------------------------------ таблица Σ
TABLE = [((1, 1, 1), 60.00), ((1, 0, 0), 36.87), ((1, 1, 1), 38.21), ((1, 1, 0), 38.94),
         ((1, 1, 0), 50.48), ((1, 0, 0), 22.62), ((1, 1, 1), 27.80), ((2, 1, 0), 48.19),
         ((1, 0, 0), 28.07), ((1, 1, 0), 26.53), ((1, 1, 1), 21.79), ((1, 0, 0), 16.26),
         ((1, 1, 0), 31.59), ((1, 0, 0), 43.60)]


def check_table():
    """Точный угол: tg(θ/2) = (n/m)·√N, N = h² + k² + l² (подбор m, n)."""
    print('Проверка таблицы РСУ (ось, угол -> Σ):')
    for ax, th in TABLE:
        N = sum(c * c for c in ax)
        cand = [(abs(math.degrees(2 * math.atan(n / m * math.sqrt(N))) - th), m, n)
                for m in range(1, 12) for n in range(1, m + 1) if math.gcd(m, n) == 1]
        err, m, n = min(cand)
        th_ex = math.degrees(2 * math.atan(n / m * math.sqrt(N)))
        S = sigma_of(rot(ax, math.radians(th_ex)))
        print(f'  {ax}  {th_ex:6.2f}°  Σ = {S}  (Брэндон: Δθ = {15 / math.sqrt(S):.1f}°)')


def limits():
    print('Предельные углы разворота:')
    for ax in [(1, 0, 0), (1, 1, 0), (1, 1, 1)]:
        lim = 0
        for th in np.arange(0.05, 180, 0.05):
            if disorientation(rot(ax, math.radians(th))) < th - 1e-6:
                break
            lim = th
        print(f'  {ax}: {lim:.2f}°')
    # максимум по осям внутри стандартного треугольника
    best = (0, None)
    for a in np.linspace(0, 1, 41):
        for b in np.linspace(0, a, int(41 * a) + 1):
            ax = np.array([1.0, a, b])
            lo, hi = 40.0, 70.0
            for _ in range(40):
                mid = (lo + hi) / 2
                if disorientation(rot(ax, math.radians(mid))) >= mid - 1e-9:
                    lo = mid
                else:
                    hi = mid
            if lo > best[0]:
                best = (lo, ax / np.linalg.norm(ax))
    print(f'  максимум: {best[0]:.2f}° вокруг {np.round(best[1] / best[1].max(), 3)}')


# ------------------------------------------------------------ Σ5: рисунок
def csl5():
    th = 2 * math.atan(1 / 3)                     # 36,87°
    c0, s0 = math.cos(th), math.sin(th)
    for v in ((1, 2), (-2, 1)):                   # узлы РСУ принадлежат повёрнутой решётке
        w = (c0 * v[0] + s0 * v[1], -s0 * v[0] + c0 * v[1])
        assert all(abs(x - round(x)) < 1e-9 for x in w)
    L = ['% создано tools/gen_boundaries.py: решётка совпадающих узлов Σ5 [001]']
    W = 3.25
    # (а) дихроматическая картина
    L.append('\\begin{scope}')
    L.append(f'\\clip (-{W},-{W}) rectangle ({W},{W});')
    c, s = math.cos(th), math.sin(th)
    # ячейки РСУ: базис (1,2), (-2,1) в решётке 1 (поворот против часовой стрелки)
    for i in range(-3, 4):
        for j in range(-3, 4):
            o = i * np.array([1, 2]) + j * np.array([-2, 1])
            L.append(f'\\draw[fmGreen!70!black, line width=0.5pt] ({o[0]},{o[1]}) -- ++(1,2) -- ++(-2,1);')
    L.append(f'\\fill[fmGreen!15] (0,0) -- (1,2) -- (-1,3) -- (-2,1) -- cycle;')
    L.append(f'\\draw[fmGreen!70!black, line width=1.1pt] (0,0) -- (1,2) -- (-1,3) -- (-2,1) -- cycle;')
    # узлы ПРН (DSC) внутри выделенной ячейки: базис (0,2; 0,4), (-0,4; 0,2)
    cell = np.array([[1, 2], [-2, 1]], float).T
    for i in range(-10, 11):
        for j in range(-10, 11):
            q = i * np.array([0.2, 0.4]) + j * np.array([-0.4, 0.2])
            u, v = np.linalg.solve(cell, q)
            if -1e-9 <= u <= 1 + 1e-9 and -1e-9 <= v <= 1 + 1e-9:
                L.append(f'\\fill[fmGreen!60!black] ({q[0]:.3f},{q[1]:.3f}) circle (0.045);')
    for i in range(-6, 7):
        for j in range(-6, 7):
            x, y = i, j
            if abs(x) <= W + 0.2 and abs(y) <= W + 0.2:
                L.append(f'\\fill[fmBlue] ({x},{y}) circle (0.09);')
            X, Y = c * i - s * j, s * i + c * j
            if abs(X) <= W + 0.2 and abs(Y) <= W + 0.2:
                L.append(f'\\draw[fmRed, line width=0.7pt] ({X:.4f},{Y:.4f}) circle (0.15);')
            if abs(X - round(X)) < 1e-6 and abs(Y - round(Y)) < 1e-6 and \
                    abs(X) <= W + 0.2 and abs(Y) <= W + 0.2:
                L.append(f'\\draw[black, line width=1.1pt] ({X:.4f},{Y:.4f}) circle (0.25);')
    L.append('\\end{scope}')
    L.append(f'\\draw[atomLine] (-{W},-{W}) rectangle ({W},{W});')
    # (б) симметричная граница наклона (310)
    L.append('\\begin{scope}[shift={(7.6,0)}]')
    L.append(f'\\clip (-{W},-{W}) rectangle ({W},{W});')
    a = th / 2
    atoms = []
    for i in range(-8, 9):
        for j in range(-8, 9):
            for sgn, side in ((-1, -1), (1, 1)):         # зерно A: x<0, поворот −θ/2
                ca, sa = math.cos(sgn * a), math.sin(sgn * a)
                X, Y = ca * i - sa * j, sa * i + ca * j
                if side * X >= -1e-9 and abs(X) <= W + 0.3 and abs(Y) <= W + 0.3:
                    atoms.append((X, Y, side))
    uniq = []
    for X, Y, sd in atoms:
        if any((X - u[0]) ** 2 + (Y - u[1]) ** 2 < 1e-8 for u in uniq):
            continue
        uniq.append((X, Y, sd))
    L.append(f'\\fill[fmBlue!7] (-{W},-{W}) rectangle (0,{W});')
    L.append(f'\\fill[fmRed!7] (0,-{W}) rectangle ({W},{W});')
    for X, Y, sd in uniq:
        if abs(X) < 1e-6:
            L.append(f'\\node[atom=2.3mm, fill=atomI] at ({X:.4f},{Y:.4f}) {{}};')
            L.append(f'\\draw[black, line width=1.1pt] ({X:.4f},{Y:.4f}) circle (0.25);')
        else:
            st = 'atom' if sd < 0 else 'atomB'
            L.append(f'\\node[{st}=2.3mm] at ({X:.4f},{Y:.4f}) {{}};')
    L.append('\\end{scope}')
    L.append(f'\\draw[atomLine] (7.6-{W},-{W}) rectangle (7.6+{W},{W});')
    L.append(f'\\draw[fmGray, dashed, line width=0.6pt] (7.6,-{W}) -- (7.6,{W});')
    with open(os.path.join(OUT, 'gb-csl5.tex'), 'w') as f:
        f.write('\n'.join(L) + '\n')


# ------------------------------------------------------------ энергия границы
def energy():
    rs = lambda x: x * (1 - math.log(x)) if 0 < x < 1 else (0.0 if x == 0 else 1.0)
    # (угол, относительная энергия в провале, ширина провала по Брэндону)
    cusps = [(0, 0.0, 15), (180, 0.0, 15), (70.53, 0.05, 15 / math.sqrt(3)),
             (109.47, 0.55, 15 / math.sqrt(3)), (50.48, 0.55, 15 / math.sqrt(11)),
             (129.52, 0.75, 15 / math.sqrt(11)), (38.94, 0.80, 15 / math.sqrt(9)),
             (141.06, 0.80, 15 / math.sqrt(9))]
    rows = []
    for th in np.linspace(0, 180, 1801):
        E = 1.0
        for t0, e0, w in cusps:
            E = min(E, e0 + (1 - e0) * rs(abs(th - t0) / w))
        rows.append((th, E))
    with open(os.path.join(OUT, 'gb-energy.dat'), 'w') as f:
        f.write('# создано tools/gen_boundaries.py: схема E(θ) для границ наклона <110> (E/E_пл)\n')
        f.write('theta E\n')
        for th, E in rows:
            f.write(f'{th:.2f} {E:.4f}\n')


# ------------------------------------------------------------ Маккензи
def mackenzie(N=400000, seed=7):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=(N, 4)); q /= np.linalg.norm(q, axis=1)[:, None]
    # операции симметрии как кватернионы
    def mat2quat(M):
        w = math.sqrt(max(0, 1 + np.trace(M))) / 2
        if w > 1e-6:
            return np.array([w, (M[2, 1] - M[1, 2]) / (4 * w), (M[0, 2] - M[2, 0]) / (4 * w),
                             (M[1, 0] - M[0, 1]) / (4 * w)])
        x = math.sqrt(max(0, 1 + M[0, 0] - M[1, 1] - M[2, 2])) / 2
        y = math.sqrt(max(0, 1 - M[0, 0] + M[1, 1] - M[2, 2])) / 2
        z = math.sqrt(max(0, 1 - M[0, 0] - M[1, 1] + M[2, 2])) / 2
        # знаки по недиагональным элементам
        if x > 1e-6:
            y = math.copysign(y, M[0, 1] + M[1, 0]); z = math.copysign(z, M[0, 2] + M[2, 0])
        elif y > 1e-6:
            z = math.copysign(z, M[1, 2] + M[2, 1])
        return np.array([0, x, y, z])
    S = np.array([mat2quat(M) for M in OPS])
    # cos(θ/2) = max |<q, s>|  (угол поворота q·s⁻¹)
    w = np.max(np.abs(q @ S.T), axis=1)
    ang = np.degrees(2 * np.arccos(np.clip(w, -1, 1)))
    print(f'Маккензи: максимум {ang.max():.2f}°, среднее {ang.mean():.1f}°')
    h, e = np.histogram(ang, bins=np.arange(0, 64, 1.0), density=True)
    with open(os.path.join(OUT, 'gb-mackenzie.dat'), 'w') as f:
        f.write('# создано tools/gen_boundaries.py: распределение Маккензи (плотность, 1/град)\n')
        f.write('theta p\n')
        for k in range(len(h)):
            f.write(f'{(e[k] + e[k + 1]) / 2:.1f} {h[k]:.5f}\n')
    # контроль: аналитический хвост — плотность у 45°
    return ang


if __name__ == '__main__':
    check_table(); limits(); csl5(); energy(); mackenzie()
    print('готово: gb-csl5.tex, gb-energy.dat, gb-mackenzie.dat')
