#!/usr/bin/env python3
"""Генерирует TikZ-рисунок малоугловой границы наклона (стенки краевых
дислокаций) из суперпозиции упругих полей смещений краевых дислокаций
(изотропная упругость, ν = 0,3). Результат: figures/tilt-boundary.tex.

Запуск из каталога fizmat-tom1:  python3 tools/gen_tilt_boundary.py
"""
import math, sys

NU = 0.3
B = 1.0            # вектор Бюргерса = период решётки
H = 4              # расстояние между дислокациями в стенке (в периодах)
NDISL = 15         # число дислокаций в стенке (нечётное, симметрично)
XW, YW = 6.3, 6.3  # полуширина окна рисунка
R_ATOM = 0.2

def disp(x, y):
    """Смещение от одной краевой дислокации в начале координат; разрез по x<0."""
    r2 = x * x + y * y
    th = math.atan2(y, x)
    ux = B / (2 * math.pi) * (th + x * y / (2 * (1 - NU) * r2))
    uy = -B / (2 * math.pi) * ((1 - 2 * NU) / (4 * (1 - NU)) * math.log(r2)
                               + (x * x - y * y) / (4 * (1 - NU) * r2))
    return ux, uy

cores = [(0.0, H * (k - NDISL // 2)) for k in range(NDISL)]
atoms = []
for i in range(-12, 12):
    for j in range(-12, 12):
        x0, y0 = i + 0.5, j + 0.5
        ux = uy = 0.0
        for (cx, cy) in cores:
            dx, dy = disp(x0 - cx, y0 - cy)
            ux += dx; uy += dy
        x, y = x0 + ux, y0 + uy
        atoms.append((x, y))
# центрируем по вертикали (логарифмические члены дают общий сдвиг)
my = sum(a[1] for a in atoms if abs(a[0]) < 3) / max(1, sum(1 for a in atoms if abs(a[0]) < 3))
atoms = [(x, y - my) for (x, y) in atoms]
inwin = [a for a in atoms if abs(a[0]) < XW + 1 and abs(a[1]) < YW + 1]

out = []
out.append("% Сгенерировано tools/gen_tilt_boundary.py — не редактировать вручную")
out.append("\\begin{scope}")
out.append("\\clip (%.2f,%.2f) rectangle (%.2f,%.2f);" % (-XW, -YW, XW, YW))
# связи: пары атомов на расстоянии < 1.25
bonds = []
for a in range(len(inwin)):
    xa, ya = inwin[a]
    for b in range(a + 1, len(inwin)):
        xb, yb = inwin[b]
        d2 = (xa - xb) ** 2 + (ya - yb) ** 2
        if d2 < 1.25 ** 2:
            bonds.append("(%.3f,%.3f)--(%.3f,%.3f)" % (xa, ya, xb, yb))
out.append("\\draw[atomLine!55, line width=0.45pt] " + " ".join(bonds) + ";")
for (x, y) in inwin:
    out.append("\\fill[atomA, draw=atomLine, line width=0.35pt] (%.3f,%.3f) circle (%.2f);" % (x, y, R_ATOM))
out.append("\\end{scope}")
for (cx, cy) in cores:
    cyy = cy - my
    if abs(cyy) < YW - 0.3:
        out.append("\\dislT[fmRed]{(%.3f,%.3f)}{0}" % (cx, cyy - 0.05))
# опорные точки для размерных линий: две соседние дислокации стенки
out.append("\\coordinate (TB0) at (0,%.3f);" % (0 - my))
out.append("\\coordinate (TB1) at (0,%.3f);" % (H - my))
open("figures/tilt-boundary.tex", "w", encoding="utf-8").write("\n".join(out) + "\n")
print("atoms:", len(inwin), "bonds:", len(bonds), "theta(deg)=%.1f" % math.degrees(B / H))
