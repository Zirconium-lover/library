#!/usr/bin/env python3
"""Рисунки к п. 3.4 (диффузия), рассчитанные численно.

  figures/dif-erf-*.dat      — профили C/C0 = ½[1 + erf(x/2√Dt)] и
                               разложение на гауссовы «слои»
  figures/dif-matano.dat     — профиль взаимной диффузии при D(C) = D0·e^{βC}
                               (численное решение), figures/dif-matano.tex — x_M
  figures/dif-fisher.tex     — изоконцентраты в модели Фишера
  figures/dif-walk.tex       — двумерное случайное блуждание

Запуск: python3 tools/gen_diffusion.py"""
import math, os, random
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'figures')


def save(name, text):
    with open(os.path.join(OUT, name), 'w', encoding='utf-8') as f:
        f.write('% сгенерировано tools/gen_diffusion.py — не править вручную\n' + text)


def save_dat(name, cols, header):
    rows = zip(*cols)
    with open(os.path.join(OUT, name), 'w', encoding='utf-8') as f:
        f.write(header + '\n' + '\n'.join(' '.join(f'{v:.5f}' for v in r) for r in rows) + '\n')


# ---------------------------------------------------------------- erf-профили
x = np.linspace(-4, 4, 161)
cols, hdr = [x], ['x']
for L in (0.5, 1.0, 2.0):           # L = 2√(Dt) в условных единицах
    cols.append(0.5 * (1 + np.array([math.erf(v / L) for v in x])))
    hdr.append(f'L{int(L*10)}')
save_dat('dif-erf.dat', cols, ' '.join(hdr))

# разложение на «слои»: правая половина (x > 0) разбита на полоски ширины 0,5,
# каждая расплывается в гауссиану; сумма даёт erf-профиль (L = 1)
L = 1.0
cols, hdr = [x], ['x']
total = np.zeros_like(x)
for k, a in enumerate(np.arange(0.25, 6.0, 0.5)):
    g = 0.5 / (L * math.sqrt(math.pi)) * np.exp(-((x - a) / L) ** 2)
    total += g
    if k < 6:
        cols.append(g); hdr.append(f'g{k}')
cols.append(total); hdr.append('sum')
save_dat('dif-erf-sum.dat', cols, ' '.join(hdr))

# ---------------------------------------------------------------- Больцман–Матано
# ∂C/∂t = ∂/∂x(D(C)∂C/∂x), D = D0·exp(βC); C = 0 слева, 1 справа
N, Lx = 801, 8.0
xs = np.linspace(-Lx / 2, Lx / 2, N)
dx = xs[1] - xs[0]
C = (xs > 0).astype(float)
D0, beta = 0.05, 3.0
t, tend = 0.0, 1.0
while t < tend:
    D = D0 * np.exp(beta * C)
    Dh = 0.5 * (D[1:] + D[:-1])
    dt = min(0.2 * dx * dx / Dh.max(), tend - t)
    flux = -Dh * (C[1:] - C[:-1]) / dx
    C[1:-1] -= dt * (flux[1:] - flux[:-1]) / dx
    t += dt
# плоскость Матано: ∫ (x - x_M) dC = 0  ->  x_M = ∫ x dC
dC = np.diff(C)
xm = float(np.sum(0.5 * (xs[1:] + xs[:-1]) * dC) / np.sum(dC))
sel = (xs > -3) & (xs < 3)
save_dat('dif-matano.dat', [xs[sel], C[sel]], 'x C')
save('dif-matano.tex', f'\\def\\MatanoX{{{xm:.4f}}}\n')
print('Матано: x_M =', round(xm, 4))

# ---------------------------------------------------------------- модель Фишера
# C/C0 = exp(-b·y)·erfc(x/2√(D_l t)); в единицах 2√(D_l t) = 1, b = 0,35
b = 0.35
lines = []
for c in (0.6, 0.35, 0.18, 0.08, 0.03):
    pts = []
    for xv in np.linspace(0.0, 1.6, 120):
        e = math.erfc(xv)
        if e <= c:
            break
        y = math.log(e / c) / b
        pts.append((xv, y))
    if len(pts) < 2:
        continue
    pts.append((pts[-1][0] + 1e-3, 0.0))
    right = ' '.join(f'({px:.3f},{-py:.3f})' for px, py in pts)
    left = ' '.join(f'({-px:.3f},{-py:.3f})' for px, py in reversed(pts))
    lines.append(f'\\draw[fisher] plot[smooth] coordinates {{{left} {right}}};')
save('dif-fisher.tex', '\n'.join(lines) + '\n')
print('Фишер: ok')

# ---------------------------------------------------------------- случайное блуждание
random.seed(7)
p = [(0, 0)]
for _ in range(400):
    dxy = random.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
    p.append((p[-1][0] + dxy[0], p[-1][1] + dxy[1]))
path = ' -- '.join(f'({a},{c})' for a, c in p)
save('dif-walk.tex', f'\\draw[walk] {path};\n\\coordinate (walkEnd) at {p[-1]};\n')
print('блуждание: конец', p[-1], 'R =', round(math.hypot(*p[-1]), 1), '√N =', 20)
