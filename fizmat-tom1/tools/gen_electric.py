#!/usr/bin/env python3
"""Рисунки к п. 3.5 (электрические свойства), рассчитанные численно.

  figures/el-bg.dat — сопротивление по формуле Блоха–Грюнайзена
      ρ_ф(T) ∝ (T/θ)^5 ∫_0^{θ/T} x^5 / ((e^x − 1)(1 − e^−x)) dx,
      нормированное на ρ_ф(θ), и полное ρ = ρ_ост + ρ_ф (правило Маттиссена)

Запуск: python3 tools/gen_electric.py"""
import math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')


def J5(y, n=2000):
    if y <= 0:
        return 0.0
    x = np.linspace(1e-6, y, n)
    f = x**5 / ((np.exp(x) - 1) * (1 - np.exp(-x)))
    return float(np.trapezoid(f, x))


def rho(t):          # t = T/θ
    return t**5 * J5(1 / t)


r1 = rho(1.0)
ts = np.concatenate([np.linspace(0.01, 0.2, 60), np.linspace(0.21, 1.2, 60)])
rows = []
for t in ts:
    r = rho(t) / r1
    rows.append((t, r, r + 0.02, r + 0.08, t / (4 * r1)))   # чистый; с ρ_ост; высокотемпературная асимптота
with open(os.path.join(OUT, 'el-bg.dat'), 'w') as f:
    f.write('t r r1 r2 lin\n' + '\n'.join(' '.join(f'{v:.6f}' for v in row) for row in rows) + '\n')
# проверка асимптот: при t>>1 ρ ≈ t/4 (в ед. K/Aθ), при t<<1 ρ ≈ 124,4 t^5
print('J5(∞) ≈', round(J5(60), 2), '(5!·ζ(5) = 124,43)')
print('ρ(1)/ρ(2) =', round(r1 / rho(2.0), 3), ' (линейный закон дал бы 0,5)')
