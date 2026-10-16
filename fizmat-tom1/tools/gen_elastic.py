#!/usr/bin/env python3
"""Рисунки к п. 3.8 (упругие свойства).

  figures/el-young.dat — модуль Юнга простых веществ (ГПа) в зависимости
      от атомного номера; округлённые справочные значения для
      поликристаллов при комнатной температуре; столбец p — период
      (строки с nan разрывают линию между периодами)
  figures/el-snoek.dat — пик Снука в тантале с кислородом при двух
      частотах крутильного маятника (0,31 и 1,2 Гц):
      Q⁻¹ = 2 Q⁻¹_max ωτ/(1 + ω²τ²), τ = τ_0 exp(H/RT),
      H = 106 кДж/моль, τ_0 = 10⁻¹⁴ с
  figures/el-spectrum.dat — схематический спектр релаксации металла при
      комнатной температуре (по К. Зинеру): сумма дебаевских пиков
      Q_i·2u/(1+u²), u = ν/ν_i

Запуск: python3 tools/gen_elastic.py"""
import math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')

# ---------------------------------------------------------------- модуль Юнга
# Z: (символ, E, ГПа)
young = {
    3: ('Li', 4.9), 4: ('Be', 287), 11: ('Na', 10), 12: ('Mg', 45), 13: ('Al', 70),
    14: ('Si', 150), 19: ('K', 3.5), 20: ('Ca', 20), 21: ('Sc', 74), 22: ('Ti', 116),
    23: ('V', 128), 24: ('Cr', 279), 25: ('Mn', 198), 26: ('Fe', 211), 27: ('Co', 209),
    28: ('Ni', 200), 29: ('Cu', 130), 30: ('Zn', 108), 32: ('Ge', 103), 37: ('Rb', 2.4),
    38: ('Sr', 16), 39: ('Y', 64), 40: ('Zr', 88), 41: ('Nb', 105), 42: ('Mo', 329),
    44: ('Ru', 447), 45: ('Rh', 380), 46: ('Pd', 121), 47: ('Ag', 83), 48: ('Cd', 50),
    49: ('In', 11), 50: ('Sn', 50), 51: ('Sb', 55), 55: ('Cs', 1.7), 56: ('Ba', 13),
    57: ('La', 37), 58: ('Ce', 34), 59: ('Pr', 37), 60: ('Nd', 41), 62: ('Sm', 50),
    63: ('Eu', 18), 64: ('Gd', 55), 65: ('Tb', 56), 66: ('Dy', 61), 67: ('Ho', 65),
    68: ('Er', 70), 69: ('Tm', 74), 70: ('Yb', 24), 71: ('Lu', 69), 72: ('Hf', 78),
    73: ('Ta', 186), 74: ('W', 411), 75: ('Re', 463), 76: ('Os', 560), 77: ('Ir', 528),
    78: ('Pt', 168), 79: ('Au', 78), 81: ('Tl', 8), 82: ('Pb', 16), 83: ('Bi', 32),
}
period = lambda z: 2 if z <= 10 else 3 if z <= 18 else 4 if z <= 36 else 5 if z <= 54 else 6
with open(os.path.join(OUT, 'el-young.dat'), 'w') as f:
    f.write('Z E p\n')
    prev = None
    for z in sorted(young):
        p = period(z)
        if prev is not None and p != prev:
            f.write(f'{z - 0.5} nan {p}\n')     # nan — разрыв линии между периодами
        f.write(f'{z} {young[z][1]} {p}\n')
        prev = p

# ---------------------------------------------------------------- пик Снука Ta–O
R, H, tau0 = 8.314, 106e3, 1e-14
T = np.linspace(330, 500, 171)
tau = tau0 * np.exp(H / (R * T))
rows = []
for nu in (0.31, 1.2):
    wt = 2 * math.pi * nu * tau
    rows.append(2 * 4.0 * wt / (1 + wt**2))      # Q⁻¹·10³, максимум 4·10⁻³
with open(os.path.join(OUT, 'el-snoek.dat'), 'w') as f:
    f.write('T q1 q2\n')
    for t, a, b in zip(T, *rows):
        f.write(f'{t:.1f} {a:.4f} {b:.4f}\n')
Tpk = [float(T[np.argmax(r)]) for r in rows]
Hcalc = R * math.log(1.2 / 0.31) / (1 / Tpk[0] - 1 / Tpk[1])
print('пики Снука при T =', [round(t, 1) for t in Tpk], 'К; H по двум частотам =',
      round(Hcalc / 1e3, 1), 'кДж/моль')

# ---------------------------------------------------------------- спектр релаксации
peaks = [(1e-6, 1.0), (1e-4, 2.2), (1e-2, 1.4), (1.0, 3.0), (1e2, 1.8), (1e5, 1.2)]
lg = np.linspace(-8, 7, 400)
nu = 10**lg
q = sum(Q * 2 * (nu / n0) / (1 + (nu / n0)**2) for n0, Q in peaks)
with open(os.path.join(OUT, 'el-spectrum.dat'), 'w') as f:
    f.write('lg q\n')
    for a, b in zip(lg, q):
        f.write(f'{a:.4f} {b:.4f}\n')
