#!/usr/bin/env python3
"""Рисунки к п. 3.7 (тепловые свойства), рассчитанные численно.

  figures/th-cv.dat — теплоёмкость по Дебаю и Эйнштейну:
      C_D/3R = 3 t^3 ∫_0^{1/t} x^4 e^x/(e^x−1)^2 dx,  t = T/θ_D;
      C_E/3R = (θ_E/T)^2 e^{θ_E/T}/(e^{θ_E/T}−1)^2 при θ_E = 0,75 θ_D
      (типичное соотношение, чтобы кривые совпадали при высоких T);
      закон T^3: C/3R = (4π^4/5) t^3
  figures/th-diatomic.dat — дисперсия двухатомной цепочки (m = 1, M = 2,5,
      жёсткость K = 1, расстояние между соседями a = 1):
      ω² = K(1/m+1/M) ∓ K·√((1/m+1/M)² − 4 sin²(ka)/(mM)), 0 ≤ ka ≤ π/2
  figures/th-dos.dat, th-thetaD.dat — фононный спектр ГЦК-решётки с
      центральными силами между ближайшими соседями (расчёт по сетке
      волновых векторов) в сравнении со спектром Дебая той же площади и
      того же низкочастотного хода; эффективная θ_D(T), при которой
      формула Дебая воспроизводит теплоёмкость этого спектра
  figures/th-kappa-diel.dat — решёточная теплопроводность диэлектрика
      λ = C υ l/3, 1/l = 1/L + 1/l_U + 1/l_т,
      l_U = l_0 (θ/T) e^{θ/3T} (процессы переброса), l_т = 0,3 мм·(30 К/T)^4
      (рэлеевское рассеяние на точечных дефектах)
      (кристаллы размером L = 10 мм и 1 мм);
      θ = 1040 K, υ = 6 км/с, n = 1,17·10^29 м^-3 (по данным для сапфира)
  figures/th-kappa-cu.dat — теплопроводность металла по правилу
      сложения тепловых сопротивлений: W = b/T + W_∞ (T/T_1)²/(1+(T/T_1)²),
      W_∞ = 1/400 м·К/Вт, T_1 = 40 К (медь), b — для двух уровней чистоты
  figures/th-morse.dat, th-morse-seg.tex — потенциал Морзе
      U = D[(1 − e^{−a(r−r_0)})² − 1] и отрезки «амплитуды колебаний»
      при нескольких энергиях; середина отрезка
      r_ср = r_0 − ln(−E/D)/(2a)

Запуск: python3 tools/gen_thermal.py"""
import math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')


def save(name, header, rows, fmt='{:.6g}'):
    with open(os.path.join(OUT, name), 'w') as f:
        f.write(header + '\n')
        for r in rows:
            f.write(' '.join(fmt.format(v) for v in r) + '\n')


def debye(t):
    """C_V/3R по Дебаю при t = T/θ_D."""
    y = min(1.0 / t, 60.0)          # выше 60 подынтегральная функция ничтожна
    x = np.linspace(1e-8, y, 4000)
    f = x**4 * np.exp(x) / np.expm1(x)**2
    return 3 * t**3 * float(np.trapezoid(f, x))


def einstein(t, ratio=0.75):
    u = ratio / t
    if u > 300:
        return 0.0
    return u * u * math.exp(u) / math.expm1(u)**2


# ---------------------------------------------------------------- C(T)
ts = np.concatenate([np.linspace(0.01, 0.3, 90), np.linspace(0.31, 1.6, 70)])
rows = [(t, debye(t), einstein(t), min(1.3, 4 * math.pi**4 / 5 * t**3)) for t in ts]
save('th-cv.dat', 't D E T3', rows, '{:.6f}')
print('C_D(1) =', round(debye(1.0), 4), ' C_D(0.1) =', round(debye(0.1), 5),
      ' T^3-закон при 0,1:', round(4 * math.pi**4 / 5 * 0.001, 5))

# ---------------------------------------------------------------- двухатомная цепочка
m, M, K = 1.0, 2.5, 1.0
s = 1 / m + 1 / M
ka = np.linspace(0, math.pi / 2, 121)
root = np.sqrt(np.maximum(s * s - 4 * np.sin(ka)**2 / (m * M), 0))
ac = np.sqrt(K * (s - root))
op = np.sqrt(K * (s + root))
save('th-diatomic.dat', 'k ac op', zip(ka, ac, op), '{:.5f}')
print('щель: ω_ак,max =', round(ac[-1], 4), '= √(2K/M) =', round(math.sqrt(2 * K / M), 4),
      '; ω_опт,min =', round(op[-1], 4), '= √(2K/m) =', round(math.sqrt(2 * K / m), 4))

# ---------------------------------------------------------------- спектр ГЦК
# ближайшие соседи ГЦК (a = 1): (±1/2, ±1/2, 0) и перестановки
nbrs = []
for i in (-0.5, 0.5):
    for j in (-0.5, 0.5):
        nbrs += [(i, j, 0), (i, 0, j), (0, i, j)]
nbrs = np.array(nbrs)
unit = nbrs / np.linalg.norm(nbrs, axis=1)[:, None]
proj = np.einsum('ni,nj->nij', unit, unit)          # d̂ d̂ᵀ
rng = np.random.default_rng(1)
nk = 400000
kv = rng.random((nk, 3)) * 4 * math.pi            # куб обратного пространства 4π/a
phase = 1 - np.cos(kv @ nbrs.T)                    # (nk, 12)
D = np.einsum('kn,nij->kij', phase, proj)          # K/M = 1
w2 = np.linalg.eigvalsh(D)
nu = np.sqrt(np.maximum(w2, 0)).ravel()
numax = nu.max()
nu /= numax
bins = np.linspace(0, 1.0, 81)
hist, edges = np.histogram(nu, bins=bins)
cent = 0.5 * (edges[1:] + edges[:-1])
g = hist / (nk * (edges[1] - edges[0]))           # нормировка: ∫g dν = 3
# низкочастотный ход берём точно, из скоростей звука по всем направлениям:
# при малых k D(k) ≈ Σ_d (k·d)²/2 · d̂d̂ᵀ, число колебаний ниже ω на атом
# N(ω) = Ω ω³ Σ_i⟨υ_i⁻³⟩/(6π²), Ω = a³/4; у Дебая N(ω) = 3(ω/ω_D)³
dirs = rng.normal(size=(20000, 3))
dirs /= np.linalg.norm(dirs, axis=1)[:, None]
Ds = np.einsum('kn,nij->kij', 0.5 * (dirs @ nbrs.T)**2, proj)
vel = np.sqrt(np.linalg.eigvalsh(Ds))
inv3 = float(np.mean(np.sum(vel**-3, axis=1)))
nuD = (18 * math.pi**2 / (0.25 * inv3))**(1 / 3) / numax
gD = np.where(cent < nuD, 9 * cent**2 / nuD**3, 0.0)
rows = [(0.0, 0.0, 0.0)] + list(zip(cent, g, gD))
save('th-dos.dat', 'nu g gD', rows, '{:.5f}')
print('ГЦК: ν_D/ν_max =', round(nuD, 3), '; ∫g =', round(float(np.sum(g) * (edges[1] - edges[0])), 3))

# θ_D(T): решаем C_D(T/θ) = C_модель(T); температуры в единицах hν_max/k
def c_model(T):
    x = nu / T
    x = x[x < 200]
    e = np.exp(x)
    return float(np.sum(x * x * e / (e - 1)**2)) / (nk * 3)   # в единицах 3R

theta0 = nuD                                     # θ_D(0) в тех же единицах
rows = []
for T in np.concatenate([np.linspace(0.012, 0.1, 23), np.linspace(0.11, 0.8, 24)]):
    c = c_model(T)
    lo, hi = 0.2, 3.0                              # θ ищем бисекцией
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if debye(T / mid) > c:
            lo = mid
        else:
            hi = mid
    rows.append((T / theta0, 0.5 * (lo + hi) / theta0))
save('th-thetaD.dat', 't th', rows, '{:.5f}')
print('θ_D(T)/θ_D(0): min =', round(min(r[1] for r in rows), 3), ' при высоких T:', round(rows[-1][1], 3))

# ---------------------------------------------------------------- λ диэлектрика
kB = 1.380649e-23
theta, v, n = 1040.0, 6.0e3, 1.17e29
Ts = np.logspace(0, math.log10(600), 160)
C = np.array([3 * n * kB * debye(T / theta) for T in Ts])     # Дж/(м³·К)
l0 = 0.67e-9
lU = l0 * (theta / Ts) * np.exp(theta / (3 * Ts))
lam = {}
for L in (1e-2, 1e-3):
    leff = 1 / (1 / L + 1 / lU + (Ts / 30)**4 / 3e-4)
    lam[L] = C * v * leff / 3
save('th-kappa-diel.dat', 'T big small', zip(Ts, lam[1e-2], lam[1e-3]), '{:.5g}')
i = int(np.argmax(lam[1e-2])); j = int(np.argmax(lam[1e-3]))
print('диэлектрик: max(L=10 мм) =', round(lam[1e-2][i]), 'при', round(Ts[i], 1),
      'К; max(L=1 мм) =', round(lam[1e-3][j]), 'при', round(Ts[j], 1), 'К; λ(300 К) =',
      round(float(np.interp(300, Ts, lam[1e-2])), 1))

# ---------------------------------------------------------------- λ металла
Winf, T1 = 1 / 400, 40.0
Ts = np.logspace(0, math.log10(400), 160)
Wph = Winf * (Ts / T1)**2 / (1 + (Ts / T1)**2)
pure = 1 / (1 / (1500 * Ts) + Wph)
dirty = 1 / (1 / (150 * Ts) + Wph)
save('th-kappa-cu.dat', 'T pure dirty', zip(Ts, pure, dirty), '{:.5g}')
i = int(np.argmax(pure)); j = int(np.argmax(dirty))
print('металл: max(чистый) =', round(pure[i]), 'при', round(Ts[i], 1), 'К; max(грязный) =',
      round(dirty[j]), 'при', round(Ts[j], 1), 'К')

# ---------------------------------------------------------------- потенциал Морзе
Dm, am, r0 = 1.0, 1.6, 1.0
r = np.linspace(0.62, 3.2, 200)
U = Dm * ((1 - np.exp(-am * (r - r0)))**2 - 1)
harm = -Dm + Dm * am**2 * (r - r0)**2              # гармоническое приближение
save('th-morse.dat', 'r U H', zip(r, U, harm), '{:.5f}')
with open(os.path.join(OUT, 'th-morse-seg.tex'), 'w') as f:
    f.write('% создано tools/gen_thermal.py: уровни энергии и середины отрезков\n')
    for k, E in enumerate((-0.85, -0.65, -0.45, -0.25), start=1):
        sq = math.sqrt(1 + E / Dm)
        xl = r0 - math.log(1 + sq) / am
        xr = r0 - math.log(1 - sq) / am
        xm = r0 - math.log(-E / Dm) / (2 * am)
        f.write(f'\\draw[fmBlue, thick] (axis cs:{xl:.4f},{E}) -- (axis cs:{xr:.4f},{E});\n')
        f.write(f'\\fill[fmRed] (axis cs:{xm:.4f},{E}) circle (1.6pt);\n')
        f.write(f'\\node[font=\\scriptsize, anchor=west] at (axis cs:{xr:.4f},{E}) {{$T_{k}$}};\n')
print('Морзе: середины при E = −0,85…−0,25:',
      [round(r0 - math.log(-E) / (2 * am), 3) for E in (-0.85, -0.65, -0.45, -0.25)])

# ---------------------------------------------------------------- C_p железа и β-латуни (схемы)
# Плавные аппроксимации справочных данных (Дж/(моль·К)), без скрытых теплот:
# α-Fe: фон 22,5 + 0,0115(T−298) + магнитный λ-пик у T_C = 1042 К;
# γ-Fe (1184–1665 К), δ-Fe (1665–1809 К), жидкость (> 1809 К).
TC = 1042.0
rows = []
for T in np.concatenate([np.linspace(298, 1030, 60), np.linspace(1032, 1042, 6)]):
    rows.append((T, 22.5 + 0.0115 * (T - 298) + 50 / (1 + (TC - T) / 60)**1.1))
for T in np.linspace(1042.5, 1184, 40):
    rows.append((T, 22.5 + 0.0115 * (T - 298) + 50 / (1 + (T - TC) / 25)))
save('th-fe-alpha.dat', 'T C', rows, '{:.4g}')
save('th-fe-gamma.dat', 'T C', [(T, 33.9 + 0.008 * (T - 1184)) for T in (1184, 1665)], '{:.4g}')
save('th-fe-delta.dat', 'T C', [(T, 41.5 + 0.008 * (T - 1665)) for T in (1665, 1809)], '{:.4g}')
save('th-fe-liq.dat', 'T C', [(1809, 46.0), (1900, 46.0)], '{:.4g}')
# β-латунь CuZn: λ-аномалия упорядочения у T_K = 741 К (на моль атомов)
TK = 741.0
rows = [(T, 24.5 + 0.007 * (T - 300) + 14 / (1 + (TK - T) / 40)) for T in
        np.concatenate([np.linspace(300, 730, 50), np.linspace(731, 741, 6)])]
rows += [(T, 24.5 + 0.007 * (T - 300) + 3 / (1 + (T - TK) / 15)) for T in np.linspace(741.5, 900, 40)]
save('th-cuzn.dat', 'T C', rows, '{:.4g}')
