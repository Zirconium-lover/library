#!/usr/bin/env python3
"""Эталоны для notes/nonlocal-benchmarks.md (1D-брусок).

A. Peerlings et al. 1996, разд. 5: аналитика градиентной модели
   «идеального повреждения» D = 1 - kappa_i/kappa (ур. 45-59).
B. Дискретная длина рассеяния l_d = 1/alpha_h(x_s,x_s) (Bazant & Jirasek 2002,
   ур. 63) для интегрального (exp(-(d/l)^2), радиус обрезки R) и градиентного
   (линейные КЭ, сосредоточенная масса, ē - l^2 ē'' = e) операторов.
C. Полуаналитика бруска «локальная идеальная пластичность + повреждение по
   нелокальному Δε_p» (класс модели ccx-arch2): F(Δ), энергия.
Запуск: python3 bar1d_nonlocal.py   (нужны numpy, scipy)
"""
import numpy as np
from scipy.optimize import brentq
from math import erf, sqrt, pi

# ---------------- A. Peerlings 1996 -----------------------------------------
def peerlings_state(beta, w, L=100.0, l=10.0, alpha=0.1, ki=1e-4, c=1.0):
    """По beta и w строит решение от торца к центру; возвращает невязку
    условия гладкости в x=l/2 и ΔL (ур. 53)."""
    sc = np.sqrt(c)
    gam = (1 - beta) / (1 - alpha) - 1
    kb, kg = np.sqrt(beta / c), np.sqrt(gam / c)
    x1, x2, x3 = l / 2, w / 2, L / 2
    A = beta * ki / np.cosh((x2 - x3) / sc)          # зона 3: (1-b)ki + A ch((x-L/2)/sc)
    s2 = A * np.sinh((x2 - x3) / sc) / sc            # наклон в x=w/2
    B1 = 0.5 * (ki + s2 / kb) * np.exp(-kb * x2)     # зона 2
    B2 = 0.5 * (ki - s2 / kb) * np.exp(kb * x2)
    v1 = B1 * np.exp(kb * x1) + B2 * np.exp(-kb * x1)
    s1 = kb * (B1 * np.exp(kb * x1) - B2 * np.exp(-kb * x1))
    res = np.arctan2(-s1, v1 * kg) - kg * x1         # tan(kg x1) = -s1/(v1 kg)
    C = np.hypot(v1, s1 / kg)
    I1 = C * np.sin(kg * x1) / kg
    I2 = (B1 * (np.exp(kb * x2) - np.exp(kb * x1)) - B2 * (np.exp(-kb * x2) - np.exp(-kb * x1))) / kb
    I3 = (1 - beta) * ki * (x3 - x2) + A * sc * (np.sinh(0.0) - np.sinh((x2 - x3) / sc))
    return res, 2 * (I1 + I2 + I3), C, gam

def peerlings_perfect(dL, L=100.0, l=10.0, alpha=0.1, ki=1e-4, c=1.0):
    """Ветвь l < w (так в статье). beta убывает с ростом dL: решение есть
    при beta in (beta_inf, beta_l), где beta_l — момент w = l, а beta_inf —
    асимптота (D(0) -> 1, нагрузка -> плато (1-beta_inf) E ki)."""
    kw = dict(L=L, l=l, alpha=alpha, ki=ki, c=c)
    def w_of(beta):
        f = lambda w: peerlings_state(beta, w, **kw)[0]
        ws = np.linspace(l * (1 + 1e-12), L * 0.999, 8000)
        fs = f(ws)
        idx = np.where(np.sign(fs[:-1]) != np.sign(fs[1:]))[0]
        return brentq(f, ws[idx[0]], ws[idx[0] + 1], xtol=1e-14) if len(idx) else np.nan
    g = lambda b: peerlings_state(b, w_of(b), **kw)[1] - dL
    bs = np.linspace(1e-3, alpha * 0.999, 2000)
    gs = np.array([g(b) for b in bs])
    ok = np.isfinite(gs)
    for i in range(len(bs) - 1):
        if ok[i] and ok[i + 1] and np.sign(gs[i]) != np.sign(gs[i + 1]):
            beta = brentq(g, bs[i], bs[i + 1], xtol=1e-15)
            w = w_of(beta)
            _, _, C, gam = peerlings_state(beta, w, **kw)
            return dict(beta=beta, gamma=gam, w=w, C=C, D0=1 - ki / C,
                        sigma_over_Eki=1 - beta)
    return None

# ---------------- B. дискретная l_d -----------------------------------------
def ld_integral(h, ell, R=2.0):
    k = np.arange(-int(R * ell / h) - 1, int(R * ell / h) + 2)
    d = np.abs(k) * h
    wts = np.where(d <= R * ell, np.exp(-(d / ell) ** 2), 0.0)
    return h * wts.sum()                     # = h / alpha_h(0)

def ld_integral_edge(h, ell, R=2.0):
    """Источник в крайнем элементе (свободный торец): вес нормируется по
    имеющимся соседям, alpha_h(0) = w(0)h / sum_{k>=0} w(kh) h."""
    k = np.arange(0, int(R * ell / h) + 2)
    d = k * h
    return h * np.where(d <= R * ell, np.exp(-(d / ell) ** 2), 0.0).sum()

def ld_gradient(h, ell, n=None, edge=False):
    """Линейные КЭ, сосредоточенная масса; источник в одном элементе,
    возврат к элементу — среднее двух узлов (как damgradient в ccx-arch2)."""
    if n is None:
        n = int(round(40 * ell / h)) | 1
    nn = n + 1
    M = np.full(nn, h); M[0] = M[-1] = h / 2
    K = np.zeros((nn, nn))
    for e in range(n):
        K[e:e+2, e:e+2] += ell ** 2 / h * np.array([[1, -1], [-1, 1]])
    A = np.diag(M) + K
    s = 0 if edge else n // 2
    f = np.zeros(nn); f[s] += h / 2; f[s + 1] += h / 2   # int N_a e dx, e=1 в элементе s
    eb = np.linalg.solve(A, f)
    dps = 0.5 * (eb[s] + eb[s + 1])
    return h / dps

# ---------------- C. брусок класса ccx-arch2 --------------------------------
def bar_ccx_class(E=1.0e5, sy=100.0, uf=0.1, Lbar=10.0, ell=0.5, kind='integral',
                  na=2001):
    """Непрерывный предел h->0, трещина в середине, вдали от торцов.
    D(x) = (a/uf)*phi(x), phi = w(x)/w(0) — нормированное ядро
    (integral: exp(-x^2/l^2); gradient: exp(-|x|/l)).
    sigma = sy*(1-a/uf);  Δ = a + sigma*int dx/(E(1-D(x))).
    Возвращает массивы a, sigma, Δ и работу W = int sigma dΔ (на ед. площади)."""
    x = np.linspace(-Lbar / 2, Lbar / 2, 20001)
    phi = np.exp(-(x / ell) ** 2) if kind == 'integral' else np.exp(-np.abs(x) / ell)
    a = np.linspace(0.0, uf * (1 - 1e-9), na)
    sig = sy * (1 - a / uf)
    comp = np.array([np.trapezoid(1.0 / (E * (1 - ai / uf * phi)), x) for ai in a])
    Delta = a + sig * comp
    W = sy ** 2 * Lbar / (2 * E) + np.trapezoid(sig, Delta)  # + упругая ветвь до пика
    return a, sig, Delta, W

if __name__ == '__main__':
    print('A. Peerlings 1996, L=100, l=10, alpha=0.1, E=2e4, ki=1e-4, c=1 mm^2')
    for dL in (0.0095, 0.01, 0.02, 0.05):
        r = peerlings_perfect(dL)
        if r is None:
            print(f'  dL={dL}: вне ветви l<w'); continue
        print(f'  dL={dL:5.3f} mm: sigma/(E ki)={r["sigma_over_Eki"]:.6f} '
              f'(F={r["sigma_over_Eki"]*2e4*1e-4*10:.5f} N), w={r["w"]:.4f} mm, '
              f'gamma={r["gamma"]:.6f}, D(0)={r["D0"]:.5f}, eps(0)={(1+r["gamma"])*r["C"]:.6f}')
    print('B. l_d/l, дискретно (непрерывно: integral sqrt(pi)*erf(R/l)=%.5f при R=2l, gradient 2)'
          % (sqrt(pi) * erf(2.0)))
    for hl in (2.0, 1.0, 0.5, 0.25, 0.125, 0.0625):
        print(f'  h/l={hl:6.4f}: integral {ld_integral(hl,1.0):.5f}  gradient {ld_gradient(hl,1.0):.5f}'
              f'   | у торца: integral {ld_integral_edge(hl,1.0):.5f}  gradient {ld_gradient(hl,1.0,edge=True):.5f}')
    print('C. брусок класса ccx-arch2 (E=1e5, sy=100, uf=0.1, Lbar=10, l=0.5)')
    for kind in ('integral', 'gradient'):
        a, s, D, W = bar_ccx_class(kind=kind)
        Gf = 100 * 0.1 / 2
        dmin = np.min(np.diff(D))
        print(f'  {kind}: W={W:.5f} (G_f={Gf}), W/G_f={W/Gf:.5f}, '
              f'Δ при пике={D[0]:.5e}, Δ при разрыве={D[-1]:.5f}, snap-back={"да" if dmin<0 else "нет"}')
        for q in (0.25, 0.5, 0.75):
            i = int(q * (len(a) - 1))
            print(f'     a/uf={a[i]/0.1:.2f}: sigma={s[i]:.3f}, Δ={D[i]:.6f}')
