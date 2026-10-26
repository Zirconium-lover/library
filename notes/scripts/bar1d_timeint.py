#!/usr/bin/env python3
"""Ошибка интегрирования по времени на эталоне Б1 (notes/time-integration.md).

1D-брусок: локальная идеальная пластичность + повреждение по нелокальному
приращению eps_p (интегральное ядро exp(-(d/l)^2), R=2l), D += l_d*d(eps_bar)/u_f,
l_d = sqrt(pi)*l. Один ослабленный элемент (sigma_y ниже на 2 %).
Две схемы:
  BE   — неявная: D_{n+1} согласовано с приращением шага n+1 (неподвижная точка);
  LAG  — как в ccx-arch2: внутри шага D~ = D_n + (l_d/u_f)*A*deps_n*(dt_{n+1}/dt_n)
         (запаздывающий темп, idtscale=1), после сходимости D_{n+1} пересчитывается
         по фактическому приращению шага (фиксация), напряжение шага — с D~.
         Это схема типа IMPL-EX (Oliver et al. 2008) для переменной повреждения.
Печатает W/G_f для N шагов по перемещению и порядок сходимости.
"""
import numpy as np
from scipy.optimize import brentq

E, SY, UF, LB, ELL, H = 1.0e5, 100.0, 0.1, 10.0, 0.5, 0.05
NE = int(round(LB / H)); S = NE // 2
LD = np.sqrt(np.pi) * ELL
sy = np.full(NE, SY); sy[S] *= 0.98
xc = (np.arange(NE) + 0.5) * H
d = np.abs(xc[:, None] - xc[None, :])
Wt = np.where(d <= 2 * ELL, np.exp(-(d / ELL) ** 2), 0.0)
A = Wt / Wt.sum(axis=1, keepdims=True)          # eps_bar_k = sum_j A_kj eps_j
K = LD / UF

def run(N, scheme, Dend=0.15, Ds=None, dtscale=True):
    ep = np.zeros(NE); D = np.zeros(NE)
    dep_prev = np.zeros(NE)
    if Ds is None:
        Ds = np.linspace(0, Dend, N + 1)
    N = len(Ds) - 1; sig_old = 0.0; W = 0.0; sig_hist = [0.0]; dt_prev = None
    for n in range(1, N + 1):
        Dl = Ds[n]
        if scheme == 'LAG':
            r = (Ds[n] - Ds[n - 1]) / dt_prev if (dtscale and dt_prev) else 1.0
            Dt = np.minimum(D + K * (A @ dep_prev) * r, 1 - 1e-12)
            comp = np.sum(H / (E * (1 - Dt)))
            st = (Dl - H * ep.sum()) / comp
            cap = np.min(sy * (1 - Dt)); k = np.argmin(sy * (1 - Dt))
            dep = np.zeros(NE)
            if st > cap:
                sig = cap
                dep[k] = (Dl - H * ep.sum() - cap * comp) / H
            else:
                sig = st
            D = np.minimum(D + K * (A @ dep), 1.0)
        else:  # BE: только ослабленный элемент течёт; ищем x = deps_S
            def state(x):
                Dn = np.minimum(D + K * A[:, S] * x, 1 - 1e-12)
                return Dn
            comp0 = np.sum(H / (E * (1 - D)))
            st = (Dl - H * ep.sum()) / comp0
            dep = np.zeros(NE)
            if st > sy[S] * (1 - D[S]):
                def f(x):
                    Dn = state(x)
                    s = sy[S] * (1 - Dn[S])
                    return H * (ep.sum() + x) + s * np.sum(H / (E * (1 - Dn))) - Dl
                xmax = (1 - D[S]) / (K * A[S, S])
                if f(xmax * (1 - 1e-12)) < 0:
                    x = xmax; sig = 0.0
                else:
                    x = brentq(f, 0.0, xmax * (1 - 1e-12), xtol=1e-16)
                    sig = sy[S] * (1 - state(x)[S])
                dep[S] = x
                D = np.minimum(D + K * A[:, S] * x, 1.0)
            else:
                sig = st
        ep += dep; dep_prev = dep; dt_prev = Ds[n] - Ds[n - 1]
        if D[S] >= 1.0 - 1e-12:
            sig = 0.0
        W += 0.5 * (sig + sig_old) * (Ds[n] - Ds[n - 1]); sig_old = sig
        sig_hist.append(sig)
    return W / (SY * UF / 2), np.max(sig_hist)

if __name__ == '__main__':
    print(f'h={H}, l={ELL}, l_d={LD:.4f}, NE={NE}; G_f={SY*UF/2}')
    res = {}
    for sch in ('BE', 'LAG'):
        prev = None
        for N in (50, 100, 200, 400, 800, 1600):
            w, smax = run(N, sch)
            res[(sch, N)] = w
            rate = '' if prev is None else f'  ошибка/пред.={abs(w-ref)/abs(prev-ref):.3f}' if False else ''
            print(f'{sch:3s} N={N:5d}: W/G_f={w:.6f}  max sigma={smax:.3f}')
            prev = w
    ref = res[('BE', 1600)]
    print('Разность LAG - BE(N=1600):')
    for N in (50, 100, 200, 400, 800, 1600):
        print(f'  N={N:5d}: {res[("LAG", N)] - ref:+.6f}')
    # переменный шаг: чередование шагов 1 и 1/10 (имитация cutback)
    rng = np.random.default_rng(1)
    dts = np.where(rng.random(400) < 0.3, 0.1, 1.0)
    Dsv = np.concatenate([[0.0], np.cumsum(dts)]) / dts.sum() * 0.15
    for sc in (True, False):
        w, _ = run(0, 'LAG', Ds=Dsv, dtscale=sc)
        print(f'LAG, переменный шаг (30% шагов в 10 раз меньше), dtscale={sc}: W/G_f={w:.6f}, отличие от BE {w-ref:+.6f}')
