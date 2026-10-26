#!/usr/bin/env python3
"""Эталоны для notes/cohesive-uc6.md.

1. uc6_law      — закон UC6 как в ccx-arch2 (cohesive_uc6.f), для сверки по точкам.
2. bk_law       — смешанная мода B-K (Turon et al. 2006, ур. 22-53) с модо-зависимой
                  Ks = Kn*(GIc/GIIc)*(Ts0/Tn0)^2 (условие Turon et al. 2010/2018),
                  история h = max lambda/lambda0(B); касательная + проверка разностями.
3. viscous_ramp — точное решение ОДУ Дюво-Лионса dDv/dt=(Db-Dv)/mu при постоянной
                  скорости раскрытия (интегральная показательная функция) против неявного Эйлера.
4. bar_snapback — стержень + когезионный слой: условие snap-back.
5. energies     — энергии одиночного элемента, в том числе остаток gmin.
"""
import numpy as np
from scipy.special import expi
from scipy.integrate import quad

# ---------------- 1. UC6 (как в коде) ---------------------------------------
def uc6_law(dl, dmax0, dv0, Kn, Tn0, Ts0, Gc, gmin=0.0, mu=0.0, dt=1.0):
    d0 = Tn0 / Kn; df = 2 * Gc / Tn0; beta = (Ts0 / Tn0) ** 2
    deff = np.sqrt(max(dl[0], 0) ** 2 + beta * (dl[1] ** 2 + dl[2] ** 2))
    dmax = max(dmax0, deff)
    db = 0.0 if dmax <= d0 else (1.0 if dmax >= df else df * (dmax - d0) / (dmax * (df - d0)))
    a = dt / (mu + dt) if mu > 0 else 1.0
    dv = min(1.0, max(dv0, dv0 + a * (db - dv0)))
    g = max(gmin, 1 - dv)
    T = np.array([Kn * dl[0] if dl[0] < 0 else g * Kn * dl[0], g * Kn * beta * dl[1], g * Kn * beta * dl[2]])
    return T, dmax, dv

# ---------------- 2. B-K (Turon) --------------------------------------------
class BK:
    def __init__(s, Kn, Tn0, Ts0, GIc, GIIc, eta):
        s.Kn, s.Tn0, s.Ts0, s.GIc, s.GIIc, s.eta = Kn, Tn0, Ts0, GIc, GIIc, eta
        s.Ks = Kn * (GIc / GIIc) * (Ts0 / Tn0) ** 2      # Turon 2010 / 2018
        s.kap = s.Ks / Kn
        s.l0I = Tn0 / Kn                                  # lambda0 в моде I
        s.m = 2 * GIc / (Kn * s.l0I ** 2)                 # lambda_f/lambda0, не зависит от B
        assert s.m > 1, 'нужно 2 GIc Kn / Tn0^2 > 1'
    def lam(s, d):
        return np.sqrt(max(d[0], 0) ** 2 + s.kap * (d[1] ** 2 + d[2] ** 2))
    def B(s, d):
        l = s.lam(d)
        return 0.0 if l == 0 else s.kap * (d[1] ** 2 + d[2] ** 2) / l ** 2
    def lam0(s, B):   # = l0I*sqrt(Gc(B)/GIc)
        return s.l0I * np.sqrt(1 + (s.GIIc / s.GIc - 1) * B ** s.eta)
    def dmg(s, h):
        if h <= 1: return 0.0
        if h >= s.m: return 1.0
        return s.m * (h - 1) / (h * (s.m - 1))
    def traction(s, d, h0):
        """d — скачок (n, s1, s2); h0 — история. Возвращает T, h, касательную."""
        Kd = np.array([s.Kn, s.Ks, s.Ks])
        l = s.lam(d); B = s.B(d)
        htr = l / s.lam0(B) if l > 0 else 0.0
        h = max(h0, htr)
        D = s.dmg(h)
        T = (1 - D) * Kd * d
        if d[0] < 0:
            T[0] = s.Kn * d[0]                        # контакт: полная Kn (Turon, ур. 18)
        C = np.diag((1 - D) * Kd)
        if d[0] < 0: C[0, 0] = s.Kn
        if htr >= h0 and 1 < h < s.m and l > 0:     # нагружение
            dDdh = s.m / ((s.m - 1) * h ** 2)
            dn = max(d[0], 0.0); ssq = d[1] ** 2 + d[2] ** 2
            dl_dd = np.array([dn / l, s.kap * d[1] / l, s.kap * d[2] / l])
            dB_dd = np.array([-2 * s.kap * ssq * dn / l ** 4,
                              2 * s.kap * d[1] * dn ** 2 / l ** 4,
                              2 * s.kap * d[2] * dn ** 2 / l ** 4])
            l0 = s.lam0(B)
            dl0_dB = 0.0 if B == 0 and s.eta < 1 else \
                s.l0I ** 2 * (s.GIIc / s.GIc - 1) * s.eta * (B ** (s.eta - 1) if B > 0 else (1.0 if s.eta == 1 else 0.0)) / (2 * l0)
            dh_dd = dl_dd / l0 - l / l0 ** 2 * dl0_dB * dB_dd
            dD = dDdh * dh_dd
            Tu = Kd * d
            if d[0] < 0: Tu[0] = 0.0
            C -= np.outer(Tu, dD)
        return T, h, C

def check_bk_tangent():
    m = BK(Kn=1e5, Tn0=50.0, Ts0=80.0, GIc=0.3, GIIc=1.0, eta=2.0)
    rng = np.random.default_rng(0); worst = 0.0
    for _ in range(200):
        d = rng.normal(size=3) * m.l0I * rng.uniform(1.2, 0.9 * m.m)
        d[0] = abs(d[0]) + 1e-3 * m.l0I
        T, h, C = m.traction(d, 0.0)
        if not (1 < h < m.m): continue
        Cfd = np.zeros((3, 3)); e = 1e-7 * np.linalg.norm(d)
        for j in range(3):
            dp = d.copy(); dp[j] += e; dm = d.copy(); dm[j] -= e
            Cfd[:, j] = (m.traction(dp, 0.0)[0] - m.traction(dm, 0.0)[0]) / (2 * e)
        worst = max(worst, np.abs(C - Cfd).max() / np.abs(Cfd).max())
    return worst, m

def bk_pure_mode_energy(m):
    """Работа до разрушения при пропорциональном нагружении для мод I, II и B=0.5."""
    out = {}
    for name, dirv in (('I', np.array([1.0, 0, 0])), ('II', np.array([0, 1.0, 0])),
                       ('B=0.5', np.array([1.0, 1 / np.sqrt(m.kap), 0]))):
        dirv = dirv / np.linalg.norm(dirv)
        tmax = 1.2 * m.m * m.lam0(1.0) / m.lam(dirv)
        ts = np.linspace(0, tmax, 200001); h = 0.0; W = 0.0; Tprev = np.zeros(3); dprev = np.zeros(3)
        for t in ts[1:]:
            d = t * dirv; T, h, _ = m.traction(d, h)
            W += 0.5 * (T + Tprev) @ (d - dprev); Tprev, dprev = T, d
        B = m.B(dirv)
        out[name] = (W, m.GIc + (m.GIIc - m.GIc) * B ** m.eta, B)
    return out

# ---------------- 3. вязкость: точное решение и неявный Эйлер ----------------
def viscous_ramp(Kn, Tn0, Gc, mu, v, t_end, nsteps):
    d0 = Tn0 / Kn; df = 2 * Gc / Tn0; c = df / (df - d0)
    t0, tf = d0 / v, df / v
    def Db(t):
        x = v * t
        return 0.0 if x <= d0 else (1.0 if x >= df else c * (1 - d0 / x))
    def Dv_exact(t):
        if t <= t0: return 0.0
        tt = min(t, tf)
        # (1/mu) int_{t0}^{tt} e^{-(tt-s)/mu} c (1 - t0/s) ds
        if tt / mu < 600:
            I = np.exp(-tt / mu) * (expi(tt / mu) - expi(t0 / mu))
        else:   # без переполнения: int_{t0}^{tt} e^{-(tt-s)/mu}/s ds
            a = max(t0, tt - 60 * mu)
            I = quad(lambda s_: np.exp(-(tt - s_) / mu) / s_, a, tt, limit=200)[0]
        val = c * (1 - np.exp(-(tt - t0) / mu)) - c * (t0 / mu) * I
        if t > tf:
            val = 1 - (1 - val) * np.exp(-(t - tf) / mu)
        return val
    ts = np.linspace(0, t_end, nsteps + 1); dt = ts[1] - ts[0]
    Dv = 0.0; err = 0.0; W_be = 0.0; Tp = 0.0
    for t in ts[1:]:
        Dv = Dv + dt / (mu + dt) * (Db(t) - Dv)
        err = max(err, abs(Dv - Dv_exact(t)))
        T = (1 - Dv) * Kn * v * t; W_be += 0.5 * (T + Tp) * v * dt; Tp = T
    W_ex = quad(lambda t: (1 - Dv_exact(t)) * Kn * v * t * v, 0, t_end, points=[t0, tf], limit=400)[0]
    Tmax_ex = max((1 - Dv_exact(t)) * Kn * v * t for t in np.linspace(0, t_end, 20001))
    return dict(err_max=err, W_exact=W_ex, W_BE=W_be, Gc=Gc, Tmax_exact=Tmax_ex, Dv=Dv_exact)

if __name__ == '__main__':
    w, m = check_bk_tangent()
    print(f'B-K: Ks/Kn={m.kap:.4f}, m=lf/l0={m.m:.2f}; макс. отн. ошибка касательной против разностей: {w:.2e}')
    for k, (W, G, B) in bk_pure_mode_energy(m).items():
        print(f'  мода {k:6s}: B={B:.3f}  W={W:.6f}  Gc(B) по B-K={G:.6f}')
    print('Вязкость (Дюво-Лионс), Kn=1e5, Tn0=50, Gc=0.3 (d0=5e-4, df=1.2e-2), v=1:')
    for r_ in (0.01, 0.03, 0.1, 0.3, 1.0):
        mu = r_ * 5e-4
        for n in (200, 2000):
            r = viscous_ramp(1e5, 50.0, 0.3, mu, 1.0, 1.2e-2 * 1.5 + 10 * mu, n)
            print(f'  mu*v/d0={r_:5.2f}, N={n:5d}: max|Dv-Dv_точн|={r["err_max"]:.2e}, '
                  f'W_точн/Gc={r["W_exact"]/0.3:.4f}, W_Эйлер/Gc={r["W_BE"]/0.3:.4f}, Tmax/Tn0={r["Tmax_exact"]/50:.4f}')
