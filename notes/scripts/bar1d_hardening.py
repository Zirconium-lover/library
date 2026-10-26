#!/usr/bin/env python3
"""1D-аналог бруска test/mech/nl_objectivity (mkcross.py): локальная J2-пластичность
с упрочнением, повреждение по нелокальному приращению eps_p, запаздывающая схема ccx.

sigma_y = 400 + 800*eps_p (до eps_p=0.05, далее 440), E=2e5, Rice-Tracey eps_f=1.65*0.01*exp(-0.5),
u_f=0.02, брусок 6, надрез-«шатёр»: A(x)=(1-0.1*(1-|x-3|))^2 при |x-3|<1, A_min=0.81.
L в законе = max(h, l_d); l_d = 2*ell (GRADIENT) или sqrt(pi)*ell (INTEGRAL); LOCAL: L=h.
Вязкость D_v: Эйлер, eta (по умолчанию 1e-3, время 0..1, u_end=0.3).
Печатает W_post/(G_f*A_min) и разложение рассеяния.
"""
import sys
import numpy as np

E, SY0, HM, EPH, UF, LB, UEND = 2e5, 400.0, 800.0, 0.05, 0.02, 6.0, 0.30
EPSF = 1.65 * 0.01 * np.exp(-0.5)

def sy(ep):
    return SY0 + HM * np.minimum(ep, EPH)

def ep_of(seff):   # обратная к sy на упрочняющемся участке
    return (seff - SY0) / HM

def operator(xc, h, A, ell, mode):
    n = len(xc); V = A * h
    if mode == 'INTEGRAL':
        d = np.abs(xc[:, None] - xc[None, :])
        W = np.where(d <= 2 * ell, np.exp(-(d / ell) ** 2), 0.0) * V[None, :]
        return W / W.sum(axis=1, keepdims=True)
    # GRADIENT: линейные КЭ, сосредоточенная масса с площадью, среднее узлов
    nn = n + 1
    M = np.zeros(nn); K = np.zeros((nn, nn))
    for e in range(n):
        M[e] += V[e] / 2; M[e + 1] += V[e] / 2
        K[e:e+2, e:e+2] += ell ** 2 * A[e] / h * np.array([[1, -1], [-1, 1]])
    G = np.linalg.inv(np.diag(M) + K)
    P = np.zeros((nn, n))
    for e in range(n):
        P[e, e] += V[e] / 2; P[e + 1, e] += V[e] / 2
    nod = G @ P                     # узлы <- элементы
    Q = np.zeros((n, nn))
    for e in range(n):
        Q[e, e] = Q[e, e + 1] = 0.5
    return Q @ nod

def run(mode, ell, n, nsteps=3000, eta=1e-3):
    h = LB / n; xc = (np.arange(n) + 0.5) * h
    dd = np.abs(xc - 3.0)
    A = np.where(dd < 1.0, (1 - 0.1 * (1 - dd)) ** 2, 1.0)
    if mode == 'LOCAL':
        Op = np.eye(n); L = np.full(n, h)
    else:
        Op = operator(xc, h, A, ell, mode)
        ld = 2 * ell if mode == 'GRADIENT' else np.sqrt(np.pi) * ell
        L = np.full(n, max(h, ld))
    ep = np.zeros(n); om = np.zeros(n); D = np.zeros(n); Dv = np.zeros(n)
    dep_prev = np.zeros(n)
    us = np.linspace(0, UEND, nsteps + 1); dt = 1.0 / nsteps
    F_old = 0.0; W = 0.0; curve = [(0.0, 0.0)]
    Wpl = np.zeros(n)                       # пластическая работа по элементам
    beta = dt / (eta + dt) if eta > 0 else 1.0
    for s in range(1, nsteps + 1):
        # предиктор D~ (запаздывающий темп; шаг постоянный -> r=1)
        dbar = Op @ dep_prev
        Dt = D.copy()
        for k in range(n):
            if om[k] >= 1:
                Dt[k] = min(1.0, D[k] + L[k] * dbar[k] / UF)
            elif om[k] + dbar[k] / EPSF > 1:
                frac = (om[k] + dbar[k] / EPSF - 1) * EPSF
                Dt[k] = min(1.0, L[k] * frac / UF)
        Dvt = np.minimum(1.0, Dv + beta * (Dt - Dv))
        g = np.maximum(1e-9, 1 - Dvt)
        def u_of(F):
            seff = F / (A * g)
            epn = np.maximum(ep, np.where(seff > sy(ep), np.minimum(ep_of(seff), np.inf), ep))
            return np.sum(h * (epn + seff / E)), epn
        Fcap = np.min(A * g * (SY0 + HM * EPH))
        u_t = us[s]
        ucap, _ = u_of(Fcap * (1 - 1e-12))
        if ucap <= u_t:
            F = Fcap; _, epn = u_of(F)
            k = np.argmin(A * g); seff = F / (A * g)
            rest = u_t - np.sum(h * (epn + seff / E))
            epn[k] += rest / h
        else:
            lo, hi = 0.0, Fcap
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                if u_of(mid)[0] < u_t: lo = mid
                else: hi = mid
            F = 0.5 * (lo + hi); _, epn = u_of(F)
        dep = epn - ep
        seff_mid = 0.5 * (F + F_old) / (A * g)
        Wpl += seff_mid * dep * A * h * g      # работа (1-D)*sigma_eff на dep
        ep = epn
        # фиксация: фактическое нелокальное приращение
        dbar = Op @ dep
        for k in range(n):
            if om[k] >= 1:
                D[k] = min(1.0, D[k] + L[k] * dbar[k] / UF)
            else:
                om_new = om[k] + dbar[k] / EPSF
                if om_new > 1:
                    D[k] = min(1.0, L[k] * (om_new - 1) * EPSF / UF)
                om[k] = om_new
        Dv = Dvt
        dep_prev = dep
        W += 0.5 * (F + F_old) * (us[s] - us[s - 1]); F_old = F
        curve.append((u_t, F))
        if F < 1e-6 * 400 and s > 10:
            break
    c = np.array(curve); ip = np.argmax(c[:, 1])
    Wpost = np.trapezoid(c[ip:, 1], c[ip:, 0])
    # ширина пластической зоны после пика: элементы с заметным приростом eps_p
    return dict(W=W, Wpost=Wpost, Fmax=c[ip, 1], ep=ep, D=D, xc=xc, h=h, Wpl=Wpl, A=A)

if __name__ == '__main__':
    GF, AM = 408.0 * UF / 2, 0.81
    cases = [('LOCAL', 0.5, 48), ('GRADIENT', 0.5, 48), ('GRADIENT', 0.5, 96),
             ('INTEGRAL', 0.5, 48), ('INTEGRAL', 0.5, 96),
             ('GRADIENT', 0.25, 96), ('INTEGRAL', 0.25, 96),
             ('GRADIENT', 1.0, 48), ('INTEGRAL', 1.0, 48)]
    for eta in (1e-3, 0.0):
        print(f'eta={eta}')
        for mode, ell, n in cases:
            r = run(mode, ell, n, eta=eta)
            pz = np.sum(r['ep'] > 0.0101) * r['h']
            print(f'  {mode:8s} ell={ell:4.2f} h={r["h"]:.4f}: W_post/(G_f A)={r["Wpost"]/(GF*AM):.3f}  '
                  f'Fmax={r["Fmax"]:.2f}  ширина зоны eps_p>eps_f: {pz:.3f}  max eps_p={r["ep"].max():.3f}')
