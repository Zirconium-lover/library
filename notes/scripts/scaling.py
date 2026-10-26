#!/usr/bin/env python3
"""Безразмерные группы модели ccx-arch2 для двух масштабов (notes/scaling-levers.md).

Единицы: мм, МПа, Н/мм (= кДж/м^2), псевдовремя шага = 1.
Параметры «disk» — из test/s3rad/m12_s3rad_gc24_w.inp (ccx-arch2) и измерений его сетки;
«cell» — из fractrure-tests/gen/mkdeck.py по умолчанию (Kn 1e8, gmin 2e-8, деформация 0.25),
толщины гидридов из geom/fn50/stats.json; h_m ячейки — оценка по числу тетраэдров.
Запуск: python3 scaling.py  (или импорт и group(dict)).
"""
import math

DISK = dict(name='disk 4 mm', D=4.0, grip=1.0, E_m=108800., s_m=600., uf_m=0.015,
            E_h=130000., s_h=700., uf_h=0.004, Kn=2e5, Tn0=300., Ts0=260., Gc=2.4,
            gmin=1e-5, mu=1e-4, t_h=0.090, a_h=3.0, h_m=0.133, h_h=0.056, h_f=0.036,
            c_damp=2e-4, ell=None)
CELL = dict(DISK, name='cell 150 um', D=0.150, grip=0.0375, Kn=1e8, gmin=2e-8,
            t_h=0.82e-3, a_h=0.10, h_m=4e-3, h_h=0.8e-3, h_f=0.6e-3)

def group(p):
    g = {}
    E_m, E_h = p['E_m'], p['E_h']
    Gf_m = p['s_m'] * p['uf_m'] / 2; Gf_h = p['s_h'] * p['uf_h'] / 2
    lch_m = E_m * Gf_m / p['s_m'] ** 2; lch_h = E_h * Gf_h / p['s_h'] ** 2
    lch_i = E_h * p['Gc'] / p['Tn0'] ** 2
    d0 = p['Tn0'] / p['Kn']; df = 2 * p['Gc'] / p['Tn0']
    v = p['grip']                   # скорость захвата за единицу псевдовремени
    dmax = p['grip']                # грубая верхняя оценка раскрытия
    g['G_f матрицы, Н/мм'] = Gf_m
    g['l_ch матрицы = E G_f/s^2, мм'] = lch_m
    g['D/l_ch,m (хрупкость Bažant)'] = p['D'] / lch_m
    g['h_m/l_ch,m'] = p['h_m'] / lch_m
    g['G_f гидрида, Н/мм'] = Gf_h
    g['l_ch гидрида, мм'] = lch_h
    g['t_h/l_ch,h'] = p['t_h'] / lch_h
    g['l_ch интерфейса = E_h Gc/Tn0^2, мм'] = lch_i
    g['a_h/l_ch,i (длина гидрида)'] = p['a_h'] / lch_i
    g['N_e = l_ch,i/h_f (Turon 2007, M=1)'] = lch_i / p['h_f']
    g['alpha = Kn t_h/E_h (Turon: >=50)'] = p['Kn'] * p['t_h'] / E_h
    g['Kn h_m/E_m (Fan&Tadmor)'] = p['Kn'] * p['h_m'] / E_m
    g['delta0, мм'] = d0
    g['delta_f, мм'] = df
    g['delta_f/delta0'] = df / d0
    g['gmin Kn dmax/Tn0'] = p['gmin'] * p['Kn'] * dmax / p['Tn0']
    g['gmin Kn dmax^2/(2Gc)'] = p['gmin'] * p['Kn'] * dmax ** 2 / (2 * p['Gc'])
    g['mu v/delta0 (UC6, без локального усиления)'] = p['mu'] * v / d0
    g['D s_m/(E_m uf_m) (snap-back конструкции, >1 опасно)'] = p['D'] * p['s_m'] / (E_m * p['uf_m'])
    g['c h_m u/s_m (демпфер/внутр. сила)'] = p['c_damp'] * p['h_m'] * v / p['s_m']
    if p.get('ell'):
        g['ell/h_m'] = p['ell'] / p['h_m']; g['ell/t_h'] = p['ell'] / p['t_h']
    return g

if __name__ == '__main__':
    a, b = group(DISK), group(CELL)
    print(f'{"группа":55s} {DISK["name"]:>14s} {CELL["name"]:>14s} {"cell/disk":>10s}')
    for k in a:
        r = b[k] / a[k] if a[k] else float('nan')
        print(f'{k:55s} {a[k]:14.4g} {b[k]:14.4g} {r:10.3g}')
