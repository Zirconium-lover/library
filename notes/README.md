# Указатель: заметки и статьи библиотеки по механикам ccx-arch2

## Заметки (каталог notes/)

| файл | о чём | скрипт |
|---|---|---|
| nonlocal-benchmarks.md | класс модели ccx (B&J 2002, ур. 61–63); эталон Б1 (брусок); дискретные $\ell_d$; Peerlings 1996 (аналитика); требования к сетке; артефакты | scripts/bar1d_nonlocal.py, scripts/bar1d_hardening.py |
| time-integration.md | запаздывающая схема = IMPL-EX, $O(\Delta t)$; контроль шага (Titscher 2019); вязкость: точные решения | scripts/bar1d_timeint.py |
| cohesive-uc6.md | UC6: свойства и замечания ($g_{\min}$, схема интегрирования, оси); одиночный элемент; вязкость (Ei); B-K с касательной; $K_n$/сетка/прочность; DCB/ENF/MMB; поворот осей | scripts/cohesive_ref.py |
| bulk-damage.md | Rice–Tracey, $\eta$; эталоны одной точки; член упрочнения $Hu_f^2/(6L)$; полоса; Jirásek & Bauer; Bažant & Oh; Doghri | — |
| path-following.md | ограничения Gutiérrez/Verhoosel/Börjesson и др.; ложное рассеяние секущей формы при пластичности; эталоны со snap-back | — |
| zr-hydride-data.md | $K_{IC}$, $E$, прочность δ-гидрида; граница; локусы; масштабы; оценки $G_c$, $l_{cz}$, $K_n$; пробелы | — |

## Статьи библиотеки по темам

**Нелокальная и градиентная регуляризация**
- `peerlings1996.pdf` — неявная градиентная модель, 1D-брусок с аналитикой (эталон).
- `10.1061@ASCE0733-93991987113@101512.pdf` — Pijaudier-Cabot & Bažant 1987, нелокальное повреждение (динамический стержень).
- `10.1061@ASCE0733-93992002128@111119.pdf` — Bažant & Jirásek 2002, обзор; $\ell_d$ (ур. 63), сверхнелокальная форма (ур. 64), сетка.
- `seupel2018.pdf` — градиентная модель вязкого повреждения в Abaqus, $b_e\le L/4$.
- `larsson2021.pdf` — градиентное повреждение для вязкого разрушения (термомеханика).
- `10.1016@j.cma.2019.06.013.pdf` — вязкое повреждение гетерогенных материалов (нелокально).
- `An-XFEM-strain-gradient-damage-model-….pdf` — XFEM + градиент.

**Время, IMPL-EX** — `titscher2019.pdf`; `needleman1988.pdf` (скоростная зависимость как регуляризатор).

**Crack band, объёмное повреждение** — `10.1007@bf02486267.pdf` (Bažant & Oh 1983), `10.1016@j.compstruc.2012.06.006.pdf`
(Jirásek & Bauer 2012), `doghri1995.pdf`, `gettu1998.pdf` (масштабный эффект), `10.1002@nme.4393.pdf` (встроенная разрывность),
`Failure-zone-homogenization-….pdf`.

**Когезионные элементы** — `turon2006.pdf` (B-K), `turon2007.pdf` ($K$, $N_e$, снижение прочности), `On-the-calibration-….pdf`
(Linke & Lammering 2023), `fan2019.pdf` (искусственная податливость), `segurado2004.pdf` (3D-элемент, поворот осей, управление),
`park2009.pdf`/`park2011.pdf` (PPR, обзор законов), `chaboche2001.pdf`, `hamitouche2008.pdf`, `gao2004.pdf` (вязкая регуляризация),
`alfano2001.pdf` (интерфейсные элементы), `zhang2015.pdf` (intrinsic CZM), `raghupathy2006.pdf` (сходимость в фрагментации),
`foulk2010.pdf`, `10.4028@…kem.774.84.pdf`, `goodarzi2016.pdf`, `An-incremental-elastic-plastic-cohesive-….pdf`.

**Управление по пути** — `10.1002@cnm.649.pdf` (Gutiérrez 2004), `verhoosel2008.pdf`, `s00466-022-02175-w.pdf` (Börjesson 2022),
`lorentz2004.pdf`, `hellweg1998.pdf`, `alfano2003.pdf`, `An-arc-length-control-technique-….pdf` (фазовое поле).

**Фазовое поле (для сравнения)** — `msekh2015.pdf`, `klinsmann2015.pdf`, `A-finite-volume-…`, `A-phase-field-approach-…`,
`An-efficient-adaptive-length-scale-…`, `A-numerical-experiment-…splitting…`, `Correction-Assessment-…`.

**Прочее** — `10.1111@cgf.12351.pdf` (вырожденные элементы, corotational), `hansbo2015.pdf`, `ghosh2019.pdf`,
`zhang2019.pdf` (разрывные обогащения), `liu2018.pdf` (DG), `li1995.pdf`, `10.1016@j.commatsci.2007.07.027.pdf`,
`yazid2009.pdf` (обзор XFEM), `oden1978.pdf`, `janicke1999.pdf`, `yazdani1995.pdf`, `Convergence-Analysis-…`, `Numerical-analysis-…Darcy…`,
`An-implicit-coupled-…peridynamics…`, `An-anisotropic-…multiscale-damage…`, `A-finite-strain-visco-hyperelastic-…`.

**Предметная область** — `НИРС Ильясов Р.А. М25-109.pdf` (постановка, свойства Э635 и гидрида, геометрия гидридов);
`31280_…pdf` (МИФИ, «Физическое материаловедение», т. 1 — физика твёрдого тела, данных по Zr-сплавам не содержит).
