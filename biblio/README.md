# Библиографическая база: разрушение Zr-сплавов с гидридами

Постановка — `TASK.md`. Форма вывода субагентов — `AGENT_SPEC.md`.

## Файлы
- `catalog.csv` — основная таблица (колонки по TASK.md); `catalog.xlsx` — собирается `python3 biblio/build_xlsx.py`.
- `needs_fulltext.csv` — что достать из полного текста (сессия с доступом к ScienceDirect).
- `candidates_unverified.csv` — интересные работы без подтверждённого DOI (в базу не входят).
- `parts/*.csv` — сырые партии по подтемам; `python3 biblio/merge.py` сливает их в `catalog.csv`
  (дедупликация по DOI и похожему названию, стабильные id, `overlap_reorient` по каталогу пользователя).

## Метод
- OpenAlex, Crossref, doi.org, Springer, Semantic Scholar из этой среды закрыты (403 / ENOTFOUND; проверялось в начале каждого хода).
  Основной инструмент — WebSearch (лимит 200 вызовов на ход на все агенты).
- Уровни `doi_check`:
  - `PDF` — DOI напечатан в PDF издателя, лежащем в корне library (равноценно V2);
  - `CAT` — DOI из каталога пользователя по переориентации (`existing_reorientation.csv`, собран с доступом к издателям);
  - `V1` — DOI виден в выдаче поиска для этого названия (часто — в URL страницы издателя);
  - `no-doi` — книги, отчёты, диссертации с ISBN/URL/handle.
  V2 (Crossref/OpenAlex) недоступен; при открытии API записи V1 нужно перепроверить пакетно.
- DOI не собирались по шаблону (PII и т. п.); сомнительные — в `candidates_unverified.csv`.
- НИРС («НИРС Ильясов Р.А. М25-109.pdf») не использовался как источник; числа, известные только из него
  (Kubo 2013 — 710 МПа, Tondro — 30 мкм, Э635 E и σ_y), в базу не внесены, пока не найден оригинал.
- Числа в `values` — с единицами и условиями, источник в `values_source` (abstract / snippet / fulltext).
  У записей REO (`notes`: «summary_ru по названию») аннотации не сверены.

## Поисковые запросы (типовые)
«<название> doi», «<автор> <год> <журнал> doi», с `allowed_domains` (sciencedirect.com, onlinelibrary.wiley.com,
tandfonline.com, link.springer.com, researchgate.net, ascelibrary.org, elibrary.ru, cyberleninka.ru, jstage.jst.go.jp);
по подтемам — ключевые слова из TASK.md («hydride fracture toughness microcantilever», «ring compression test hydrided
cladding ductile-to-brittle», «radial hydride continuity factor», «Э110 гидрирование механические свойства» и т. п.).

## Статистика (на момент последней сборки)
Всего 309. Направление 1 — 157, 2 — 138, 3 — 42.
doi_check: V1 172 (55 %), CAT 75 (24 %), PDF 51 (16 %), no-doi 11 (3 %).
Приоритет: A 74, B 139, C 96.
Пересечение с каталогом пользователя — 76; PDF в library — 61.

| подтема | записей |
|---|---|
| 1.1 | 23 |
| 1.2 | 27 |
| 1.3 | 45 |
| 1.4 | 32 |
| 1.5 | 8 |
| 1.6 | 46 |
| 1.7 | 15 |
| 1.8 | 17 |
| 2.1 | 17 |
| 2.2 | 14 |
| 2.3 | 39 |
| 2.4 | 28 |
| 2.5 | 58 |
| 2.6 | 20 |
| 2.7 | 20 |
| 2.8 | 10 |
| 3.1 | 17 |
| 3.2 | 6 |
| 3.3 | 10 |
| 3.4 | 17 |
| 3.5 | 7 |
