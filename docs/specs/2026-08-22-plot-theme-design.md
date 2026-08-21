# Design: Blue Theme & Cleanup for ZFC Plot Figures

## Problem Statement

`zfc --plot` currently emits figures with raw matplotlib defaults: pure-white
background, four-sided black axes, no grid, default small typography, and a
two-color scheme of all-black data plus red dotted reference lines. The result
is visually flat and inconsistent across the 13 figure types, which hurts
readability in reports and papers.

## Proposed Solution

Introduce an internal blue theme in `src/zfc/plot.py` (no new runtime
dependencies): a module-level color palette plus a `_apply_theme()` helper that
sets `matplotlib.rcParams` once at import. Refactor each figure builder to draw
from the palette instead of hard-coded `'black'`/`'red'`, de-spine axes, add a
light y-grid, and normalize typography while keeping every figure's data
content and semantics unchanged.

Warm orange is reserved for "reference / key" lines (identity diagonal, the
zero-fc/zfc axes, and the regression fit line) so they stay visually distinct
from the blue data body.

## Decisions (locked via Q&A)

| Question | Decision |
| --- | --- |
| Restyle depth | Full restyle (theme + de-spine + grid + typography), not recolor-only |
| Reference-line color | Keep a warm accent (`#FF8C00`) |
| Implementation approach | Internal theme layer in `plot.py` (no `.mplstyle`, no seaborn) |
| CLI surface | No new `--theme`/`--style` flag; theme is fixed internally |

## Detailed Design

### Palette (module-level constants in `src/zfc/plot.py`)

| Role | Value | Use |
| --- | --- | --- |
| `BLUE` | `#007BFF` | primary data color: scatter points, histogram bodies, box fills |
| `BLUE_DEEP` | `#004999` | structural/secondary series: titles, axis labels, box outlines, ctrl curve |
| `BLUE_MID` | `#66ADFF` | second data series (exp curve in lorenz) |
| `BLUE_LIGHT` | `#E8F1FE` | light fill (histogram faces, patches) |
| `ACCENT` | `#FF8C00` | reference/key lines (identity diagonal, 0 axes, regression fit) |
| `GRID` | `#D6DEE8` | grid lines |
| `TEXT_MUTED` | `#5A6B7B` | secondary tick text |

### Theme helper

- `_apply_theme()` runs once at module import (keeps `matplotlib.use('Agg')`
  first — AGENTS.md hard-constraint #3).
- Sets: base font size to 11, `axes.grid(axis='y')` with `GRID` at a thin
  weight, outward ticks, `figure.dpi` 150 (and matching `savefig.dpi`), and a
  consistent default `axes.edgecolor`/`text.color` derived from `BLUE_DEEP`.
- `_despine(ax)` removes top/right spines and keeps left/bottom.

### Per-figure changes (13 figure types, `zfc_plot` contract unchanged)

- `lorenz_curve` — ctrl curve `BLUE_DEEP`, exp curve `BLUE_MID`; identity
  diagonal `ACCENT` (dotted); Gini-index title in `BLUE_DEEP`; grid + de-spine.
- `counts_boxplot` — `patch_artist=True` with `BLUE_LIGHT` faces and
  `BLUE_DEEP` edges; whiskers/medians in `BLUE_DEEP`; de-spine.
- `normcount_scatter` / `lfc_normcount_scatter` / `zlfc_normcount_scatter` —
  scatter `BLUE` (alpha 0.15, `edgecolor='none'`); marginal histograms
  `BLUE` faces; identity diagonal / 0-axis lines `ACCENT`; de-spine.
- `lfcstd_lfc_scatter` — scatter `BLUE`; 0-axis `ACCENT`.
- `lfcstd_ctrlmean_figure` — scatter `BLUE`; regression line `ACCENT`
  (dotted); slope/intercept title in `BLUE_DEEP`.
- `zlfc_hist` (sgrna + gene) — histogram `BLUE` faces; 0 vertical line
  `ACCENT`.
- `zlfc_rra_scatter` / `zlfc_rrafdr_scatter` (sgrna + gene) — scatter `BLUE`;
  normalized y-axis label wording (`-log10(RRA)`, `-log10(RRA FDR)`).
- `zfc_plot` — unchanged signature and call order; only palette/theme driven.

### Invariants (do not change)

- `zfc_plot` and every figure-builder signature; output filenames; result-table
  columns (AGENTS.md hard-constraint #2).
- `matplotlib.use('Agg')` and no interactive backends / `plt.show()`
  (hard-constraint #3).
- Runtime dependencies stay `{numpy, scipy, pandas, matplotlib, scikit-learn}`
  — no seaborn (hard-constraint #5).

## Implementation Order

1. Add palette constants + `_apply_theme()` + `_despine()` to `plot.py`.
2. Refactor the 13 figure builders to use the palette and helpers.
3. `uvx ruff check` and `uv run pytest` (must stay green, incl. `test_plot.py`).
4. Re-run `zfc --plot` on the three `tmp/` count files and confirm each output
   folder still has 13 PNG + 13 PDF with the new look.

## Success Criteria

- `uvx ruff check` clean; `uv run pytest` green.
- All 13 figure types render without error and share the blue/orange theme.
- Three re-runs (`A549_High`, `A549_HL`, `A549_Low`) each produce 13 PNG + 13
  PDF in their own folder, consistent across samples.

## Out of Scope

- A user-facing `--theme`/`--style` selector.
- Changing figure data content, statistics, or output columns/filenames.
- Adding seaborn or any other runtime dependency.
- Changing the `zfc_plot` public signature.

## Open Questions

- Exact orange shade — settled on `#FF8C00`; `#FFA726` is the fallback if the
  added contrast is preferred later.