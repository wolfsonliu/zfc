# Architecture & Module Map

When to read: navigating the repo, adding a module, or tracing data flow.

## Top-level layout

- `zfc/` — the importable package (all Python logic).
- `bin/` — three executable entry points (one Python CLI, two bash scripts).
- `doc/` — paper assets (LaTeX, figures, bib). Not runtime; don't edit here
  when changing behavior.
- `README.md` — user documentation (algorithm, usage, examples, citations).
- `setup.py` — legacy packaging (superseded by `pyproject.toml` in the uv migration).
- `LICENSE` — GPLv3 or later (per `setup.py` classifiers).

## Package: `zfc/`

- `zfc/__init__.py` — public re-exports + `__version__`. Public API:
  `read_raw_count`, `write_bar_result`, `write_sg_result`, `write_g_result`,
  `zfoldchange`, `zfc_plot`.
- `zfc/programio.py` — I/O only. `read_raw_count` (input TSV → canonical
  DataFrame) and the three `write_*` output writers. Holds the column lists
  that define the output contract (see `docs/data-format.md`).
- `zfc/statsfunc.py` — self-contained statistical helpers: normalization
  (`df_normalization` + normfactor functions), small-count, ECDF, p-value
  adjustment (BH / Bonferroni), rank aggregation (robust + mean), gini. Pure
  functions; no CLI knowledge.
- `zfc/zfoldchange.py` — `zfoldchange()`, the pipeline that orchestrates the
  seven algorithm steps using `statsfunc` helpers. Returns
  `(bar_df, sg_df, g_df, train_data, reg)`.
- `zfc/plot.py` — matplotlib figure builders + `zfc_plot()` which saves all
  figures. Sets the Agg backend at import.

## `bin/` entry points

- `bin/zfc` — Python CLI (argparse). Parses args → `read_raw_count` →
  `zfoldchange` → `write_*` → optional `zfc_plot`. The canonical entry point.
- `bin/library_count_sgrna` — bash; counts sgRNAs from paired-end FASTQ against
  a library, emitting `gene guide barcode count`.
- `bin/library_count_sgrna_with_barcode` — bash; same, but also extracts the
  iBAR barcode from reads.

## Data flow

```
input.tsv ──read_raw_count──▶ DataFrame[gene, guide, barcode, ctrl, exp]
   └─ zfoldchange() ─▶ (bar_df, sg_df, g_df, train_data, reg)
        ├─ write_bar_result ─▶ <outprefix>_barcode.txt
        ├─ write_sg_result   ─▶ <outprefix>_sgrna.txt
        ├─ write_g_result    ─▶ <outprefix>_gene.txt
        └─ zfc_plot()        ─▶ <outprefix>_*.png / *.pdf
```

## Coupling notes (change these together)

- Output columns are defined once in `programio.py` and consumed by
  `plot.py` (`RRA_Score_down/up`, `zlfc`, `lfc`, `lfc_std`, ...). Never change
  a column name in only one place.
- The `normalization` string flows CLI → `zfoldchange` → `df_normalization`;
  keep the three in sync (see the `upper_quartile`/`quantile` note in
  `docs/algorithm.md`).
- The bash count scripts and `read_raw_count` agree on `gene/guide/barcode`
  (+ `count`) naming.

## Where to add code

- New statistical helper → `zfc/statsfunc.py`.
- New pipeline step → `zfc/zfoldchange.py`.
- New input parsing / output format → `zfc/programio.py`.
- New figure → `zfc/plot.py` + register it in `zfc_plot()`.
- New top-level command → `bin/` (later registered as a `[project.scripts]` entry).