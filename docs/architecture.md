# Architecture & Module Map

When to read: navigating the repo, adding a module, or tracing data flow.

## Top-level layout

- `src/zfc/` — the importable package (all Python logic).
- `doc/` — paper assets (LaTeX, figures, bib). Not runtime; don't edit here
  when changing behavior.
- `README.md` — user documentation (algorithm, usage, examples, citations).
- `pyproject.toml` — build config (hatchling), project metadata, `zfc` entry point, ruff config.
- `LICENSE` — GPLv3 or later (per `pyproject.toml`).

## Package: `src/zfc/`

- `src/zfc/__init__.py` — public re-exports + `__version__`. Public API:
  `read_raw_count`, `write_bar_result`, `write_sg_result`, `write_g_result`,
  `zfoldchange`, `zfc_plot`.
- `src/zfc/cli.py` — the `zfc` CLI: `main()` holds the argparse definition and
  orchestration.
- `src/zfc/programio.py` — I/O only. `read_raw_count` (input TSV → canonical
  DataFrame) and the three `write_*` output writers. Holds the column lists
  that define the output contract (see `docs/data-format.md`).
- `src/zfc/statsfunc.py` — self-contained statistical helpers: normalization
  (`df_normalization` + normfactor functions), small-count, ECDF, p-value
  adjustment (BH / Bonferroni), rank aggregation (robust + mean), gini. Pure
  functions; no CLI knowledge.
- `src/zfc/zfoldchange.py` — `zfoldchange()`, the pipeline that orchestrates the
  seven algorithm steps using `statsfunc` helpers. Returns
  `(bar_df, sg_df, g_df, train_data, reg)`.
- `src/zfc/plot.py` — matplotlib figure builders + `zfc_plot()` which saves all
  figures. Sets the Agg backend at import.

## Entry points

- `zfc` — the `[project.scripts]` entry point (implemented in
  `src/zfc/cli.py::main`), run via `uv run zfc` or the installed `zfc` command.
- `zfc-count` — the `[project.scripts]` entry point (implemented in
  `src/zfc/fastqcount.py::main`); a fast, gzip-aware FASTQ counter for sgRNA /
  iBAR reads.

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
- `zfc-count` output and `read_raw_count` agree on `gene/guide/barcode`
  (+ `count`) naming.

## Where to add code

- New statistical helper → `src/zfc/statsfunc.py`.
- New pipeline step → `src/zfc/zfoldchange.py`.
- New input parsing / output format → `src/zfc/programio.py`.
- New figure → `src/zfc/plot.py` + register it in `zfc_plot()`.
- New top-level command → add it under `src/zfc/` and register it in
  `[project.scripts]` (see `src/zfc/fastqcount.py` + the `zfc-count` entry).