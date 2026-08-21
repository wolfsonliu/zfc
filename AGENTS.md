# AGENTS.md

Routing guide for AI agents (and humans) working on **ZFC**. Keep this file
short: it points to topic docs under `docs/` and reads only what the current
task needs. User-facing documentation (algorithm, usage, examples, citations)
lives in `README.md`; paper assets live in `doc/` — do not duplicate them here.

## Project Overview

ZFC is a Python 3 command-line tool that computes z-scored log2 fold changes
(ZLFC) with p-values, FDR, and rank-aggregation scores for pooled CRISPR
screening count data — with or without iBAR barcodes and/or replicates. It
reads a tab-separated count table and writes per-barcode, per-sgRNA, and
per-gene result tables plus optional plots.

## Quick Start

- `uv sync`            # create .venv; install deps + the zfc package + CLI entry points
- `uvx ruff check`     # lint (add `--fix` to auto-fix what it can)
- `uv run pytest`      # run tests
- `uv run zfc --help`  # smoke-test the installed CLI

## Topic Docs

Read a doc only when its condition matches.

- **Development workflow** — `docs/development.md`
  Read before any code change: uv setup/lint/test/run, adding deps, versioning, style.
- **Architecture & module map** — `docs/architecture.md`
  Read when navigating the repo, adding modules, or tracing data flow.
- **Algorithm & statistics** — `docs/algorithm.md`
  Read before touching fold-change/stats math, null distributions, or rank aggregation.
- **Data formats & CLI contract** — `docs/data-format.md`
  Read before touching I/O, the `zfc` CLI, or the `bin/` counting scripts.

## Hard Constraints (non-negotiable)

1. **The input count table is tab-separated with exactly the header
   `gene guide barcode ctrl exp`.**
   When: any I/O, CLI, or validation change · Why: `read_raw_count`
   (`src/zfc/programio.py`) and `src/zfc/cli.py` hard-code these names · Expires: if a
   configurable/aliased schema is introduced deliberately.

2. **Do not change result-table columns or output filenames without updating
   every consumer.**
   When: touching `programio.py` write helpers or `src/zfc/cli.py` · Why: `src/zfc/plot.py`
   and downstream scripts index these columns (`zlfc`, `RRA_Score_*`,
   `Mean_Rank_*`, ...) · Expires: if an explicit, versioned output schema is adopted.

3. **Keep matplotlib headless.** `src/zfc/plot.py` must keep `matplotlib.use('Agg')`
   and never add interactive backends or `plt.show()`.
   When: any plotting change · Why: runs on servers/CI without a display ·
   Expires: if interactive/GUI plotting becomes a deliberate feature.

4. **Version is single-sourced from `src/zfc/__init__.py::__version__`
   (read by hatchling via `pyproject.toml`).** Current: `0.1.7`.
   When: every release/version change · Why: `zfc --version` reads
   `__version__` and packaging reads the same value · Expires: if the version
   moves elsewhere, update `pyproject.toml` accordingly.

5. **Do not add a runtime dependency without updating the dependency spec.**
   Runtime deps are declared in `pyproject.toml` (`[project].dependencies`):
   numpy, scipy, pandas, matplotlib, scikit-learn (import name `sklearn`).
   When: adding any new import · Why: keeps installs reproducible · Expires:
   remains while these are the deps.

6. **Lint and tests must be clean before a change is considered done.**
   `uvx ruff check` and `uv run pytest` must pass.
   When: before committing or reporting completion · Why: defined toolchain for
   this project (uv + ruff) · Expires: if the toolchain is replaced by an agreed alternative.

7. **`--normalization` names and behavior are public contract.**
   Accepts `{total, median, upper_quartile, median_ratio, none}`; changes must
   stay backward-compatible or be versioned.
   When: touching normalization code or the CLI · Why: users rely on the names;
   note the known `upper_quartile`/`quantile` bug and stale README wording are
   documented in `docs/algorithm.md` · Expires: remains while the CLI is public.

8. **The `--outprefix` directory must already exist; the CLI does not create it.**
   When: running or scripting `zfc` · Why: `zfc --help` states this and
   `write_*` (via `to_csv`) does not make parent dirs · Expires: if the CLI
   changes to create directories.