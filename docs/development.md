# Development Workflow

When to read: before making any code change. Covers how to set up, lint, test,
and run the project, plus dependency/versioning/style conventions.

## Setup

- `uv sync` creates `.venv`, installs dependencies + the `zfc` package + the
  `zfc` CLI entry point (declared in `pyproject.toml`).
  - Why (source): packaging is `pyproject.toml` (hatchling) + uv.
  - When: first checkout, and after any dependency change.
  - Expires: if the packaging tool is replaced.

The two bash counting scripts (`bin/library_count_sgrna`,
`bin/library_count_sgrna_with_barcode`) are repo-only tools; run them as
`./bin/...` — they are not installed onto PATH.

## Lint

- `uvx ruff check` — lint.
- `uvx ruff check --fix` — auto-fix what it can, then re-check.
  - Why: uv + ruff is the chosen lint toolchain (per project decision).
  - When: before committing, on every change.
  - Expires: if an alternative is adopted.

Ruff settings live in `pyproject.toml` (`[tool.ruff]`, rule set pinned to
`E4, E7, E9, F`). Keep the tree "ruff clean"; do not add `# noqa` unless it
carries a one-line reason.

## Tests

- `uv run pytest`.
  - Why: pytest is the standard runner; adopted as part of the uv migration.
  - When: after code changes; in CI once wired.
  - Expires: if a different runner is adopted.

A suite exists under `tests/` (10 tests): CLI smoke, normalization, I/O
contract, and a small end-to-end `zfoldchange` run. Add small-fixture tests
alongside new or changed behavior; the contracts to lock are in
`docs/data-format.md` and `docs/algorithm.md`.

## Run the CLI

- `uv run zfc --help` (installed entry point).
- `library_count_sgrna` and `library_count_sgrna_with_barcode` are bash scripts
  run from the repo as `./bin/...`; see `docs/data-format.md` for their I/O.

## Adding dependencies

- `uv add <package>` (uv updates `pyproject.toml` + lockfile).
- Runtime deps: numpy, scipy, pandas, matplotlib, scikit-learn (`sklearn` is
  the import name; the legacy PyPI alias is `sklearn`, but prefer
  `scikit-learn` in `pyproject.toml`).
  - Why: the scientific stack used in `src/zfc/*.py`.
  - When: whenever an import is added.
  - Expires: remains until the dependency list changes by design.

## Versioning

Bump `src/zfc/__init__.py::__version__` (current: `0.1.7`); hatchling reads it for
packaging and `zfc --version` reads it at runtime.
  - Why: version is single-sourced from one place.
  - When: every release.
  - Expires: if the version moves elsewhere, update `pyproject.toml` accordingly.

## Style & hygiene

- Prefer plain, dependency-light functions; give math helpers a one-line
  docstring citing their reference (existing pattern in `src/zfc/statsfunc.py`).
- Delete dead code rather than leaving it (see "Known issues" in
  `docs/algorithm.md`).
- Only ruff is configured for linting; do not introduce unrelated reformatting
  churn in files a change doesn't touch.
