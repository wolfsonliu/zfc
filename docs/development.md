# Development Workflow

When to read: before making any code change. Covers how to set up, lint, test,
and run the project, plus dependency/versioning/style conventions.

## Setup

- `uv sync` creates `.venv`, installs dependencies + the `zfc` package + CLI
  entry points (`zfc`, `library_count_sgrna`, `library_count_sgrna_with_barcode`).
  - Why (source): the project is migrating from `setup.py` /
    `python3 setup.py install` to a `pyproject.toml` + uv layout.
  - When: first checkout, and after any dependency change.
  - Expires: if the packaging tool is replaced.

> Status: `pyproject.toml` is not yet present. Writing it — with `[project]`,
> `[project.scripts]` (move the three `bin/` entry points), a dev dependency
> group (pytest, ruff), and `[tool.ruff]` — is the first improvement task.
> Until then the legacy install is `python3 setup.py install`.

## Lint

- `uvx ruff check` — lint.
- `uvx ruff check --fix` — auto-fix what it can, then re-check.
  - Why: uv + ruff is the chosen lint toolchain (per project decision).
  - When: before committing, on every change.
  - Expires: if an alternative is adopted.

Ruff settings belong in the future `pyproject.toml` under `[tool.ruff]`. Keep
the tree "ruff clean"; do not add `# noqa` unless it carries a one-line reason.

## Tests

- `uv run pytest`.
  - Why: pytest is the standard runner; adopted as part of the uv migration.
  - When: after code changes; in CI once wired.
  - Expires: if a different runner is adopted.

Status: no test suite exists yet. Add tests under `tests/` that lock the I/O
contract (`docs/data-format.md`) and the algorithm (`docs/algorithm.md`) with
small fixtures. Priority coverage: the `upper_quartile` normalization branch,
output column stability, and a tiny end-to-end `zfoldchange` run.

## Run the CLI

- `uv run zfc --help` (or `uv run bin/zfc --help` before entry points exist).
- `library_count_sgrna` and `library_count_sgrna_with_barcode` are bash scripts
  consumed from the shell; see `docs/data-format.md` for their I/O.

## Adding dependencies

- `uv add <package>` (uv updates `pyproject.toml` + lockfile).
- Runtime deps: numpy, scipy, pandas, matplotlib, scikit-learn (`sklearn` is
  the import name; the legacy PyPI alias is `sklearn`, but prefer
  `scikit-learn` in `pyproject.toml`).
  - Why: the scientific stack used in `zfc/*.py`.
  - When: whenever an import is added.
  - Expires: remains until the dependency list changes by design.

## Versioning

Bump `zfc/__init__.py::__version__` and `setup.py::version` together
(current: `0.1.7`).
  - Why: `bin/zfc --version` reads `__version__`; packaging reads `setup.py`.
  - When: every release.
  - Expires: when `pyproject.toml` single-sources the version.

## Style & hygiene

- Prefer plain, dependency-light functions; give math helpers a one-line
  docstring citing their reference (existing pattern in `zfc/statsfunc.py`).
- Delete dead code rather than leaving it (see "Known issues" in
  `docs/algorithm.md`).
- Only ruff is configured for linting; do not introduce unrelated reformatting
  churn in files a change doesn't touch.

## Migration checklist (first improvement task)

1. Add `pyproject.toml` with `[project]`, `[project.scripts]`, dev group
   (pytest, ruff), and `[tool.ruff]`.
2. `uv lock` / `uv sync`; confirm `uv run zfc --help` works.
3. Add a minimal test suite; make `uv run pytest` and `uvx ruff check` green.
4. Remove `setup.py` only after the entry points and version are single-sourced.