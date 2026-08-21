# Design: Switch the Project to uv

## Problem Statement

The project is packaged with legacy setuptools (`setup.py` +
`python3 setup.py install`): there is no lockfile, no lint configuration, and
no test suite. Improvement development will use uv:

- `uv sync` — create `.venv`, install deps + the package + CLI entry points
- `uvx ruff check` — lint (`--fix` auto-fixes what it can)

We need a `pyproject.toml`-based uv layout so the environment, the `zfc`
package, and its CLI entry point install reproducibly, with lint and a test
runner wired in.

## Proposed Solution

Adopt uv with the `hatchling` build backend. Add a `pyproject.toml` that
declares project metadata, a `zfc` console-script entry point, and a `dev`
dependency group (`pytest`, `ruff`). Single-source the version from
`zfc/__init__.py`. Add a thin `Makefile` wrapper and a minimal pytest smoke
suite. Keep the two bash counting scripts as repo-only tools (not installed).

## Decisions (locked via Q&A)

| Question | Decision |
| --- | --- |
| Migration scope | Minimal — no bash-script porting |
| `requires-python` | `>=3.11` |
| Makefile wrapper | Yes (`make setup` / `make test` / `make check`) |
| Test suite | Minimal pytest smoke suite |
| Bash scripts on PATH | Repo-only, not installed |

## Detailed Design

### Packaging — `pyproject.toml`

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "zfc"
dynamic = ["version"]
description = "Calculate z-scored log2 fold change for pooled CRISPR screening count data."
readme = "README.md"
requires-python = ">=3.11"
license = "GPL-3.0-or-later"
dependencies = ["numpy", "scipy", "pandas", "matplotlib", "scikit-learn"]

[project.scripts]
zfc = "zfc.cli:main"

[dependency-groups]
dev = ["pytest", "ruff"]

[tool.hatch.version]
path = "zfc/__init__.py"

[tool.ruff]
target-version = "py311"
```

- The dependency is `scikit-learn` (the imported package name is `sklearn`,
  which `zfc/zfoldchange.py` imports via `from sklearn import linear_model`).
- Version is single-sourced from `zfc/__init__.py::__version__`, superseding
  `AGENTS.md` hard-constraint #4 (previously two places).

### CLI entry point — `zfc/cli.py`

Move the argparse definition and orchestration currently in `bin/zfc` into a
`zfc/cli.py::main()` function (same arguments and defaults). Register it as
`[project.scripts] zfc = "zfc.cli:main"`.

Keep `bin/zfc` as a two-line shim (`from zfc.cli import main; main()`) so
direct invocation keeps working.

### `Makefile`

```make
.PHONY: setup test lint check
setup: ; uv sync
test:  ; uv run pytest
lint:  ; uvx ruff check
check: lint test
```

### Tests — `tests/` (pytest, smoke only)

- `tests/test_cli.py` — `--help` and `--version` exit cleanly.
- `tests/test_normalization.py` — `df_normalization` accepts all five documented
  method names (`total`, `median`, `upper_quartile`, `median_ratio`, `none`)
  and returns the same shape.
- `tests/test_io.py` — `read_raw_count` yields the canonical columns/index;
  `write_*` produce files with the expected columns.
- `tests/test_zfoldchange.py` — end-to-end `zfoldchange()` on a synthetic
  fixture (~1000 barcodes) asserts the returned tuple's shapes/columns.

Fixture note: `zfoldchange` bins `ctrl` into 200 bins for the LFC-std linear
model, so a tiny fixture can yield an empty training set and crash; the
smoke fixture must be large enough to populate bins.

### Housekeeping

- Remove `setup.py` (replaced by `pyproject.toml`).
- Update `.gitignore` (add `.venv/`, `.pytest_cache/`, `.ruff_cache/`).
- Commit `uv.lock`.
- Update `README.md` Installation (clone → `uv sync`); keep the PyPI section.
- Update `AGENTS.md` + `docs/*.md` to match the post-migration reality
  (single-sourced version, drop the legacy install note, mark the migration
  checklist done, note the CLI lives in `zfc/cli.py`).

## Implementation Order

1. Add `pyproject.toml` and `zfc/cli.py`; remove `setup.py`.
2. `uv sync` (generates and locks `uv.lock`); verify `uv run zfc --help`.
3. Add `Makefile`, `tests/`, ruff config, and `.gitignore` updates.
4. `uvx ruff check --fix` then `uv run pytest` until green.
5. Update `README.md`, `AGENTS.md`, and `docs/*.md` to match.

## Success Criteria

- `uv sync` succeeds and installs the `zfc` entry point.
- `uv run zfc --help` and `uv run zfc --version` work.
- `bin/zfc` shim still works.
- `uvx ruff check` is clean.
- `uv run pytest` is green.
- `make setup`, `make test`, `make check` all work.
- `uv.lock` is committed and `setup.py` is removed.

## Out of Scope

- Porting the two bash scripts to Python entry points.
- Fixing the `upper_quartile` normalization bug (tracked in `docs/algorithm.md`).
- Deleting the dead `np.math.factorial` code in `zfc/statsfunc.py`
  (not imported, so not a blocker under NumPy 2).
- CI/CD and PyPI-publish configuration.
- Pinning exact dependency versions (loose floor only for now).

## Open Questions

- Exact ruff rule selection / line length (start with defaults +
  `target-version = "py311"`).
- `license` field form (SPDX string `GPL-3.0-or-later`; fall back to
  `{ text = ... }` if the toolchain complains).