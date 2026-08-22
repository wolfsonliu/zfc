# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.1] - 2026-08-22

### Added
- Six new result figures (each saved as both PNG and PDF):
  - `*_gene_volcano` / `*_sgrna_volcano` — volcano plots (ZLFC vs. `-log10`
    adjusted p-value, significant hits highlighted).
  - `*_gene_ma` — MA plot (`log2` mean normalized count vs. ZLFC).
  - `*_gene_top_hits` — ranked bar plot of top depleted/enriched genes.
  - `*_sgrna_p_hist` / `*_gene_p_hist` — permutation p-value distribution
    histograms.

### Changed
- Restyled all figures with a blue theme (`#007BFF` / `#004999` family) plus a
  warm orange (`#FF8C00`) reserved for reference/key lines; de-spined axes,
  light y-grid, 150 dpi output, and consistent typography. Numeric results and
  output file names are unchanged.
- Updated `docs/data-format.md` to document the new plot outputs.
- Added PyPI project metadata to `pyproject.toml` (project URLs, classifiers,
  keywords).

### Fixed
- Constrained the p-value histograms to their actual `[0, 0.5]` support, so the
  figure no longer leaves its right half empty.

## [0.2.0] - 2026-08-21

- Initial uv/hatchling packaging, `zfc` CLI entry point, lint (ruff) and test
  (pytest) wiring. (Version bump and known-issue resolution.)