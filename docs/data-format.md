# Data Formats & CLI Contract

When to read: before touching I/O, the `zfc` CLI, the `bin/` counting scripts,
or validation logic.

## Input count table (`zfc`)

Tab-separated file with a header row naming exactly these five columns:

```
gene	guide	barcode	ctrl	exp
```

- `gene`, `guide`, `barcode` — string identifiers.
- `ctrl`, `exp` — raw counts (numeric; may be fractional for merged replicates).
- `barcode` may equal `guide` when there is no iBAR (one barcode per guide).

Source: `src/zfc/programio.py::read_raw_count` and `src/zfc/cli.py` help. Columns are
hard-coded and renamed to these canonical names on read.

## Input scenarios

1. iBAR, no replicates — `barcode` holds the iBAR.
2. No iBAR, replicates — format replicates as barcode entries (`R_A`, `R_B`, ...).
3. iBAR + replicates — merge as `<iBAR>-<replicate>` in the `barcode` column.
4. Neither — copy `guide` into `barcode`.

(Full worked examples: `README.md`.)

## Output files (written in the `--outprefix` directory)

- `<outprefix>_barcode.txt` — columns: `ctrl, exp, fc, lfc, lfc_std, zlfc, p`.
- `<outprefix>_sgrna.txt` — columns: `zlfc, count, p, p_adj, RRA_Score_down,
  RRA_Score_down_adj, RRA_Score_up, RRA_Score_up_adj, Mean_Rank_down,
  Mean_Rank_up`; sorted ascending by `zlfc`.
- `<outprefix>_gene.txt` — same columns as the sgRNA file, at gene level.

Column lists are defined once in `src/zfc/programio.py`
(`write_bar_result`/`write_sg_result`/`write_g_result`) — that file is the
single source of truth; plots read from these columns. The `--outprefix`
directory must exist before analysis (the CLI does not create it).

## Plot outputs (with `--plot`)

`<outprefix>_lorenz`, `counts_boxplot`, `normcount_scatter`,
`lfc_normcount_scatter`, `zlfc_normcount_scatter`, `lfc_std_scatter`,
`sgrna_zlfc_hist`, `gene_zlfc_hist`, `sgrna_zlfc_rra_scatter`,
`gene_zlfc_rra_scatter`, `sgrna_zlfc_rra_fdr_scatter`,
`gene_zlfc_rra_fdr_scatter`, `gene_volcano`, `sgrna_volcano`, `gene_ma`,
`gene_top_hits`, `sgrna_p_hist`, `gene_p_hist`, `lfcstd_ctrlmean` — each
saved as `.png` and `.pdf`. See `src/zfc/plot.py::zfc_plot`.

## FASTQ counting (`zfc-count`)

The legacy bash tools `bin/library_count_sgrna` and
`bin/library_count_sgrna_with_barcode` are deprecated in favor of the Python
`zfc-count` CLI (module `src/zfc/fastqcount.py`):

```
zfc-count --mode {guide,barcode} -1 R1.fastq[.gz] -2 R2.fastq[.gz] [-l library.txt]
```

- Without `-l` → raw counts: `guide<TAB>count` (guide mode) or
  `guide<TAB>barcode<TAB>count` (barcode mode); not directly consumable by
  `zfc` (no `gene` column).
- With `-l` → `gene<TAB>guide<TAB>barcode<TAB>count`, matching the `zfc`
  input shape (counts become ctrl/exp after pivoting samples). Guide mode
  emits `barcode == guide`.
- Library (guide mode): `gene <tab> guide`; (barcode mode):
  `gene <tab> guide <tab> barcode`. A header row beginning with `gene` is
  skipped.
- gzip input is auto-detected by magic bytes (no extension required).
- Patterns default to the iBAR design; overridable via `-p`/`-a`/`-b`.
  Fuller details: `docs/specs/2026-08-24-fastqcount-design.md`.

## Normalization option names

`--normalization` accepts: `total`, `median`, `upper_quartile`, `median_ratio`,
`none`. (The former `upper_quartile`/`quantile` branch bug was fixed in
`0.2.0` — see `docs/algorithm.md`.)

## Rules

- Header/column names above are a public contract; changing them breaks
  `read_raw_count` and all downstream output/plots.
  - Why: names are hard-coded in `programio.py` and `cli.py`.
  - When: any I/O or schema change.
  - Expires: if a configurable/aliased schema is introduced.

- Result-table columns must stay stable even if internal DataFrames change.
  - Why: `plot.py` and downstream users index by name.
  - When: any output change.
  - Expires: if an explicit, versioned output schema is adopted.