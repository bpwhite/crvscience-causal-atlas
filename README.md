# Causal Atlas (`causal-atlas` v0.2.0)

![Causal-Atlas Headless Workbench](validation/headless_screenshot.png)

**Causal Atlas** is an independent, clean-room native Rust desktop GUI, terminal UI, and CLI workbench for causal biological pathway, upstream regulator, mechanistic network, and perturbation signature inference backed by a bundled **SQLite 3 WAL caching engine** and open-source biomedical data streams.

- **User Guide**: [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md)
- **Scientific Methods**: [`docs/SCIENTIFIC_METHODS.md`](docs/SCIENTIFIC_METHODS.md)
- **Headless Validation**: [`validation/README.md`](validation/README.md)

## Clean-Room & Bundled Architecture

- **100% Independent Implementation**: Built from peer-reviewed statistical and graph-theoretic publications ([Krämer et al., 2014, PMC3928520](https://pmc.ncbi.nlm.nih.gov/articles/PMC3928520/); [Benjamini & Hochberg, 1995](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x); [Lamb et al., 2006](https://doi.org/10.1126/science.1132939)) without referencing proprietary SDKs, APIs, or commercial knowledge bases.
- **Zero External Runtime Dependencies**:
  - Enforces `#![forbid(unsafe_code)]` across `ca-core`, `ca-streams`, and all `ca-app` UI/export modules (`ca-app/src/main.rs` uses `#![deny(unsafe_code)]` with a single bounded `AttachConsole` call on Windows CLI invocations).
  - Statically bundles SQLite 3 (`rusqlite` with `bundled` feature) in **WAL (`Write-Ahead Logging`)** mode — no external system `libsqlite3`, R, Python, Java, or OpenSSL required.
  - Uses pure-Rust TLS (`ureq` + `rustls`) and pure-Rust `.zip` / `.gz` archive extraction (`zip` + `flate2`) for batch GPML and GMT synchronization.

## Scientific Summary

`ca-core` implements five integrated causal inference engines (detailed mathematically in [`docs/SCIENTIFIC_METHODS.md`](docs/SCIENTIFIC_METHODS.md)):

1. **Canonical Pathway Over-Representation & MAP Simulation (`canonical.rs`)**: One-sided right-tail Fisher's exact enrichment with Benjamini-Hochberg FDR adjusted over the **full tested pathway count `m`** (not only the overlapping subset), combined with signed pathway Activation Z-score and iterative Molecule Activity Predictor (`simulate_map_in_kg`) propagation.
2. **Upstream Regulator Analysis & Mechanistic Networks (`ura.rs`)**: Krämer et al. (2014) signed Activation Z-scores and skewness-adjusted Z-scores excluding unsigned (`Ambiguous` / `expected_dir == 0`) members from Z-score numerators/denominators, total-hypothesis BH-FDR, target concordance evaluated against `z.signum()` (so downregulated targets of inhibited regulators with `z <= -2.0` score `consistent = true`), and transitive-triangle Mechanistic Networks.
3. **Multi-Hop Causal Network Analysis (`cna.rs`)**: Bounded shortest-path search (`K = 1..=3`) across signaling (`E_a`) and transcriptional (`E_t`) edges with `log2` in-degree-binned Monte Carlo permutation p-values sampled via unbiased **partial Fisher-Yates** shuffling.
4. **Downstream Effects & Regulator Effects (`dea.rs`)**: Signed Gene Ontology / Open Targets biological process and disease enrichment rendered as a hierarchical TreeMap, bridged with consistent URA targets into 3-tier Regulator Effects cascades.
5. **Rank-Based Signed Perturbation Connectivity (`perturb.rs`)**: Two-sample Kolmogorov-Smirnov signed CMap / L1000 enrichment over finite descending `log2FC` ranks with explicit `NaN` / `Inf` guards, scoring therapeutic **Reversers** (`<= -25`) vs. **Mimickers** (`>= +25`).
6. **High-Precision Tail Statistics (`stats.rs`)**: Exact `ln(k!)` hypergeometric Fisher tests and `erfc_approx` combining Abramowitz & Stegun 7.1.26 (`|x| <= 5`) with an incomplete-gamma **continued-fraction tail** (`|x| > 5`) preserving ~15-digit relative accuracy deep into extreme p-value tails.

## Controls & Keyboard Shortcuts

Text-input-guarded hotkeys are ignored while typing in the Omnibox or table filters:

| Shortcut | Scope | Action |
| :--- | :--- | :--- |
| `Ctrl+F` / `Cmd+F` | Global | Focus the **Gene / Protein Omnibox** and open autocomplete |
| `1` – `5` | Global (outside text inputs) | Switch across the 5 workbench pillars (`Cohort & Volcano`, `Pathway Atlas`, `Causal Regulators`, `Phenotype & Drug Match`, `Knowledge Databases`) |
| `F` or `0` | Global (outside text inputs) | Fit the **Pathway Atlas** camera to the active diagram (`Fit All`) |
| `R` | Global (outside text inputs) | Clear all active **MAP** node clamps and drug-rescue pins |
| `Ctrl+E` / `Cmd+E` | Global (outside text inputs) | Export vector `.svg` (active Pathway MAP, or Downstream TreeMap when on that sub-tab) |
| `Ctrl+O` / `Cmd+O` | Global (outside text inputs) | Open expression CSV/TSV dataset import dialog |
| `Ctrl+S` / `Cmd+S` | Global (outside text inputs) | Export full `CoreAnalysisReport` JSON dialog |
| `F1` | Global | Toggle the **About Causal Atlas** dialog |
| `Esc` | Global (layered) | Dismiss Omnibox dropdown → About modal → Pathway Full Screen |
| Scroll / Drag / Click | Pathway & Volcano | Cursor-anchored zoom (`15%–350%`), canvas pan, drag-to-pin nodes, click volcano points or table symbols to jump to Omnibox & Pathway Atlas |

GUI productivity features include **3 Cutoff Presets** (`Stringent`: `|log2FC| >= 1.5, FDR <= 0.01`; `Standard`: `1.0, 0.05`; `Exploratory`: `0.5, 0.10`), **sortable & filterable** Gene/URA/CNA tables with **CSV export**, subtle **subcellular compartment band fills**, a **toggleable Mini-Map**, and **edge-sign & SVG activity colorbar legends**.

## CLI Reference (`causal-atlas`)

Global flags (accepted before or after subcommands):
- `--db <PATH>`: Path to the SQLite 3 WAL cache database (default: `causal_atlas.sqlite3`)
- `--screenshot <OUT_PNG>`: Render a `1600x1000` offscreen wgpu screenshot PNG without opening a window
- `--screenshot-view <VIEW>`: Select screenshot view (`1` = Comparison & Reversers, `2` = Pathway Atlas, `3` = Causal Regulators; default: `1`)

Subcommands:
- `gui` *(default)*: Launch the native cross-platform desktop GUI (`eframe` / `egui`)
- `download-all [--force] [--pathways-only] [--local-file <PATH>]`: Batch-sync WikiPathways GPML, Reactome GMT, and regulatory streams into SQLite (or import a local `.gpml`/`.gmt`/`.tsv`/`.zip`/`.gaf.gz`)
- `tui [--input <PATH>]`: Launch the interactive SSH-compatible terminal UI (`ratatui`)
- `analyze [--input <PATH>] [--output <PATH>] [--fc-cutoff <F64>] [--fdr-cutoff <F64>]`: Run headless Core Analysis (`--output` defaults to `causal_atlas_report.json`, `--fc-cutoff` `1.0`, `--fdr-cutoff` `0.05`)
- `export-svg [--input <PATH>] [--out-dir <DIR>]`: Export `downstream_effects_treemap.svg` and `<pathway_id>_pathway_map.svg` (default `--out-dir`: `.`)

```bash
# 1. Run all unit & headless validation tests across the workspace
cargo test --workspace

# 2. Launch the native desktop GUI backed by causal_atlas.sqlite3
cargo run --release -p ca-app

# 3. Batch-download all WikiPathways GPML + Reactome GMT + regulatory streams
cargo run --release -p ca-app -- --db causal_atlas.sqlite3 download-all

# 4. Run cached CLI Core Analysis on a DESeq2 / edgeR TSV file
cargo run --release -p ca-app -- --db causal_atlas.sqlite3 analyze --input examples/tnf_huvec_deg.tsv --output report.json

# 5. Launch the interactive headless Terminal UI (SSH / workstation compatible)
cargo run --release -p ca-app -- --db causal_atlas.sqlite3 tui

# 6. Export vector SVG Canonical Pathway MAP & Downstream TreeMap diagrams
cargo run --release -p ca-app -- --db causal_atlas.sqlite3 export-svg --out-dir exports
```

## Headless Validation Metrics

Run the deterministic dry-lab validation suite:

```bash
cargo test -p ca-core --test headless_validation -- --nocapture
```

Observed results on the `HUVEC + TNF Inflammatory Cascade` (`28` measured genes, `24` active DEGs: `19` up, `5` down, `mean_sign_skew = +0.583`) and contrast cohorts (`validation/causal_atlas_validation_summary.json`):

- **T1 (DESeq2 Wald & BH-FDR Calibration)**: `bh_monotonicity_verified = true`, `wald_two_sided_p_calibration_max_abs_diff_log10 = 7.11e-15`, excluding housekeeping/unrelated controls (`GAPDH`, `ACTB`, `TFF1`, `GREB1`).
- **T2 (Krämer URA & Fisher Enrichment)**: `13` scored regulators across `18` tested regulons (`5` significant activated, `1` significant inhibited at `FDR <= 0.05`), led by `TNF` (`Z = +4.47`, `FDR = 4.11e-10`), `NFKB1` (`Z = +4.00`, `FDR = 3.08e-7`), `RELA` (`Z = +3.16`, `FDR = 5.37e-4`), and `DEXAMETHASONE` (`Z = -2.64`, `FDR = 7.62e-3`); `MCF-7 + Beta-Estradiol` contrast activates `ESR1` (`Z = +4.00`) while `TNF` stays unresolved (`Z = -1.73`).
- **T3 (Multi-Hop CNA & Rank-Based CMap)**: `25` causal network hypotheses for TNF (`TNF` `K=2`, `Z = +4.69`, `net_bias_p = 9.99e-4`); recovers transcriptionally silent `SOST` (`log2FC = -0.35`, inferred `Z = -4.12`, `net_bias_p = 9.99e-4`) in Cohort 3; scores `30` perturbation signatures (`12` reversers, `14` mimickers), led by reversers `BMS-345541` (`-92.4`), `KLF2 Lentiviral Overexpression` (`-91.7`), and `TPCA-1` (`-90.0`).

Artifacts: [`validation/causal_atlas_validation_summary.json`](validation/causal_atlas_validation_summary.json), [`validation/upstream_regulators_ranked.csv`](validation/upstream_regulators_ranked.csv), and [`validation/README.md`](validation/README.md).

## Citation & CRV Science Apps Hub

- **CRV Science Apps Hub**: <https://crvscience.com/apps>
- **Preprint (`APP-PREPRINT-04`)**: [Causal-Atlas: Tracing How Genes and Drugs Turn Biological Pathways On and Off](https://crvscience.com/apps/preprints/causal-atlas-pathway-perturbation)
- **Concept DOI**: [`10.5281/zenodo.24222376`](https://doi.org/10.5281/zenodo.24222376)

```bibtex
@article{white2026causalatlas,
  author  = {Bryan P. White},
  title   = {Causal-Atlas: Tracing How Genes and Drugs Turn Biological Pathways On and Off},
  year    = {2026},
  journal = {CRV Science Apps Preprints},
  note    = {APP-PREPRINT-04, v0.2.0},
  doi     = {10.5281/zenodo.24222376},
  url     = {https://crvscience.com/apps/preprints/causal-atlas-pathway-perturbation}
}
```

See [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.
