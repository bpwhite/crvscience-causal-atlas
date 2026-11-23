# Causal Atlas User Guide

This guide covers interactive desktop usage, headless terminal/CLI workflows, statistical presets, pathway modeling, table sorting/filtering, and CSV/SVG exports in **Causal Atlas** (`causal-atlas`).

---

## 1. Execution Modes (GUI, TUI, and CLI)

All modes share the same persistent **SQLite 3 WAL** database (`--db`, defaulting to `causal_atlas.sqlite3`), pre-seeded offline with curated pathways, signed causal interactions, and SigCom LINCS consensus signatures.

### Native Desktop GUI
```bash
cargo run --release -p ca-app -- --db causal_atlas.sqlite3
```
Supports **Academic Dark** and **Publication Light** themes (toggled from the top toolbar), persistent window state, and background multi-stream HTTP synchronization.

### Headless Terminal UI (SSH / Workstation)
```bash
cargo run --release -p ca-app -- --db causal_atlas.sqlite3 tui --input examples/tnf_huvec_deg.tsv
```
- Switch views with `1`–`4`, `Tab`, `Left`, or `Right` (`1: Canonical Pathways`, `2: Upstream Regulators`, `3: Multi-Hop Causal Networks`, `4: Perturbation Reversers`).
- Exit cleanly with `q` or `Esc`.

### Headless CLI Analysis & Exports
- **Batch Knowledge Sync**:
  ```bash
  causal-atlas --db causal_atlas.sqlite3 download-all
  # Optional flags: --force (ignore ETag/TTL), --pathways-only, --local-file <archive_or_table>
  ```
- **Core Analysis JSON**:
  ```bash
  causal-atlas --db causal_atlas.sqlite3 analyze --input examples/tnf_huvec_deg.tsv --output report.json --fc-cutoff 1.0 --fdr-cutoff 0.05
  ```
- **Vector SVG Export**:
  ```bash
  causal-atlas --db causal_atlas.sqlite3 export-svg --input examples/tnf_huvec_deg.tsv --out-dir exports
  ```
  Writes `downstream_effects_treemap.svg` and `<top_pathway_id>_pathway_map.svg`.
- **Offscreen GPU Screenshot**:
  ```bash
  causal-atlas --screenshot out.png --screenshot-view 2
  ```
  Captures a `1600x1000` PNG (`--screenshot-view`: `1` = Comparison & Reversers, `2` = Pathway Atlas, `3` = Causal Regulators).

---

## 2. Left Control Panel & Cutoff Presets

Use the left panel to select among built-in benchmark cohorts (`HUVEC + TNF Inflammatory Cascade`, `MCF-7 + Beta-Estradiol Signaling`, `Spinal Osteoproliferation (Multi-Hop Wnt/SOST)`), import a custom DESeq2/edgeR/limma/Scanpy `.csv` / `.tsv` file, or tune statistical cutoffs:

| Preset | `|log2FC|` | `Max FDR` | Intended Use |
| :--- | :---: | :---: | :--- |
| **Stringent** | `>= 1.50` | `<= 0.01` | High-specificity core transcriptional drivers with minimal false positives |
| **Standard** | `>= 1.00` | `<= 0.05` | Balanced 2-fold expression and 5% Benjamini-Hochberg FDR default |
| **Exploratory** | `>= 0.50` | `<= 0.10` | Broader sensitivity for subtle or early signaling cascades |

You can also fine-tune `|log2FC|` (`0.25–3.0`), `Max FDR` (`0.0001–0.25`, logarithmic scale), and `CNA Hop K` (`1–3`) directly via sliders. Every parameter combination is SHA-256 cached in SQLite `analysis_result_cache` for sub-millisecond cohort switching.

---

## 3. Global Omnibox (`Ctrl+F`) & Keyboard Shortcuts

- **Omnibox (`Ctrl+F` / `Cmd+F`)**: Type any HGNC gene symbol, receptor, transcription factor, or drug/ligand. Pressing `Enter` or clicking an autocomplete row jumps straight to that molecule in **Pathway Atlas** (with a pulsing gold selection halo and smooth camera pan) or highlights it in **Cohort & Volcano** if unmapped to a pathway.
- **Context Strip**: Displays the currently focused molecule's molecular class, subcellular compartment, observed `log2FC` / FDR, and clickable pathway badges.
- **Global Hotkeys** (guarded so typing inside any text box never triggers them):
  - `1`–`5`: Switch across the 5 workbench pillars
  - `F` or `0`: Fit the Pathway Atlas camera (`Fit All`)
  - `R`: Clear all active MAP clamps (`Reset MAP Clamps`)
  - `Ctrl+E` / `Cmd+E`: Export vector `.svg` (active Pathway MAP, or Downstream TreeMap if currently viewing the TreeMap sub-tab)
  - `Ctrl+O`: Import expression CSV/TSV
  - `Ctrl+S`: Export Core Analysis JSON report
  - `F1`: Toggle About dialog
  - `Esc`: Layered dismissal (Omnibox dropdown → About window → Pathway Full Screen)

---

## 4. The Five Workbench Pillars

### Pillar 1: Cohort & Volcano (`1`)
- **Interactive Volcano Plot**: Displays `log2FC` (X-axis) against `-log10(BH-FDR)` (Y-axis) with numeric tick marks, axis titles, labeled cutoff guides (`±|log2FC|` and `FDR <= ...`), and collision-avoidant DEG labels prioritized by significance.
  - **Hover** over any point to inspect exact `log2FC` and `BH-FDR`.
  - **Click** any point to jump directly to that gene in the Omnibox and Pathway Atlas.
- **Measured Gene Table**: Filter by symbol substring and direction (`All`, `DEGs`, `Up (+)`, `Down (-)`), click any column header (`Symbol`, `log2FC`, `P-value`, `BH-FDR`) to sort ascending/descending (`▲`/`▼`), click any gene symbol to focus it, or click **Export CSV** to save the filtered & sorted table.

### Pillar 2: Pathway Atlas (`2`)
- **Ranked Pathway Drawer**: Filter canonical pathways by name or ID; bars encode `-log10(FDR)` and Activation Z polarity.
- **Subcellular Pathway Canvas**: Renders **5 shaded subcellular compartments** (`Extracellular Space`, `Plasma Membrane`, `Cytoplasm`, `Mitochondria / Organelle`, `Nucleus` with DNA double-helix promoter track), distinct entity glyphs (`Kinase`, `Phosphatase`, `Enzyme`, `Receptor`, `TF`, `Cytokine`, `Chemical/Drug`), `(P)` phosphorylation badges, and obstacle-deflected cubic Bezier causal edges.
- **Camera & Layout Controls**: Scroll to zoom (`15%–350%`) anchored at the cursor, drag background to pan, drag any node to pin it (`Unpin (N)` clears pins), click **1-Hop Ego Filter** to isolate the selected node's local neighborhood, click **Re-Layout** to re-solve the deterministic 120-iteration compartment force layout, or toggle **Hide Mini-Map** / **Show Mini-Map** (toolbar button or `[×]` inside the mini-map header).
- **Edge-Sign Legend**: Displayed at bottom-left (`Activating`, `Inhibiting`, `Binding`, `Inconsistent`).
- **MAP Clamping & Drug Rescue**: Click any node and clamp it to **Activate (+1)**, **Inhibit (-1)**, or **Unclamp**; press `R` to reset all clamps. Click a **Drug Rescue** pill (or **Rescue** in Pillar 4) to apply a compound's target and reversed-DEG clamps simultaneously across the pathway.

### Pillar 3: Causal Regulators (`3`)
- **Upstream Regulators (URA) Sub-Tab**:
  - Filter regulators by symbol, name, or type and by state (`All`, `Activated`, `Inhibited`); sort by any column (`Regulator`, `Type`, `State`, `Act Z`, `Bias Z`, `Overlap P`, `Targets`); click **Export CSV** to save the filtered/sorted URA table.
  - Inspect **Target Evidence** below the table (showing `Consistent (+)` vs. `Opposite (-)` evaluated relative to the regulator's signed `Activation Z`) and explore the interactive zoomable/pannable **Mechanistic Network** cascade on the right.
- **Multi-Hop Causal Networks (CNA) Sub-Tab**:
  - Filter by root regulator, name, type, intermediate relay, or state (`All`, `Activated`, `Inhibited`); sort by `Master Root`, `Depth`, `State`, `Act Z`, `Overlap P`, or `Net-Bias P`; click **Export CSV** to export the CNA table.
  - Click any row to inspect shortest signed multi-hop causal paths (`K = 1..3`) from root regulators (even when transcriptionally unchanged, such as `SOST` in Cohort 3) to measured dataset DEGs.

### Pillar 4: Phenotype & Drug Match (`4`)
- **Multi-Cohort Comparison & Drug Reversers**: Compare Activation Z-scores across all loaded cohorts (`Upstream Regulators`, `Canonical Pathways`, `Downstream Effects`) alongside ranked **SigCom LINCS** perturbation connectivity scores (`Reverser <= -25` vs. `Mimicker >= +25`), with 1-click **Rescue** into Pathway Atlas.
- **Downstream Effects TreeMap**: Interactive hierarchical GO / Open Targets TreeMap sized by `-log10(P)` and colored by Activation Z, with **Export TreeMap SVG**.
- **Regulator Effects (3-Tier Flow)**: Zoomable causal flow explorer linking Upstream Regulators → shared consistent DEGs → Downstream Biological Processes / Diseases.

### Pillar 5: Knowledge Databases (`5`)
Monitor SQLite WAL file footprint, cached pathway/coordinate/edge counts, analysis cache hit/miss telemetry, live sync logs, and trigger full or per-stream ETag/SHA-256 updates (`WikiPathways GPML`, `Reactome GMT`, `OmniPath`, `CollecTRI`, `SIGNOR 3.0`, `Pathway Commons`, `STRING v12`, `GO GAF`, `Open Targets v4`, `SigCom LINCS`), local archive imports, `.ca-kb.json.gz` snapshots, cache clear, or `SQLite VACUUM`.

---

## 5. Exports (CSV, SVG, and JSON)

- **CSV Tables**: Available directly from the **Measured Gene Table**, **URA**, and **CNA** panels; exports obey your active filter and sort order.
- **Vector SVG (`Ctrl+E` or toolbar)**:
  - **Pathway MAP SVG**: Exports the subcellular compartment bands, `74x30` entity nodes, Bezier causal edges clipped cleanly outside `[cx±37, cy±15]` node boxes so `#arrow-act` / `#arrow-inh` arrowheads never collide with node borders, and a bottom **MAP Activity Colorbar & Edge Legend**.
  - **Downstream TreeMap SVG**: Exports the hierarchical GO / Open Targets TreeMap.
- **JSON Report (`Ctrl+S` or `analyze` CLI)**: Exports the complete `CoreAnalysisReport`.
