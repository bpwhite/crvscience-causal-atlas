# Headless Validation (`crvscience-causal-atlas` v0.2.0)

Deterministic, dry-lab validation suite for [`ca-core`](../crates/ca-core) executed via compiled Rust integration tests without GUI or network side effects.

## Command

```bash
cargo test --workspace
# Or run only the headless validation integration suite:
cargo test -p ca-core --test headless_validation -- --nocapture
```

## Test Coverage

1. **T1 (`t1_deseq2_wald_and_bh_fdr_calibration_on_huvec_tnf_cohort`)**: Parses `examples/tnf_huvec_deg.tsv` (`parse_expression_table`) and evaluates the built-in `HUVEC + TNF Inflammatory Cascade` cohort (`28` measured genes), verifying Benjamini-Hochberg step-up FDR monotonicity (`benjamini_hochberg`), Wald Z-to-two-sided `p`-value roundtrip precision (`z_to_two_sided_p`), and exclusion of non-significant housekeeping/unrelated controls (`GAPDH`, `ACTB`, `TFF1`, `GREB1`).
2. **T2 (`t2_kramer_causal_z_and_fisher_upstream_regulators_ranking`)**: Runs Krämer et al. (2014) signed Upstream Regulator Analysis (`run_upstream_regulators`) + Fisher's exact enrichment across the seed causal knowledge graph, confirming significant activation of `TNF`, `NFKB1`, and `RELA`, specificity against the `MCF-7 + Beta-Estradiol Signaling` contrast cohort, and writing all scored upstream regulators to [`upstream_regulators_ranked.csv`](upstream_regulators_ranked.csv).
3. **T3 (`t3_signed_mechanistic_paths_and_cmap_perturbation_reversal`)**: Constructs the TNF transitive-triangle Mechanistic Network (`build_mechanistic_network`), executes multi-hop Causal Network Analysis (`run_causal_networks` with 1,000 degree-preserving Monte Carlo permutations, including recovering transcriptionally silent root `SOST` in Cohort 3), and scores CMap/L1000 therapeutic reversers vs. mimickers (`run_perturbation_match`).

## Generated Validation Artifacts

- [`causal_atlas_validation_summary.json`](causal_atlas_validation_summary.json) — Structured JSON telemetry for T1–T3 (mirrored to `app-preprints/releases/causal-atlas-pathway-perturbation-validation-telemetry.json`).
- [`upstream_regulators_ranked.csv`](upstream_regulators_ranked.csv) — Ranked Upstream Regulator Analysis table (`rank`, `activation_z`, `bias_corrected_z`, Fisher `p`, BH `fdr`, mechanistic triangle summaries).
