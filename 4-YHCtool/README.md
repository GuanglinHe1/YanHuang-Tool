# YHCtool

A configurable haplogroup population-analysis pipeline. All bundled sample identifiers, populations, regions, coordinates and STR values are synthetic and have no connection to real participants or locations.

## Run

Install R and the packages `yaml`, `readxl`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `ggrepel`, `FactoMineR`, `vegan`, `pheatmap`, `ggalluvial`, `scales`, `ape`, `rstatix`, `readr`, `ggbeeswarm`, `ggsignif`, `maps`, and `viridis`. The pipeline checks dependencies through normal R loading and never installs packages during analysis. `Rscript` must be on PATH for configuration loading; each step also declares its analysis executable under `tools.rscript`.

```bash
bash script/run_pipeline.sh
Rscript script/check_pipeline.R
```

The runner works from any working directory. To run one step, use `bash pipe/2-frequency.sh`, for example. Each step reads only its matching YAML in `conf/`. Run prerequisites before their consumers.

## Steps

| Step | Behavior |
| --- | --- |
| 1-preprocess | Validate sample identifiers and required columns, clean haplogroups, optionally equalize population sample sizes |
| 2-frequency | Counts, frequencies, diversity and effective diversity, prevalence, Fisher enrichment, haplogroup spread and bipartite edges |
| 3-distance | TVD, Jensen-Shannon distance, Nei D, DA, Bray-Curtis, Euclidean, Manhattan, optional Y-STR Rst, Mantel and PERMANOVA |
| 4-visualization | Distance heatmaps, hierarchical and NJ trees, MDS, PCA, CA, stacked bars, population-series and multilayer alluvial plots |
| 5-bootstrap | Resampled frequencies and diversity, pairwise tests, per-haplogroup plots and correlation-ordered facets |
| 6-geography | Cosine and JS similarity, consistency, geographic neighbours, heatmaps, metric and edge maps, all-pairs and target-population geographic Sankey plots |

## Inputs and configuration

`input/samples.tsv` contains `sample_id`, `haplogroup_full`, and `population`. Preprocessing also accepts CSV, TXT, XLSX and XLS, with the existing unambiguous column aliases `sample`/`id`, `haplogroup`/`hap`, and `pop`/`group`. Missing or ambiguous fields and duplicate sample IDs fail validation. Excel uses its first sheet. Additional columns are retained.

`input/populations.tsv` has one row per population, with `latitude`, `longitude`, and optional `region`. Coordinates must be complete and valid for geographic plotting. Optional `paths.population_order` in step 6 accepts a text file with one population per line. `input/ystr.tsv` contains `population` and numeric STR loci. Set `paths.metadata` or `paths.ystr` to an empty string to omit that optional input. A configured but missing file is an error. Geographic analyses require population coverage matching the frequency table.

Edit per-step `paths`, `tools`, `runtime`, and `analysis` values in `conf/`. Paths resolve relative to the project root. `analysis.levels` is a quoted comma-separated list of haplogroup prefix lengths. Unavailable lengths are skipped. `analysis.seed` fixes sampling and permutation seeds. `runtime.jobs` defaults to 1; current R stages execute serially. `runtime.overwrite` defaults to false and prevents replacing populated result directories. Set it to true explicitly for reruns. Synthetic-example bootstrap sizes and permutation counts are deliberately small; configure appropriate values for real analyses.

## Outputs

Machine-readable results are TSV files written without significant-digit rounding under `output/result/<step>/`. Plots are PDF only under `output/figure/<step>/`; each has a same-name TSV containing source and rendered layer data, or the matrix used by base graphics. Layer tables use `data_source` to distinguish source data from rendered layers. Newick trees are also exported. Arial AFM metrics are registered from `data/Arial.afm`; PDF text remains editable and fonts are not embedded. The viewer needs Arial or a compatible substitute. Font ownership notices in the metric file are retained.

Intermediates belong in `temp/<step>/`; logs belong in `log/`. Logs may contain paths supplied by the person running the pipeline and should be reviewed before sharing. Bootstrap output is separated by haplogroup level. The optional Mantel calculation preserves the existing Euclidean latitude/longitude convention; geographic neighbour calculations use great-circle distance. Bootstrap pairwise tests preserve the original descriptive resampling implementation and should not be interpreted as independent biological replication.

`script/check_pipeline.R` checks synthetic frequencies against independently tabulated counts, seven distance formulas, symmetry, bootstrap level separation, PDF sidecars and R syntax. Run it after the complete example pipeline. Input files are never rewritten.
