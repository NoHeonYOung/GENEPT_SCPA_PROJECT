# Phase 1B Vanilla SCPA reproduction summary

Gate status: `PASS`

## SCPA qval convention

SCPA qval is interpreted using the package convention: larger qval = stronger multivariate pathway difference. All rankings below are descending; rank 1 has the largest qval, and qval=0 belongs at the weakest end.

## Official workflow followed

The run uses SCPA::seurat_extract on the RNA assay's existing log1p-normalized data layer (pseudocount 0.001), then SCPA::compare_pathways with the official combined metabolic pathway collection.

## Compatibility adaptation

Serialized Seurat 3.1.5 was updated in memory to 5.4.0; the source RDS was not modified.

## Expression and pathways

- Assay/layer: `RNA/data`
- Matrix dimensions: 0h=17856x500, 12h=17856x500, 24h=17856x500, reference_resting_0h=17856x500, reference_activated_24h=17856x500
- Collection: SCPA combined metabolic pathways (Hallmark, KEGG, Reactome)
- Input/analyzed/excluded pathways: 243/124/119

## A. Hour-only analyses

These four retained analyses compare all naïve CD4 cells grouped only by observed Hour. They are valid time-point grouped analyses, not an exact paper Figure 4 or two-population tutorial reproduction.

### global_0_12_24 — top pathways by descending qval

| Rank | Pathway | qval |
| --- | --- | --- |
|  1 | REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES | 12.816029 |
|  2 | REACTOME_METABOLISM_OF_POLYAMINES | 10.246066 |
|  3 | REACTOME_RESPIRATORY_ELECTRON_TRANSPORT |  9.339898 |
|  4 | HALLMARK_HEME_METABOLISM |  9.013440 |
|  5 | HALLMARK_OXIDATIVE_PHOSPHORYLATION |  8.977158 |
|  6 | HALLMARK_GLYCOLYSIS |  8.940874 |
|  7 | KEGG_OXIDATIVE_PHOSPHORYLATION |  8.759424 |
|  8 | REACTOME_TP53_REGULATES_METABOLIC_GENES |  8.723128 |
|  9 | REACTOME_METABOLISM_OF_CARBOHYDRATES |  8.686830 |
| 10 | HALLMARK_FATTY_ACID_METABOLISM |  8.650530 |

### pairwise_0_vs_12 — top pathways by descending qval

| Rank | Pathway | qval | FC |
| --- | --- | --- | --- |
|  1 | REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES | 8.567192 | -20.956764 |
|  2 | REACTOME_METABOLISM_OF_POLYAMINES | 7.644537 | -20.517201 |
|  3 | KEGG_GLYCOLYSIS_GLUCONEOGENESIS | 6.803927 |  -9.846883 |
|  4 | REACTOME_METABOLISM_OF_CARBOHYDRATES | 6.719740 | -11.487476 |
|  5 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 6.551284 | -34.579799 |
|  6 | REACTOME_GLUCONEOGENESIS | 6.551284 |  -7.183667 |
|  7 | REACTOME_TRIGLYCERIDE_METABOLISM | 6.551284 |  -2.442427 |
|  8 | HALLMARK_GLYCOLYSIS | 6.467014 | -15.467057 |
|  9 | HALLMARK_FATTY_ACID_METABOLISM | 6.298379 |  -8.356905 |
| 10 | REACTOME_BIOLOGICAL_OXIDATIONS | 6.214013 |  -4.202281 |

### pairwise_12_vs_24 — top pathways by descending qval

| Rank | Pathway | qval | FC |
| --- | --- | --- | --- |
|  1 | REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES | 6.467014 |  -6.700752 |
|  2 | KEGG_OXIDATIVE_PHOSPHORYLATION | 4.858030 | -12.052417 |
|  3 | REACTOME_RESPIRATORY_ELECTRON_TRANSPORT | 4.516437 |  -9.304964 |
|  4 | REACTOME_TP53_REGULATES_METABOLIC_GENES | 4.430804 |  -4.670702 |
|  5 | REACTOME_METABOLISM_OF_POLYAMINES | 4.000897 |  -6.382151 |
|  6 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 3.914508 | -21.410441 |
|  7 | HALLMARK_HEME_METABOLISM | 3.827961 |  -1.570246 |
|  8 | HALLMARK_GLYCOLYSIS | 3.479927 |  -5.609338 |
|  9 | REACTOME_METABOLISM_OF_NUCLEOTIDES | 3.479927 |  -4.918973 |
| 10 | REACTOME_METABOLISM_OF_STEROIDS | 3.392380 |  -1.425836 |

### pairwise_0_vs_24 — top pathways by descending qval

| Rank | Pathway | qval | FC |
| --- | --- | --- | --- |
|  1 | REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES | 9.404697 | -27.657515 |
|  2 | REACTOME_METABOLISM_OF_CARBOHYDRATES | 8.231878 | -16.252576 |
|  3 | REACTOME_METABOLISM_OF_POLYAMINES | 8.231878 | -26.899352 |
|  4 | HALLMARK_FATTY_ACID_METABOLISM | 8.148017 | -16.059864 |
|  5 | REACTOME_GLUCOSE_METABOLISM | 8.148017 | -14.848377 |
|  6 | REACTOME_METABOLISM_OF_VITAMINS_AND_COFACTORS | 8.148017 |  -8.230564 |
|  7 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 7.980252 | -55.990240 |
|  8 | KEGG_GLYCOLYSIS_GLUCONEOGENESIS | 7.980252 | -13.643516 |
|  9 | KEGG_PURINE_METABOLISM | 7.980252 | -15.014238 |
| 10 | REACTOME_GLUCONEOGENESIS | 7.980252 | -10.102883 |

## B. Official two-population reference reproduction

- Population 1: `Cell_Type=Resting AND Hour=0`
- Population 2: `Cell_Type=Activated AND Hour=24`
- Full/sampled cells: 4341/500 vs 1697/500
- Seed/assay/layer: `20260810`, `RNA/data`
- Analyzed pathways and finite qval/FC: 124/124/124

### Reference top pathways by descending qval

| Rank | Pathway | qval | FC |
| --- | --- | --- | --- |
|  1 | REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES | 10.241332 | -33.07390 |
|  2 | HALLMARK_OXIDATIVE_PHOSPHORYLATION |  9.906766 | -69.42992 |
|  3 | REACTOME_METABOLISM_OF_CARBOHYDRATES |  9.906766 | -21.86758 |
|  4 | HALLMARK_FATTY_ACID_METABOLISM |  9.823108 | -18.65359 |
|  5 | HALLMARK_GLYCOLYSIS |  9.823108 | -28.33758 |
|  6 | REACTOME_GLUCONEOGENESIS |  9.823108 | -12.47004 |
|  7 | REACTOME_METABOLISM_OF_POLYAMINES |  9.823108 | -34.29539 |
|  8 | KEGG_GLYCOLYSIS_GLUCONEOGENESIS |  9.739441 | -17.35488 |
|  9 | REACTOME_GLYCOLYSIS |  9.739441 | -17.35399 |
| 10 | REACTOME_GLUCOSE_METABOLISM |  9.572086 | -19.72361 |

### Arachidonic-acid reference rows

| Rank | Pathway | qval | FC |
| --- | --- | --- | --- |
| 56 | REACTOME_ARACHIDONIC_ACID_METABOLISM | 6.635526 | -0.8938763 |
| 69 | KEGG_ARACHIDONIC_ACID_METABOLISM | 5.876184 | -0.2206279 |

- Qualitative agreement: `QUALITATIVELY_CONSISTENT`
- High-qval pathways with modest |FC| (<=5): 5; high qval does not require large FC: `TRUE`.
- The tutorial comparison is qualitative; numerical identity is not required. High qval with modest FC supports the SCPA multivariate interpretation.
- Any numerical difference may reflect the frozen package version, seed/downsampling, and the in-memory Seurat compatibility update; parameters were not retuned.

## Paper/tutorial comparison figures

- Visualization status: `PASS`
- Files: /home/node00/nhy_python/GenePT_SCPA/data/processed/genept_scpa/phase1/figures/01_global_qval_rank.png, /home/node00/nhy_python/GenePT_SCPA/data/processed/genept_scpa/phase1/figures/02_0_vs_24_enrichment_qval.png, /home/node00/nhy_python/GenePT_SCPA/data/processed/genept_scpa/phase1/figures/03_four_analysis_qval_heatmap.png, /home/node00/nhy_python/GenePT_SCPA/data/processed/genept_scpa/phase1/figures/phase1b_paper_comparison.png, /home/node00/nhy_python/GenePT_SCPA/data/processed/genept_scpa/phase1/figures/phase1b_paper_comparison.pdf, /home/node00/nhy_python/GenePT_SCPA/data/processed/genept_scpa/phase1/figures/phase1b_figure_notes.md, /home/node00/nhy_python/GenePT_SCPA/data/processed/genept_scpa/phase1/figures/05_reference_resting0_vs_activated24_qval_fc.png
- The 01–03/composite figures describe Hour-only analyses; the separate 05 figure describes the official Resting 0 h versus Activated 24 h reference. All comparisons are qualitative only.

## Hour-only qualitative comparison targets

The official tutorial highlights arachidonic-acid metabolism as a large multivariate change that need not have large mean enrichment. The multisample tutorial highlights glycolysis-related pathways. These are review targets, not tuning targets.

### Arachidonic-related rows in 0 vs 24

2 arachidonic-related row(s) were found; significance and FC direction require GPT review.

| Rank | Pathway | qval | FC |
| --- | --- | --- | --- |
| 42 | REACTOME_ARACHIDONIC_ACID_METABOLISM | 5.622378 | -1.0038580 |
| 51 | KEGG_ARACHIDONIC_ACID_METABOLISM | 4.943220 | -0.4915287 |

### Glycolysis-related rows in global 0/12/24

3 glycolysis-related row(s) were found; qval strength requires GPT review.

| Rank | Pathway | qval |
| --- | --- | --- |
|  6 | HALLMARK_GLYCOLYSIS | 8.940874 |
| 11 | REACTOME_GLYCOLYSIS | 8.214751 |
| 18 | KEGG_GLYCOLYSIS_GLUCONEOGENESIS | 7.269192 |

## Agreement, uncertainty, and review

- This run preserves official SCPA expression extraction, pathway collection, statistics, and thresholds.
- It is not a numerical replication of paper Figure 4: this protocol groups all cells by real Hour, whereas the paper/tutorial also use Cell_Type-specific or pseudotime-milestone populations.
- Global qval alone does not identify timing or direction; review it with all three pairwise outputs.
- Phase 1B PASS additionally requires the separate Resting 0 h versus Activated 24 h reference run and its qualitative review.
- No parameter was tuned after viewing results.
- Final gate basis: `successful pipeline execution plus qualitative reference reproduction, not numerical identity`.
- Runtime versions: SCPA 1.6.2; Seurat 5.5.1; SeuratObject 5.4.0.
