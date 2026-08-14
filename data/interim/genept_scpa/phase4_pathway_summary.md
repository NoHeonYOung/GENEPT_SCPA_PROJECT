# Phase 4 pathway comparison summary

Gate status: `READY_FOR_GPT_REVIEW`

## 1. Research question

How does a pathway-specific GenePT semantic projection change pathway rankings when Vanilla and GenePT-informed branches use the same cells, pathways, and genes?

## 2. Why Phase 4 follows Phase 3

Phase 3 established that whole-cell GenePT-w preserves a detectable CD4/CD8 multivariate difference. Phase 4 localizes the comparison to curated pathways without treating the 1,536 semantic dimensions as genes.

## 3. Vanilla pathway method

For each pathway, the 500 x p log-expression matrices for CD4 and CD8 are compared with `multicross::mcm()`.

## 4. GenePT-informed pathway method

The same X_P matrices are projected as Z_P = X_P x E_P, yielding 500 x 1,536 matrices, and compared with the same SCPA-core MCM function.

## 5. Paired gene-set policy

Both branches use pathway genes present in CD4, present in CD8, and exactly mappable to the official GenePT artifact. Expression is normalized over the full transcriptome before pathway subsetting; pathway-local renormalization is prohibited.

## 6. Pathways analyzed

Input/eligible/excluded: 243/123/120. Frozen min/max paired genes: 15/500.

## 7. Rank agreement

Spearman=0.682166; Kendall=0.672931; Top-10 overlap=6; Top-20 overlap=12. These are agreement metrics, not accuracy metrics.

## 8. Largest rank shifts

`rank_delta = genept_rank - vanilla_rank`; negative values move upward after GenePT-informed projection.

### Upward shifts

| Pathway | Vanilla rank | GenePT rank | Rank delta |
| --- | ---: | ---: | ---: |
| REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA | 109 | 12 | -97 |
| REACTOME_METABOLISM_OF_POLYAMINES | 95 | 10 | -85 |
| KEGG_OXIDATIVE_PHOSPHORYLATION | 59 | 13 | -46 |
| HALLMARK_GLYCOLYSIS | 32 | 4 | -28 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 25 | 2 | -23 |
| KEGG_ARGININE_AND_PROLINE_METABOLISM | 36 | 16 | -20 |
| REACTOME_NUCLEOBASE_CATABOLISM | 29 | 11 | -18 |
| HALLMARK_HEME_METABOLISM | 33 | 22 | -11 |
| KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM | 34 | 23 | -11 |
| KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM | 35 | 24 | -11 |

### Downward shifts

| Pathway | Vanilla rank | GenePT rank | Rank delta |
| --- | ---: | ---: | ---: |
| REACTOME_TP53_REGULATES_METABOLIC_GENES | 11 | 120 | 109 |
| REACTOME_TRIGLYCERIDE_METABOLISM | 17 | 122 | 105 |
| REACTOME_SYNTHESIS_OF_IP3_AND_IP4_IN_THE_CYTOSOL | 22 | 113 | 91 |
| REACTOME_SYNTHESIS_OF_ACTIVE_UBIQUITIN_ROLES_OF_E1_AND_E2_ENZYMES | 21 | 110 | 89 |
| REACTOME_PEPTIDE_HORMONE_METABOLISM | 16 | 98 | 82 |
| REACTOME_RESPIRATORY_ELECTRON_TRANSPORT | 30 | 105 | 75 |
| REACTOME_GLUCONEOGENESIS | 3 | 75 | 72 |
| REACTOME_GLYCOGEN_METABOLISM | 9 | 78 | 69 |
| REACTOME_DISEASES_OF_CARBOHYDRATE_METABOLISM | 12 | 70 | 58 |
| REACTOME_GLYCOLYSIS | 28 | 79 | 51 |

## 9. Example pathways

| Pathway | Vanilla rank | GenePT rank | Rank delta |
| --- | ---: | ---: | ---: |
| REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES | 1 | 1 | 0 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 25 | 2 | -23 |
| KEGG_PURINE_METABOLISM | 2 | 3 | 1 |
| REACTOME_GLUCONEOGENESIS | 3 | 75 | 72 |
| HALLMARK_GLYCOLYSIS | 32 | 4 | -28 |
| REACTOME_METABOLISM_OF_CARBOHYDRATES | 4 | 20 | 16 |
| HALLMARK_FATTY_ACID_METABOLISM | 5 | 8 | 3 |
| REACTOME_SYNTHESIS_OF_PC | 7 | 5 | -2 |
| KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450 | 6 | 6 | 0 |
| REACTOME_GLYCOSAMINOGLYCAN_METABOLISM | 10 | 7 | -3 |

## 10. What the result means

The result quantifies agreement and relative pathway reordering after the semantic projection.

## 11. What the result DOES NOT mean

Smaller p-values, larger qval values, or upward ranks do not establish that GenePT is better or more accurate. The representation geometries differ.

## 12. Next gene-level analysis

The saved manifest preserves pathway genes, paired genes, feature order, embedding keys, match types, canonical cells, and preprocessing so Phase 5 can regenerate inputs. No gene masking or leave-one-gene-out analysis was run.

## 13. Semantic-control plan

Phase 6 will compare True, gene-to-embedding Permuted, and dimension-matched Random embeddings with repeated-sampling robustness. These controls were not run in Phase 4.
