# Phase 1B figure comparison notes

These figures are designed for qualitative comparison with SCPA paper Figure 4 and the official tutorials.

- SCPA qval uses the package convention: larger qval means a stronger multivariate pathway difference; qval=0 is at the weakest end.
- Panel A follows the official global qval/rank presentation and highlights glycolysis-related pathways.
- Panel B follows the official `-FC` versus `qval` presentation and highlights arachidonic-related pathways.
- Panel C is an added diagnostic heatmap comparing qval across the requested global and pairwise analyses.
- Pairwise FC is population 1 minus population 2; therefore Panel B uses `-FC`, so positive x indicates enrichment toward 24 h.
- This is not a numerical reproduction of paper Figure 4. The present protocol uses all cells grouped by Hour; the paper/tutorial also use Cell_Type-specific or pseudotime-milestone populations.
- Compare pathway rank and qualitative signal, not exact coordinates or pathway ordering.
