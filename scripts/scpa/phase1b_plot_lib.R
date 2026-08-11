phase1b_display_name <- function(pathway) {
  text <- gsub("_", " ", pathway, fixed = TRUE)
  vapply(text, function(value) paste(strwrap(value, width = 48L), collapse = "\n"), character(1))
}

phase1b_required_results <- function() {
  c("global_0_12_24", "pairwise_0_vs_12", "pairwise_12_vs_24", "pairwise_0_vs_24")
}

build_phase1b_comparison_plots <- function(results, top_n = 30L) {
  required <- phase1b_required_results()
  missing <- setdiff(required, names(results))
  if (length(missing) > 0L) {
    stop("Missing results required for Phase 1B figures: ", paste(missing, collapse = ", "))
  }

  global <- results$global_0_12_24
  global <- global[order(global$qval, decreasing = TRUE), , drop = FALSE]
  global$Rank <- seq_len(nrow(global))
  global$Target <- grepl("GLYCOLYSIS", global$Pathway, ignore.case = TRUE)
  glycolysis_label <- utils::head(global[global$Target, , drop = FALSE], 1L)

  p_global <- ggplot2::ggplot(global, ggplot2::aes(x = Rank, y = qval)) +
    ggplot2::geom_point(shape = 21, size = 2.1, stroke = 0.25, fill = "grey70") +
    ggplot2::geom_point(
      data = global[global$Target, , drop = FALSE],
      shape = 21, size = 2.7, stroke = 0.35, fill = "#D7301F"
    ) +
    ggplot2::labs(
      title = "A. Global 0/12/24 pathway rank",
      subtitle = "Higher qval = stronger; glycolysis highlighted",
      x = "Pathway rank (highest qval first)", y = "SCPA qval"
    ) +
    ggplot2::theme_classic(base_size = 11)
  if (nrow(glycolysis_label) > 0L) {
    p_global <- p_global + ggplot2::geom_text(
      data = glycolysis_label,
      ggplot2::aes(label = phase1b_display_name(Pathway)),
      hjust = 0, nudge_x = 3, size = 3, check_overlap = TRUE
    )
  }

  pair_024 <- results$pairwise_0_vs_24
  pair_024$plot_x <- -pair_024$FC
  pair_024$Target <- grepl("ARACHI", pair_024$Pathway, ignore.case = TRUE)
  pair_024$Signal <- ifelse(
    pair_024$adjPval < 0.01 & abs(pair_024$FC) > 5,
    "adjP < 0.01, |FC| > 5",
    ifelse(pair_024$adjPval < 0.01, "adjP < 0.01, |FC| <= 5", "Other")
  )
  pair_024$Signal <- factor(
    pair_024$Signal,
    levels = c("Other", "adjP < 0.01, |FC| <= 5", "adjP < 0.01, |FC| > 5")
  )
  arachidonic_label <- pair_024[pair_024$Target, , drop = FALSE]
  arachidonic_label <- utils::head(
    arachidonic_label[order(arachidonic_label$qval, decreasing = TRUE), , drop = FALSE],
    1L
  )

  p_scatter <- ggplot2::ggplot(pair_024, ggplot2::aes(x = plot_x, y = qval)) +
    ggplot2::geom_vline(xintercept = c(-5, 5), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_point(
      ggplot2::aes(fill = Signal), shape = 21, size = 2.4, stroke = 0.3
    ) +
    ggplot2::geom_point(
      data = pair_024[pair_024$Target, , drop = FALSE],
      shape = 21, size = 3, stroke = 0.4, fill = "#FF4500", color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Other" = "black",
        "adjP < 0.01, |FC| <= 5" = "#84B0F0",
        "adjP < 0.01, |FC| > 5" = "#6DBF88"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "B. Pairwise 0 h vs 24 h",
      subtitle = "Higher qval = stronger; arachidonic highlighted",
      x = "Enrichment toward 24 h (-FC)", y = "SCPA qval", fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  if (nrow(arachidonic_label) > 0L) {
    p_scatter <- p_scatter + ggplot2::geom_text(
      data = arachidonic_label,
      ggplot2::aes(label = phase1b_display_name(Pathway)),
      hjust = 0, nudge_x = 1, nudge_y = 0.15, size = 3, check_overlap = TRUE
    )
  }

  analysis_labels <- c(
    global_0_12_24 = "Global\n0/12/24",
    pairwise_0_vs_12 = "0 vs 12",
    pairwise_12_vs_24 = "12 vs 24",
    pairwise_0_vs_24 = "0 vs 24"
  )
  long <- do.call(rbind, lapply(required, function(analysis_name) {
    data.frame(
      Pathway = results[[analysis_name]]$Pathway,
      Analysis = analysis_name,
      qval = results[[analysis_name]]$qval,
      stringsAsFactors = FALSE
    )
  }))
  maximum_qval <- tapply(long$qval, long$Pathway, max, na.rm = TRUE)
  ordered_pathways <- names(sort(maximum_qval, decreasing = TRUE))
  selected_pathways <- utils::head(ordered_pathways, min(as.integer(top_n), length(ordered_pathways)))
  target_pathways <- unique(long$Pathway[grepl("GLYCOLYSIS|ARACHI", long$Pathway, ignore.case = TRUE)])
  selected_pathways <- unique(c(selected_pathways, target_pathways))
  heat <- long[long$Pathway %in% selected_pathways, , drop = FALSE]
  heat$Analysis <- factor(heat$Analysis, levels = required, labels = unname(analysis_labels[required]))
  heat$Pathway <- factor(heat$Pathway, levels = rev(selected_pathways))

  p_heatmap <- ggplot2::ggplot(heat, ggplot2::aes(x = Analysis, y = Pathway, fill = qval)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B") +
    ggplot2::scale_y_discrete(labels = phase1b_display_name) +
    ggplot2::labs(
      title = "C. qval comparison across all four analyses",
      subtitle = paste0("Top ", top_n, " by maximum qval, plus glycolysis/arachidonic targets"),
      x = NULL, y = NULL, fill = "qval"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7.5)
    )

  combined <- patchwork::wrap_plots(
    p_global, p_scatter, p_heatmap,
    design = "AB\nCC",
    heights = c(1, 1.8)
  ) + patchwork::plot_annotation(
    title = "Phase 1B Vanilla SCPA - comparison with the paper/tutorial",
    subtitle = paste(
      "Qualitative comparison only: this run uses all cells grouped by real Hour;",
      "the paper uses Cell_Type-specific and pseudotime-milestone populations."
    )
  )

  list(global_rank = p_global, pairwise_0_vs_24 = p_scatter, four_analysis_heatmap = p_heatmap, combined = combined)
}

save_phase1b_plot_atomic <- function(plot, path, width, height, dpi = 300L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  extension <- tools::file_ext(path)
  temporary <- tempfile(pattern = ".phase1b_plot_", tmpdir = dirname(path), fileext = paste0(".", extension))
  on.exit(unlink(temporary), add = TRUE)
  ggplot2::ggsave(temporary, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  if (!file.rename(temporary, path)) stop("Could not atomically write figure: ", path)
  path
}

render_phase1b_figures <- function(results, output_dir, top_n = 30L) {
  plots <- build_phase1b_comparison_plots(results, top_n = top_n)
  files <- c(
    global_rank_png = save_phase1b_plot_atomic(
      plots$global_rank, file.path(output_dir, "01_global_qval_rank.png"), 7, 5
    ),
    pairwise_0_vs_24_png = save_phase1b_plot_atomic(
      plots$pairwise_0_vs_24, file.path(output_dir, "02_0_vs_24_enrichment_qval.png"), 7, 5
    ),
    four_analysis_heatmap_png = save_phase1b_plot_atomic(
      plots$four_analysis_heatmap, file.path(output_dir, "03_four_analysis_qval_heatmap.png"), 8, 10
    ),
    combined_png = save_phase1b_plot_atomic(
      plots$combined, file.path(output_dir, "phase1b_paper_comparison.png"), 14, 13
    ),
    combined_pdf = save_phase1b_plot_atomic(
      plots$combined, file.path(output_dir, "phase1b_paper_comparison.pdf"), 14, 13
    )
  )

  notes_path <- file.path(output_dir, "phase1b_figure_notes.md")
  notes <- c(
    "# Phase 1B figure comparison notes",
    "",
    "These figures are designed for qualitative comparison with SCPA paper Figure 4 and the official tutorials.",
    "",
    "- SCPA qval uses the package convention: larger qval means a stronger multivariate pathway difference; qval=0 is at the weakest end.",
    "- Panel A follows the official global qval/rank presentation and highlights glycolysis-related pathways.",
    "- Panel B follows the official `-FC` versus `qval` presentation and highlights arachidonic-related pathways.",
    "- Panel C is an added diagnostic heatmap comparing qval across the requested global and pairwise analyses.",
    "- Pairwise FC is population 1 minus population 2; therefore Panel B uses `-FC`, so positive x indicates enrichment toward 24 h.",
    "- This is not a numerical reproduction of paper Figure 4. The present protocol uses all cells grouped by Hour; the paper/tutorial also use Cell_Type-specific or pseudotime-milestone populations.",
    "- Compare pathway rank and qualitative signal, not exact coordinates or pathway ordering."
  )
  temporary <- tempfile(pattern = ".phase1b_figure_notes_", tmpdir = output_dir)
  on.exit(unlink(temporary), add = TRUE)
  writeLines(notes, temporary)
  if (!file.rename(temporary, notes_path)) stop("Could not atomically write figure notes")
  c(files, figure_notes = notes_path)
}

build_phase1b_reference_plot <- function(result) {
  failures <- validate_scpa_result(result, pairwise = TRUE)
  if (length(failures) > 0L) {
    stop("Invalid reference result for plotting: ", paste(failures, collapse = ", "))
  }
  targets <- c(
    "REACTOME_ARACHIDONIC_ACID_METABOLISM",
    "KEGG_ARACHIDONIC_ACID_METABOLISM"
  )
  plotted <- result
  plotted$plot_x <- -plotted$FC
  plotted$Target <- plotted$Pathway %in% targets
  plotted$TargetName <- ifelse(
    plotted$Pathway == targets[[1]], "Reactome arachidonic acid metabolism",
    ifelse(plotted$Pathway == targets[[2]], "KEGG arachidonic acid metabolism", NA_character_)
  )
  ggplot2::ggplot(plotted, ggplot2::aes(x = plot_x, y = qval)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey65", linewidth = 0.3) +
    ggplot2::geom_point(shape = 21, size = 2.5, stroke = 0.3, fill = "grey65") +
    ggplot2::geom_point(
      data = plotted[plotted$Target, , drop = FALSE],
      ggplot2::aes(fill = TargetName), shape = 21, size = 3.6, stroke = 0.5, color = "black"
    ) +
    ggplot2::geom_text(
      data = plotted[plotted$Target, , drop = FALSE],
      ggplot2::aes(label = phase1b_display_name(Pathway)),
      hjust = 0, nudge_x = 0.3, nudge_y = 0.1, size = 3, check_overlap = TRUE
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Reactome arachidonic acid metabolism" = "#D7301F",
      "KEGG arachidonic acid metabolism" = "#4575B4"
    )) +
    ggplot2::labs(
      title = "Official two-population reference: Resting 0 h vs Activated 24 h",
      subtitle = "Larger qval = stronger SCPA difference; qualitative reproduction only",
      x = "Enrichment toward Activated 24 h (-FC)",
      y = "SCPA qval",
      fill = NULL,
      caption = paste(
        "Population 1: Resting 0 h; population 2: Activated 24 h.",
        "Positive x indicates enrichment toward Activated 24 h."
      )
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

render_phase1b_reference_figure <- function(result, output_dir) {
  plot <- build_phase1b_reference_plot(result)
  save_phase1b_plot_atomic(
    plot,
    file.path(output_dir, "05_reference_resting0_vs_activated24_qval_fc.png"),
    width = 8,
    height = 6
  )
}
