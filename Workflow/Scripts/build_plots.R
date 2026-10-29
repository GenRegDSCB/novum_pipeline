suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(scales)
})


#----- Parameters from snakemake@params -----#
plot_type      <- snakemake@params$type
basesize       <- as.integer(snakemake@params$basesize %||% 12)
pvalue_raw     <- snakemake@params$pvalue_threshold
pvalue_cutoff  <- if (length(pvalue_raw) == 0) NA_real_ else as.numeric(pvalue_raw)
has_pvalue_cut <- !is.na(pvalue_cutoff)
locus_tags     <- snakemake@params$locus   %||% character(0)
gene_names     <- snakemake@params$gene    %||% character(0)
protein_names  <- snakemake@params$protein %||% character(0)
variant_labels <- snakemake@params$variant_labels
per_mirna_top_n <- as.integer(snakemake@params$per_mirna_top_n %||% 12L)


#----- Read each variant CSV, tagging it with the label provided by Snakemake. trim_ws off:
#      the four duplex strings are space-padded into column alignment -----#
read_one <- function(path, label) {
  read_csv(path, show_col_types = FALSE, trim_ws = FALSE) |> mutate(Calibration = label)
}

merged <- map2(snakemake@input$annotated, variant_labels, read_one) |>
  bind_rows() |>
  mutate(Calibration = factor(Calibration, levels = c("uncalibrated", "calibrated")))


#----- Validated-hits filter: each entry is "<target>" (any miRNA) or "<target>,<miRNA>" (only that miRNA) -----#
parse_pairs <- function(entries, target_col) {
  if (length(entries) == 0) {
    return(tibble(!!target_col := character(0), miRNA_filter = character(0)))
  }
  parts <- strsplit(as.character(entries), ",", fixed = TRUE)
  if (any(lengths(parts) > 2L)) {
    bad <- entries[lengths(parts) > 2L][1L]
    stop(sprintf("build_plots %s: expected '<target>' or '<target>,<miRNA>', got: %s", target_col, bad))
  }
  target <- trimws(vapply(parts, `[[`, character(1), 1L))
  mirna  <- trimws(vapply(parts, function(p) if (length(p) > 1L) p[[2L]] else NA_character_, character(1)))
  if (any(target == ""))                stop(sprintf("build_plots %s: empty target before comma in entry: %s", target_col, entries[target == ""][1L]))
  if (any(!is.na(mirna) & mirna == "")) stop(sprintf("build_plots %s: empty miRNA after comma in entry: %s",  target_col, entries[!is.na(mirna) & mirna == ""][1L]))
  tibble(!!target_col := target, miRNA_filter = mirna)
}

locus_pairs   <- parse_pairs(locus_tags,    "locus_tag")
gene_pairs    <- parse_pairs(gene_names,    "gene_name")
protein_pairs <- parse_pairs(protein_names, "protein_name")

if (nrow(protein_pairs) > 0L && !"protein_name" %in% names(merged)) {
  stop("build_plots protein: input CSV has no `protein_name` column.")
}

#----- Keeps the rows matching a (target, miRNA) pair list; an empty list selects nothing -----#
join_filter <- function(df, pairs, by_col) {
  if (nrow(pairs) == 0L) return(df[0, , drop = FALSE])
  df |>
    inner_join(pairs, by = by_col, relationship = "many-to-many") |>
    filter(is.na(miRNA_filter) | miRNA == miRNA_filter) |>
    select(-miRNA_filter)
}

validated <- bind_rows(
  join_filter(merged, locus_pairs,   "locus_tag"),
  join_filter(merged, gene_pairs,    "gene_name"),
  join_filter(merged, protein_pairs, "protein_name")
)

# One row per (miRNA, locus_tag, Calibration) for the p-value distribution vlines
# and for the legend/scale keys. Position-based plots use the full `validated` set.
validated_best <- validated |>
  group_by(Calibration) |>
  arrange(P_value) |>
  distinct(miRNA, locus_tag, .keep_all = TRUE) |>
  ungroup()

# Empty validated set → drop the per-hit highlight layer (red points / dashed
# vlines) but still render everything else. Guards `scales::hue_pal()(0)` below,
# which would otherwise abort the whole job. Warn only when a hit was actually
# requested — an empty locus/gene/protein config means none was wanted.
has_validated <- nrow(validated_best) > 0L
requested <- c(
  if (length(locus_tags))    paste0("locus=",   toString(locus_tags)),
  if (length(gene_names))    paste0("gene=",    toString(gene_names)),
  if (length(protein_names)) paste0("protein=", toString(protein_names))
)
if (!has_validated && length(requested) > 0L) {
  warning(sprintf(
    "build_plots: validated-hits filter matched 0 rows (%s) — rendering plots without the highlight layer.",
    paste(requested, collapse = "; ")
  ), call. = FALSE)
}

#----- Auto color + label maps for the validated hits, keyed by locus_tag -----#
# Label cascade for readability: gene_name > protein_name > locus_tag.
if (!"protein_name" %in% names(validated_best)) {
  validated_best$protein_name <- NA_character_
}
hit_levels <- unique(validated_best$locus_tag)
hit_colors <- if (has_validated) setNames(scales::hue_pal()(length(hit_levels)), hit_levels) else character(0)
hit_label_table <- validated_best |>
  distinct(locus_tag, miRNA, gene_name, protein_name) |>
  mutate(
    gene_name    = na_if(as.character(gene_name),    ""),
    protein_name = na_if(as.character(protein_name), ""),
    locus_tag    = as.character(locus_tag),
    label        = sprintf("%s | %s", miRNA, coalesce(gene_name, protein_name, locus_tag))
  )
hit_labels <- setNames(hit_label_table$label, hit_label_table$locus_tag)
hit_shapes <- setNames(seq_along(hit_levels) + 14L, hit_levels)

#----- Per-hit highlight layers: the geom + its manual scale, or NULL when there
# are no validated hits (so `plot + highlight_*()` is a no-op). -----#
highlight_vlines <- function() {
  if (!has_validated) return(NULL)
  list(
    geom_vline(
      data = validated_best,
      aes(xintercept = P_value, color = locus_tag),
      linetype = "dashed"
    ),
    scale_color_manual(name = "Validated hits", values = hit_colors, labels = hit_labels)
  )
}
#----- Red points marking the validated hits; NULL when none were requested -----#
highlight_points <- function(data, ...) {
  if (!has_validated) return(NULL)
  list(
    geom_point(data = data, mapping = aes(shape = locus_tag), color = "red", ...),
    scale_shape_manual(name = "Validated hits", values = hit_shapes, labels = hit_labels)
  )
}


#----- Shared Calibration aesthetics -----#
calibration_colors <- c(uncalibrated = "#e41a1c", calibrated = "#377eb8")
calibration_labels <- c(uncalibrated = "Without Calibration", calibrated = "With Calibration")
facet_labels       <- as_labeller(calibration_labels)

#----- Applies the configured p-value cutoff, or passes the frame through when none is set -----#
apply_pvalue_cut <- function(df) {
  if (has_pvalue_cut) df |> filter(P_value <= pvalue_cutoff) else df
}

#----- Title helper: use the p-value-cut variant of a plot title when a cutoff is active -----#
cut_title <- function(plain, cut_fmt) {
  if (has_pvalue_cut) sprintf(cut_fmt, pvalue_cutoff) else plain
}


#----- Plot builders -----#
#----- P-value density per miRNA, one facet per Calibration variant -----#
plot_pvalue_distribution <- function() {
  dens <- merged |>
    group_by(Calibration) |>
    summarise(ymax = max(density(P_value, na.rm = TRUE)$y), .groups = "drop")
  ycut     <- ceiling(min(dens$ymax) * 3.5)
  ynatural <- ceiling(max(dens$ymax) * 1.05)

  ggplot(merged, aes(x = P_value)) +
    geom_density(color = "black", fill = "#d1d1d1") +
    highlight_vlines() +
    coord_cartesian(ylim = c(0, ycut)) +
    facet_wrap(~Calibration, labeller = facet_labels) +
    labs(
      x = "p-value", y = "Density",
      title = sprintf("Distribution of p-values (uncut Y-axis: %d)", ynatural)
    ) +
    theme_bw(base_size = basesize)
}

#----- P-value against relative position along the gene -----#
plot_position_pvalue <- function() {
  ggplot(merged, aes(x = Position, y = -log10(P_value))) +
    geom_point(
      data = apply_pvalue_cut(merged),
      color = "black", size = 0.7
    ) +
    highlight_points(validated) +
    geom_smooth(method = "lm", linewidth = 0.5) +
    stat_cor(method = "spearman") +
    labs(
      x = "Relative Position", y = "-log10(p-value)",
      title = cut_title("Position vs p-values", "Position vs p-values (≤ %g)")
    ) +
    theme_bw(base_size = basesize) +
    facet_wrap(~Calibration, labeller = facet_labels)
}

#----- Hybridization energy against relative position, with a trend fit -----#
plot_position_energy <- function() {
  pick <- if (any(merged$Calibration == "calibrated")) "calibrated" else "uncalibrated"
  base_data <- merged |>
    filter(Calibration == pick) |>
    mutate(facet_label = calibration_labels[[pick]])

  ggplot(base_data, aes(x = Position, y = Energy)) +
    geom_point(
      data = apply_pvalue_cut(base_data),
      color = "black", size = 0.7
    ) +
    highlight_points(validated |> filter(Calibration == pick)) +
    geom_smooth(method = "lm", linewidth = 0.5) +
    stat_cor(method = "spearman") +
    labs(
      x = "Relative Position", y = "Energy (KCal/Mol)",
      title = cut_title("Position vs Energy", "Position vs Energy (p-value ≤ %g)")
    ) +
    facet_wrap(~facet_label) +
    theme_bw(base_size = basesize)
}

#----- Volcano: hybridization energy against p-value -----#
plot_volcano <- function() {
  ggplot(merged, aes(x = Energy, y = -log10(P_value))) +
    geom_point(
      data = apply_pvalue_cut(merged),
      color = "black", size = 0.6, alpha = 0.4
    ) +
    highlight_points(validated, size = 2) +
    facet_wrap(~Calibration, labeller = facet_labels) +
    labs(
      x = "Energy (KCal/Mol)", y = "-log10(p-value)",
      title = cut_title("Volcano: Energy vs p-value", "Volcano (p ≤ %g)")
    ) +
    theme_bw(base_size = basesize)
}

#----- Empirical CDF of p-values, comparing the Calibration variants -----#
plot_pvalue_ecdf <- function() {
  data <- merged |> mutate(facet_label = "p-value ECDF vs Uniform")
  ggplot(data, aes(x = P_value, color = Calibration)) +
    stat_ecdf(geom = "step", linewidth = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(name = NULL, values = calibration_colors, labels = calibration_labels) +
    labs(
      x = "p-value", y = "Empirical CDF",
      title = "Calibration QC"
    ) +
    facet_wrap(~facet_label) +
    theme_bw(base_size = basesize)
}

#----- Per-hit calibrated-vs-uncalibrated p-value shift; NULL unless both variants ran -----#
plot_calibration_delta <- function() {
  variants_present <- as.character(unique(merged$Calibration))
  if (!all(c("uncalibrated", "calibrated") %in% variants_present)) {
    return(NULL)
  }

  paired <- merged |>
    select(miRNA, Gene, locus_tag, gene_name, protein_name, Position, Calibration, P_value) |>
    group_by(miRNA, Gene, locus_tag, gene_name, protein_name, Position, Calibration) |>
    summarise(P_value = min(P_value), .groups = "drop") |>
    pivot_wider(
      names_from  = Calibration,
      values_from = P_value
    ) |>
    filter(!is.na(uncalibrated) & !is.na(calibrated)) |>
    mutate(uncal_log = -log10(uncalibrated), cal_log = -log10(calibrated))

  validated_paired <- paired |> semi_join(validated_best, by = c("miRNA", "locus_tag"))

  paired           <- paired           |> mutate(facet_label = "Calibrated vs Uncalibrated")
  validated_paired <- validated_paired |> mutate(facet_label = "Calibrated vs Uncalibrated")

  ggplot(paired, aes(x = uncal_log, y = cal_log)) +
    geom_point(color = "black", size = 0.4, alpha = 0.25) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
    highlight_points(validated_paired, size = 2.5) +
    labs(
      x = "-log10(p-value), uncalibrated",
      y = "-log10(p-value), calibrated",
      title = "Per-hit calibration impact"
    ) +
    facet_wrap(~facet_label) +
    theme_bw(base_size = basesize)
}

#----- Top-N hits per miRNA, ranked within each Calibration variant -----#
plot_per_mirna <- function() {
  top_mirnas <- merged |>
    count(miRNA, sort = TRUE) |>
    slice_head(n = per_mirna_top_n) |>
    pull(miRNA)
  data <- merged |>
    filter(miRNA %in% top_mirnas) |>
    mutate(miRNA = factor(miRNA, levels = top_mirnas))

  ggplot(data, aes(x = Calibration, y = -log10(P_value), fill = Calibration)) +
    geom_violin(alpha = 0.55, color = NA) +
    geom_boxplot(width = 0.18, alpha = 0.85, outlier.size = 0.3) +
    facet_wrap(~miRNA, scales = "free_y") +
    scale_fill_manual(values = calibration_colors, labels = calibration_labels) +
    labs(
      x = NULL, y = "-log10(p-value)",
      title = sprintf("Per-miRNA p-value distribution (top %d by hit count)", length(top_mirnas))
    ) +
    theme_bw(base_size = basesize) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
}

#----- Density of hit positions along the gene -----#
plot_position_density <- function() {
  if (has_pvalue_cut) {
    data <- merged |>
      mutate(p_tier = factor(
        if_else(P_value <= pvalue_cutoff, "significant", "non-significant"),
        levels = c("non-significant", "significant")
      ))
    aes_used  <- aes(x = Position, fill = p_tier)
    fill_name <- sprintf("p ≤ %g", pvalue_cutoff)
    fill_vals <- c("non-significant" = "#bdbdbd", "significant" = "#377eb8")
  } else {
    data <- merged
    aes_used  <- aes(x = Position, fill = Calibration)
    fill_name <- NULL
    fill_vals <- calibration_colors
  }
  ggplot(data, aes_used) +
    geom_histogram(binwidth = 0.02, alpha = 0.8, position = "stack") +
    facet_wrap(~Calibration, labeller = facet_labels) +
    scale_fill_manual(name = fill_name, values = fill_vals) +
    labs(
      x = "Relative position along gene (0 = 5', 1 = 3')",
      y = "Hit count",
      title = "Distribution of hit positions"
    ) +
    theme_bw(base_size = basesize)
}

#----- Canonical-seed class (6mer / 7mer-A1 / 7mer-m8 / 8mer / none). The miRNA prints 3'->5',
#      so nt 1 is the rightmost miRNA-bearing column; walking left over those columns yields
#      nt 2..8 and skips target bulges -----#
classify_seed_row <- function(m_matches, m_unmatches, t_matches, t_unmatches) {
  mm <- strsplit(m_matches,   "", fixed = TRUE)[[1]]
  mu <- strsplit(m_unmatches, "", fixed = TRUE)[[1]]
  tm <- strsplit(t_matches,   "", fixed = TRUE)[[1]]
  tu <- strsplit(t_unmatches, "", fixed = TRUE)[[1]]

  n <- max(length(mm), length(mu), length(tm), length(tu))
  pad <- function(v) c(v, rep(" ", n - length(v)))
  mm <- pad(mm); mu <- pad(mu); tm <- pad(tm); tu <- pad(tu)

  m_present <- (mm != " ") | (mu != " ")
  m_cols    <- which(m_present)
  if (length(m_cols) < 8) return("none")

  nt_cols <- rev(m_cols)             # nt_cols[k] = alignment col of miRNA nt k
  paired  <- function(k) mm[nt_cols[k]] != " " && tm[nt_cols[k]] != " "

  if (!all(vapply(2:7, paired, logical(1)))) return("none")
  pos8 <- paired(8)

  t_at_1 <- if (tm[nt_cols[1]] != " ") tm[nt_cols[1]] else tu[nt_cols[1]]
  a1     <- !is.na(t_at_1) && t_at_1 == "A"

  if (pos8 &&  a1) return("8mer")
  if (pos8)        return("7mer-m8")
  if (a1)          return("7mer-A1")
  "6mer"
}

#----- Canonical-seed class composition per miRNA -----#
plot_seed_class <- function() {
  seed_levels <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "none")
  data <- merged |>
    mutate(seed_class = factor(
      mapply(
        classify_seed_row,
        miRNA_matches, miRNA_unmatches, Target_matches, Target_unmatches,
        SIMPLIFY = TRUE, USE.NAMES = FALSE
      ),
      levels = seed_levels
    ))

  ggplot(data, aes(x = seed_class, y = -log10(P_value), fill = seed_class)) +
    geom_violin(alpha = 0.55, color = NA, scale = "width") +
    geom_boxplot(width = 0.18, alpha = 0.9, outlier.size = 0.3) +
    facet_wrap(~Calibration, labeller = facet_labels) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      x = "Seed class (canonical miRNA targeting)",
      y = "-log10(p-value)",
      title = "p-value distribution by seed class"
    ) +
    theme_bw(base_size = basesize) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}


#----- Plot registry: name -> builder -----#
plot_builders <- list(
  pvalue_distribution = plot_pvalue_distribution,
  position_pvalue     = plot_position_pvalue,
  position_energy     = plot_position_energy,
  volcano             = plot_volcano,
  pvalue_ecdf         = plot_pvalue_ecdf,
  calibration_delta   = plot_calibration_delta,
  per_mirna           = plot_per_mirna,
  position_density    = plot_position_density,
  seed_class          = plot_seed_class
)


#----- Resolve the requested selection -----#
# "all" -> every plot in registered order; otherwise comma-split + trimmed.
selected <- if (identical(plot_type, "all")) {
  names(plot_builders)
} else {
  trimws(strsplit(plot_type, ",", fixed = TRUE)[[1]])
}

unknown <- setdiff(selected, names(plot_builders))
if (length(unknown) > 0) {
  stop(sprintf(
    "Unknown plots.type entries: %s. Valid: %s, all, or a comma-separated subset.",
    paste(unknown, collapse = ", "),
    paste(names(plot_builders), collapse = ", ")
  ))
}

plots <- lapply(selected, function(name) plot_builders[[name]]())
# Builders may return NULL when their preconditions aren't met (e.g. paired
# calibration_delta in a single-variant run) — drop those silently.
plots <- Filter(Negate(is.null), plots)
n     <- length(plots)

if (n == 0) {
  stop("plots: no plots to render (all selected builders skipped — check rnacalibrate.calibration_variant and plots.type).")
}


#----- Compose & save -----#
# Equal column widths: one row up to n = 3, else 3 columns with plot_spacer() padding.
if (n == 1) {
  combined   <- plots[[1]]
  out_width  <- 9.5
  out_height <- 5
} else if (n <= 3) {
  combined   <- (Reduce(`|`, plots) + plot_layout(guides = "collect")) &
                  theme(legend.position = "bottom")
  out_width  <- if (n == 2) 12 else 18
  out_height <- 5
} else {
  n_rows <- ceiling(n / 3)
  pad    <- 3L * n_rows - n
  if (pad > 0L) {
    plots <- c(plots, replicate(pad, patchwork::plot_spacer(), simplify = FALSE))
  }
  # widths = rep(1, 3) forces equal column widths; without it patchwork
  # distributes width by intrinsic content (dual-facet plots end up wider
  # than single-facet ones) and the right columns visibly shrink.
  combined   <- wrap_plots(plots, ncol = 3, widths = rep(1, 3), guides = "collect") &
                  theme(legend.position = "bottom")
  out_width  <- 18
  out_height <- 5 * n_rows
}

dir.create(dirname(snakemake@output$pdf), showWarnings = FALSE, recursive = TRUE)
ggsave(
  snakemake@output$pdf, combined,
  width = out_width, height = out_height, units = "in"
)
