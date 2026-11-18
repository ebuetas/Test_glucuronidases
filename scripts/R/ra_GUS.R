
#!/usr/bin/env Rscript

# ============================================================
# GUS relative abundance and statistical analysis
#
# Inputs:
#   1. GUS_class_x_sample_counts.tsv
#   2. QC_summary.tsv
#   3. clinical_metadata.tsv
#
# Normalization:
#   GUS class counts / Nonhuman_reads * 100
#
# Statistics:
#   Wilcoxon rank-sum test (independent groups)
#   Benjamini-Hochberg FDR correction
#
# Outputs:
#   - GUS_class_relative_abundance.tsv
#   - GUS_class_statistics.tsv
#   - GUS_class_heatmap.pdf
#   - GUS_class_boxplots.pdf
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    "Usage: Rscript relative_abundance_GUS.R ",
    "GUS_counts.tsv QC_summary.tsv clinical_metadata.tsv output_dir"
  )
}

counts_file <- args[1]
qc_file <- args[2]
clinical_file <- args[3]
output_dir <- args[4]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Read input files
# ------------------------------------------------------------

counts <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

qc <- read.delim(
  qc_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

clinical <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(counts) <- trimws(colnames(counts))
colnames(qc) <- trimws(colnames(qc))
colnames(clinical) <- trimws(colnames(clinical))

# ------------------------------------------------------------
# 2. Check required columns
# ------------------------------------------------------------

if (!"GUS" %in% colnames(counts)) {
  stop("No se encuentra la columna GUS en la matriz de recuentos.")
}

if (!all(c("Sample", "Nonhuman_reads") %in% colnames(qc))) {
  stop("QC_summary.tsv debe contener Sample y Nonhuman_reads.")
}

if (!all(c("ID", "Menopause_or_not") %in% colnames(clinical))) {
  stop("clinical_metadata.tsv debe contener ID y Menopause_or_not.")
}

# ------------------------------------------------------------
# 3. Prepare sample identifiers and check duplicates
# ------------------------------------------------------------

qc$Sample <- trimws(as.character(qc$Sample))
clinical$ID <- trimws(as.character(clinical$ID))

if (anyDuplicated(qc$Sample)) {
  stop("Hay muestras duplicadas en QC_summary.tsv.")
}

if (anyDuplicated(clinical$ID)) {
  stop("Hay muestras duplicadas en clinical_metadata.tsv.")
}

sample_ids <- colnames(counts)[-1]

if (anyDuplicated(sample_ids)) {
  stop("Hay nombres de muestra duplicados en la matriz de recuentos.")
}

# Check that all samples are present in both metadata files
missing_qc <- setdiff(sample_ids, qc$Sample)
missing_clinical <- setdiff(sample_ids, clinical$ID)

if (length(missing_qc) > 0) {
  stop(
    "Faltan estas muestras en QC_summary.tsv: ",
    paste(missing_qc, collapse = ", ")
  )
}

if (length(missing_clinical) > 0) {
  stop(
    "Faltan estas muestras en clinical_metadata.tsv: ",
    paste(missing_clinical, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 4. Normalize menopause status
# ------------------------------------------------------------
normalize_menopause <- function(x) {

  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^[:alnum:]]", "", x)

  pre_values <- c(
    "pre", "premenopausal", "premenopause",
    "premenopausica", "premenopausico",
    "no", "n", "0", "false",
    "notmenopausal", "nonmenopausal"
  )

  post_values <- c(
    "post", "postmenopausal", "postmenopause",
    "postmenopausica", "postmenopausico",
    "yes", "y", "1", "true",
    "menopausal", "si", "sí"
  )

  result <- rep(NA_character_, length(x))

  result[x %in% pre_values] <- "Pre"
  result[x %in% post_values] <- "Post"

  return(result)
}



clinical$Menopause_group <- normalize_menopause(
  clinical$Menopause_or_not
)

if (any(is.na(clinical$Menopause_group))) {

  unknown_values <- unique(
    clinical$Menopause_or_not[
      is.na(clinical$Menopause_group)
    ]
  )

  stop(
    "Hay valores de Menopause_or_not que no se pueden clasificar: ",
    paste(unknown_values, collapse = ", "),
    ". Añade estos valores a normalize_menopause()."
  )
}

# Check number of samples per group
group_table <- table(clinical$Menopause_group)

cat("\nMuestras por grupo:\n")
print(group_table)

if (!all(c("Pre", "Post") %in% names(group_table))) {
  stop("Se necesitan muestras tanto Pre como Post.")
}

# ------------------------------------------------------------
# 5. Extract non-human reads
# ------------------------------------------------------------

qc$Nonhuman_reads <- as.numeric(qc$Nonhuman_reads)

if (any(is.na(qc$Nonhuman_reads)) ||
    any(qc$Nonhuman_reads <= 0)) {
  stop("Nonhuman_reads contiene valores ausentes, no numéricos o <= 0.")
}

nonhuman <- qc$Nonhuman_reads[
  match(sample_ids, qc$Sample)
]

names(nonhuman) <- sample_ids

# ------------------------------------------------------------
# 6. Prepare count matrix
# ------------------------------------------------------------

count_matrix <- as.matrix(counts[, -1, drop = FALSE])

storage.mode(count_matrix) <- "numeric"

if (any(is.na(count_matrix)) || any(count_matrix < 0)) {
  stop("La matriz de recuentos contiene valores no válidos.")
}

rownames(count_matrix) <- counts$GUS

# ------------------------------------------------------------
# 7. Calculate relative abundance (%)
# ------------------------------------------------------------

relative_abundance <- sweep(
  count_matrix,
  MARGIN = 2,
  STATS = nonhuman,
  FUN = "/"
) * 100

# Save relative abundance matrix
relative_output <- data.frame(
  GUS = rownames(relative_abundance),
  relative_abundance,
  check.names = FALSE
)

write.table(
  relative_output,
  file = file.path(
    output_dir,
    "GUS_relative_abundance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Prepare clinical groups in matrix order
# ------------------------------------------------------------

group <- clinical$Menopause_group[
  match(sample_ids, clinical$ID)
]

names(group) <- sample_ids

# ------------------------------------------------------------
# 9. Statistical comparison: Pre vs Post
# ------------------------------------------------------------

p_values <- rep(NA_real_, nrow(relative_abundance))

statistics <- data.frame(
  GUS = rownames(relative_abundance),
  N_Pre = NA_integer_,
  N_Post = NA_integer_,
  Mean_Pre = NA_real_,
  Mean_Post = NA_real_,
  Median_Pre = NA_real_,
  Median_Post = NA_real_,
  Mean_difference_Post_minus_Pre = NA_real_,
  Fold_change_Post_vs_Pre = NA_real_,
  Log2_FC = NA_real_,
  p_value = NA_real_,
  stringsAsFactors = FALSE
)

pseudocount <- 1e-6

for (i in seq_len(nrow(relative_abundance))) {

  values <- relative_abundance[i, ]

  pre <- values[group == "Pre"]
  post <- values[group == "Post"]

  statistics$N_Pre[i] <- length(pre)
  statistics$N_Post[i] <- length(post)

  statistics$Mean_Pre[i] <- mean(pre)
  statistics$Mean_Post[i] <- mean(post)

  statistics$Median_Pre[i] <- median(pre)
  statistics$Median_Post[i] <- median(post)

  statistics$Mean_difference_Post_minus_Pre[i] <-
    mean(post) - mean(pre)

  statistics$Fold_change_Post_vs_Pre[i] <-
    (mean(post) + pseudocount) /
    (mean(pre) + pseudocount)

  statistics$Log2_FC[i] <-
    log2(statistics$Fold_change_Post_vs_Pre[i])

  # Wilcoxon rank-sum test, independent samples
  if (length(pre) >= 2 && length(post) >= 2) {

    test <- suppressWarnings(
      wilcox.test(
        post,
        pre,
        paired = FALSE,
        exact = FALSE
      )
    )

    p_values[i] <- test$p.value
  }
}

statistics$p_value <- p_values

# Multiple testing correction
statistics$p_adjusted <- p.adjust(
  statistics$p_value,
  method = "BH"
)

statistics$Significant_FDR_0.05 <- ifelse(
  is.na(statistics$p_adjusted),
  NA,
  statistics$p_adjusted < 0.05
)

# Save statistical results
write.table(
  statistics,
  file = file.path(
    output_dir,
    "GUS_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# 9B. GUS CLASS PREVALENCE: Pre vs Post
# ------------------------------------------------------------

# Minimum number of reads required to consider a GUS class present
prevalence_threshold <- 5

# Identify samples belonging to each group
pre_samples <- sample_ids[group == "Pre"]
post_samples <- sample_ids[group == "Post"]

n_pre <- length(pre_samples)
n_post <- length(post_samples)

# Initialize prevalence results
prevalence_results <- data.frame(
  GUS = rownames(count_matrix),

  N_Pre = n_pre,
  Present_Pre = NA_integer_,
  Absent_Pre = NA_integer_,
  Prevalence_Pre_pct = NA_real_,
  CI95_Pre_lower = NA_real_,
  CI95_Pre_upper = NA_real_,

  N_Post = n_post,
  Present_Post = NA_integer_,
  Absent_Post = NA_integer_,
  Prevalence_Post_pct = NA_real_,
  CI95_Post_lower = NA_real_,
  CI95_Post_upper = NA_real_,

  Prevalence_difference_Post_minus_Pre = NA_real_,
  p_value = NA_real_,
  stringsAsFactors = FALSE
)

# Calculate prevalence for each GUS class
for (i in seq_len(nrow(count_matrix))) {

  # Counts for each group
  pre_counts <- count_matrix[i, pre_samples]
  post_counts <- count_matrix[i, post_samples]

  # Presence/absence based on minimum read threshold
  pre_present <- pre_counts >= prevalence_threshold
  post_present <- post_counts >= prevalence_threshold

  # Number of positive and negative individuals
  n_pre_present <- sum(pre_present)
  n_pre_absent <- sum(!pre_present)

  n_post_present <- sum(post_present)
  n_post_absent <- sum(!post_present)

  # Prevalence percentages
  prev_pre <- 100 * n_pre_present / n_pre
  prev_post <- 100 * n_post_present / n_post

  # Exact binomial 95% confidence intervals
  ci_pre <- binom.test(
    n_pre_present,
    n_pre
  )$conf.int * 100

  ci_post <- binom.test(
    n_post_present,
    n_post
  )$conf.int * 100

  # Fisher's exact test
  contingency_table <- matrix(
    c(
      n_pre_present,
      n_pre_absent,
      n_post_present,
      n_post_absent
    ),
    nrow = 2,
    byrow = TRUE
  )

  fisher_result <- fisher.test(contingency_table)

  # Store results
  prevalence_results$Present_Pre[i] <- n_pre_present
  prevalence_results$Absent_Pre[i] <- n_pre_absent
  prevalence_results$Prevalence_Pre_pct[i] <- prev_pre
  prevalence_results$CI95_Pre_lower[i] <- ci_pre[1]
  prevalence_results$CI95_Pre_upper[i] <- ci_pre[2]

  prevalence_results$Present_Post[i] <- n_post_present
  prevalence_results$Absent_Post[i] <- n_post_absent
  prevalence_results$Prevalence_Post_pct[i] <- prev_post
  prevalence_results$CI95_Post_lower[i] <- ci_post[1]
  prevalence_results$CI95_Post_upper[i] <- ci_post[2]

  prevalence_results$Prevalence_difference_Post_minus_Pre[i] <-
    prev_post - prev_pre

  prevalence_results$p_value[i] <- fisher_result$p.value
}

# Multiple testing correction
prevalence_results$p_adjusted <- p.adjust(
  prevalence_results$p_value,
  method = "BH"
)

prevalence_results$Significant_FDR_0.05 <- (
  prevalence_results$p_adjusted < 0.05
)

# Save prevalence statistics
write.table(
  prevalence_results,
  file = file.path(
    output_dir,
    "GUS_prevalence_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nPrevalence analysis completed.\n")
cat("Minimum reads for presence:", prevalence_threshold, "\n")

# ------------------------------------------------------------
# 9C. Prevalence figure
# ------------------------------------------------------------

prev <- nrow(prevalence_results)

# PDF vectorial output
cairo_pdf(
  filename = file.path(
    output_dir,
    "GUS_prevalence.pdf"
  ),
  width = 11,
  height = max(6, prev * 0.45),
  family = "Arial",
  bg = "white"
)

par(
  mar = c(5, 11, 4, 3),
  xpd = FALSE
)

# Y-axis positions
y_positions <- seq_len(prev)

# Create empty plot
plot(
  NA,
  xlim = c(0, 100),
  ylim = c(0.5, prev + 0.5),
  xlab = "Prevalence (%)",
  ylab = "",
  yaxt = "n",
  main = paste0(
    "GUS prevalence (",
    prevalence_threshold,
    " or more reads)"
  ),
  xaxs = "i",
  yaxs = "i"
)

# Add horizontal grid
abline(
  v = seq(0, 100, by = 20),
  col = "grey90",
  lty = 1
)

# Class names on Y-axis
class_labels <- prevalence_results$GUS

# Mark significant classes with an asterisk
class_labels[
  prevalence_results$Significant_FDR_0.05
] <- paste0(
  class_labels[
    prevalence_results$Significant_FDR_0.05
  ],
  " *"
)

axis(
  side = 2,
  at = y_positions,
  labels = class_labels,
  las = 1,
  cex.axis = 0.8
)

# Plot Pre and Post prevalence with 95% confidence intervals
for (i in seq_len(prev)) {

  # Pre group
  y_pre <- y_positions[i] + 0.12

  segments(
    x0 = prevalence_results$CI95_Pre_lower[i],
    y0 = y_pre,
    x1 = prevalence_results$CI95_Pre_upper[i],
    y1 = y_pre,
    col = "steelblue",
    lwd = 2
  )

  points(
    prevalence_results$Prevalence_Pre_pct[i],
    y_pre,
    pch = 16,
    col = "steelblue",
    cex = 1.2
  )

  # Post group
  y_post <- y_positions[i] - 0.12

  segments(
    x0 = prevalence_results$CI95_Post_lower[i],
    y0 = y_post,
    x1 = prevalence_results$CI95_Post_upper[i],
    y1 = y_post,
    col = "tomato",
    lwd = 2
  )

  points(
    prevalence_results$Prevalence_Post_pct[i],
    y_post,
    pch = 16,
    col = "tomato",
    cex = 1.2
  )
}

# Legend
legend(
  "bottomright",
  legend = c("Premenopausal", "Postmenopausal"),
  col = c("steelblue", "tomato"),
  pch = 16,
  lwd = 2,
  bty = "n",
  cex = 0.9
)

# Significance explanation
mtext(
  "* FDR < 0.05",
  side = 1,
  line = 3.5,
  adj = 1,
  cex = 0.8
)

dev.off()



# ------------------------------------------------------------
# 10. Heatmap
# ------------------------------------------------------------

# Order samples by menopause group
sample_order <- c(
  sample_ids[group == "Pre"],
  sample_ids[group == "Post"]
)

heatmap_matrix <- relative_abundance[
  ,
  sample_order,
  drop = FALSE
]

# Log transformation for visualization
log_matrix <- log10(heatmap_matrix + pseudocount)

# Row-wise Z-score
z_matrix <- t(scale(t(log_matrix)))

# Replace undefined Z-scores with zero
z_matrix[!is.finite(z_matrix)] <- 0

# Limit extreme Z-scores
z_matrix[z_matrix > 2.5] <- 2.5
z_matrix[z_matrix < -2.5] <- -2.5

# Sample group colors
ordered_groups <- group[sample_order]

group_colors <- ifelse(
  ordered_groups == "Pre",
  "steelblue",
  "tomato"
)

# Heatmap colors
heat_colors <- colorRampPalette(
  c("navy", "white", "firebrick3")
)(100)

# Create PDF


cairo_pdf(
  filename = file.path(
    output_dir,
    "GUS_heatmap.pdf"
  ),
  width = max(12, length(sample_order) * 0.30),
  height = 7,
  family = "Arial",
  bg = "white"
)
# Reserve space on the right for the color scale
par(mar = c(11, 9, 4, 9))

# Draw heatmap
image(
  x = seq_len(ncol(z_matrix)),
  y = seq_len(nrow(z_matrix)),
  z = t(z_matrix),
  col = heat_colors,
  zlim = c(-2.5, 2.5),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  main = "GUS relative abundance"
)

# X-axis labels: one color per sample
for (i in seq_len(ncol(z_matrix))) {
  axis(
    side = 1,
    at = i,
    labels = colnames(z_matrix)[i],
    las = 2,
    cex.axis = 0.65,
    col.axis = group_colors[i],
    tick = FALSE
  )
}

# Y-axis labels: GUS classes
axis(
  side = 2,
  at = seq_len(nrow(z_matrix)),
  labels = rownames(z_matrix),
  las = 2,
  cex.axis = 0.8
)

box()

# Menopause group legend
legend(
  "topright",
  inset = c(-0.25, 0),
  legend = c("Premenopausal", "Postmenopausal"),
  fill = c("steelblue", "tomato"),
  xpd = TRUE,
  bty = "n",
  cex = 0.8
)

# ------------------------------------------------------------
# Color scale legend (Z-score)
# ------------------------------------------------------------

# Create a separate plotting region on the right
par(
  fig = c(0.86, 0.91, 0.20, 0.80),
  mar = c(0, 0, 0, 0),
  new = TRUE
)

gradient_values <- seq(-2.5, 2.5, length.out = 100)

image(
  x = 1,
  y = gradient_values,
  z = matrix(gradient_values, nrow = 1),
  col = heat_colors,
  zlim = c(-2.5, 2.5),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  useRaster = TRUE
)

axis(
  side = 4,
  at = c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5),
  labels = c("-2.5", "-1.5", "-0.5", "0.5", "1.5", "2.5"),
  las = 1,
  cex.axis = 0.8
)

mtext(
  "Row Z-score",
  side = 4,
  line = 2.5,
  cex = 0.8
)

dev.off()


# ------------------------------------------------------------
# 11. Boxplots of individual GUS
# ------------------------------------------------------------

# Options:
# "all"       = plot all GUS
# "significant" = plot only GUS with FDR below threshold

boxplot_mode <- "all"

# Significance threshold
boxplot_fdr_threshold <- 0.05

# Number of plots per page
plots_per_row <- 3
plots_per_column <- 4
plots_per_page <- plots_per_row * plots_per_column

# Select GUS to plot
if (boxplot_mode == "significant") {

  selected_classes <- which(
    !is.na(statistics$p_adjusted) &
    statistics$p_adjusted < boxplot_fdr_threshold
  )

} else {

  selected_classes <- seq_len(nrow(relative_abundance))

}

# Fallback if no significant GUS
if (length(selected_classes) == 0) {

  warning(
    "No GUS passed the FDR threshold. ",
    "All GUS will be plotted instead."
  )

  selected_classes <- seq_len(nrow(relative_abundance))

}

# Number of selected GUS
number_selected <- length(selected_classes)

# Number of PDF pages
number_pages <- ceiling(
  number_selected / plots_per_page
)

cat("\nGenerating GUS boxplots...\n")
cat("GUS selected:", number_selected, "\n")
cat("PDF pages:", number_pages, "\n")

# Open multipage PDF
pdf(
  file.path(
    output_dir,
    "GUS_boxplots.pdf"
  ),
  width = 12,
  height = 10,
  onefile = TRUE
)

# Generate pages
for (page in seq_len(number_pages)) {

  start_index <- (page - 1) * plots_per_page + 1

  end_index <- min(
    page * plots_per_page,
    number_selected
  )

  current_indices <- selected_classes[
    start_index:end_index
  ]

  # Layout for each page
  par(
    mfrow = c(plots_per_column, plots_per_row),
    mar = c(4, 4, 4, 1),
    oma = c(1, 1, 3, 1)
  )

  for (i in current_indices) {

    values <- relative_abundance[i, ]

    plot_data <- list(
      Pre = values[group == "Pre"],
      Post = values[group == "Post"]
    )

    boxplot(
      plot_data,
      main = rownames(relative_abundance)[i],
      ylab = "Relative abundance (%)",
      xlab = "",
      col = c("steelblue", "tomato"),
      outline = FALSE,
      cex.main = 0.8
    )

    stripchart(
      plot_data,
      vertical = TRUE,
      method = "jitter",
      add = TRUE,
      pch = 16,
      cex = 0.5
    )

    # Add FDR
    current_fdr <- statistics$p_adjusted[i]

    if (!is.na(current_fdr)) {

      mtext(
        paste0(
          "FDR = ",
          signif(current_fdr, 3)
        ),
        side = 3,
        line = 0.2,
        cex = 0.7
      )

    }

  }

  # Page title
  mtext(
    paste0(
      "GUS abundance: page ",
      page,
      " of ",
      number_pages
    ),
    outer = TRUE,
    side = 3,
    line = 1,
    cex = 1.1,
    font = 2
  )

}

dev.off()

cat(
  "\nBoxplots saved:",
  file.path(output_dir, "GUS_boxplots.pdf"),
  "\n"
)



# ------------------------------------------------------------
# 12. Summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ANALYSIS COMPLETED\n")
cat("========================================\n")

cat("Samples:", length(sample_ids), "\n")
cat("Premenopausal:", sum(group == "Pre"), "\n")
cat("Postmenopausal:", sum(group == "Post"), "\n")
cat("GUS:", nrow(relative_abundance), "\n")

cat(
  "Significant classes (FDR < 0.05): ",
  sum(statistics$p_adjusted < 0.05, na.rm = TRUE),
  "\n"
)

cat("\nOutput directory:\n", output_dir, "\n")
