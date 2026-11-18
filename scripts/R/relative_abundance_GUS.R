
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

if (!"GUS_class" %in% colnames(counts)) {
  stop("No se encuentra la columna GUS_class en la matriz de recuentos.")
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

rownames(count_matrix) <- counts$GUS_class

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
  GUS_class = rownames(relative_abundance),
  relative_abundance,
  check.names = FALSE
)

write.table(
  relative_output,
  file = file.path(
    output_dir,
    "GUS_class_relative_abundance.tsv"
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
  GUS_class = rownames(relative_abundance),
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
    "GUS_class_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



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
    "GUS_class_heatmap.pdf"
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
  main = "GUS class relative abundance"
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
# 11. Boxplots
# ------------------------------------------------------------
cairo_pdf(
  filename = file.path(
    output_dir,
    "GUS_class_boxplots.pdf"
  ),
  width = 12, 
  height = 9,
  family = "Arial",
  bg = "white"
)


number_classes <- nrow(relative_abundance)

ncol_plot <- 3
nrow_plot <- ceiling(number_classes / ncol_plot)

par(
  mfrow = c(nrow_plot, ncol_plot),
  mar = c(5, 4, 3, 1)
)

for (i in seq_len(number_classes)) {

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
    outline = FALSE
  )

  stripchart(
    plot_data,
    vertical = TRUE,
    method = "jitter",
    add = TRUE,
    pch = 16,
    cex = 0.65
  )

  if (!is.na(statistics$p_adjusted[i])) {
    mtext(
      paste0(
        "FDR = ",
        signif(statistics$p_adjusted[i], 3)
      ),
      side = 3,
      line = 0.2,
      cex = 0.75
    )
  }
}

dev.off()

# ------------------------------------------------------------
# 12. Summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ANALYSIS COMPLETED\n")
cat("========================================\n")

cat("Samples:", length(sample_ids), "\n")
cat("Premenopausal:", sum(group == "Pre"), "\n")
cat("Postmenopausal:", sum(group == "Post"), "\n")
cat("GUS classes:", nrow(relative_abundance), "\n")

cat(
  "Significant classes (FDR < 0.05): ",
  sum(statistics$p_adjusted < 0.05, na.rm = TRUE),
  "\n"
)

cat("\nOutput directory:\n", output_dir, "\n")
