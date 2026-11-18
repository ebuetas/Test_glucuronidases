
#!/usr/bin/env Rscript

# ============================================================
# REFINE GUS FIGURES
#
# Inputs:
# 1. GUS_prevalence_statistics.tsv
# 2. GUS_relative_abundance.tsv
# 3. GUS_statistics.tsv
# 4. clinical_metadata.tsv
# 5. GUS_metadata.tsv
# 6. output_dir
#
# Prevalence:
#   Include GUS with unadjusted p-value < 0.1
#   Asterisk for unadjusted p-value < 0.05
#   Display Pre/Post prevalence and 95% CI
#
# Boxplots:
#   Include GUS with unadjusted p-value < 0.05
#   Display raw p-value and FDR
#   Show GUS ID and structural class
#   12 plots per page
#
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop(
    "Usage: Rscript refine_GUS_figures.R ",
    "prevalence_statistics.tsv ",
    "relative_abundance.tsv ",
    "abundance_statistics.tsv ",
    "clinical_metadata.tsv ",
    "GUS_metadata.tsv ",
    "output_dir"
  )
}

prevalence_file <- args[1]
abundance_file <- args[2]
statistics_file <- args[3]
clinical_file <- args[4]
gus_metadata_file <- args[5]
output_dir <- args[6]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Functions
# ------------------------------------------------------------

normalize_menopause <- function(x) {

  x <- tolower(trimws(as.character(x)))
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
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
    "menopausal", "si"
  )

  result <- rep(NA_character_, length(x))

  result[x %in% pre_values] <- "Pre"
  result[x %in% post_values] <- "Post"

  result
}

# ------------------------------------------------------------
# 2. Read input files
# ------------------------------------------------------------

prevalence <- read.delim(
  prevalence_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

abundance <- read.delim(
  abundance_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

statistics <- read.delim(
  statistics_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

clinical <- read.delim(
  clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gus_metadata <- read.delim(
  gus_metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove accidental spaces in column names
colnames(prevalence) <- trimws(colnames(prevalence))
colnames(abundance) <- trimws(colnames(abundance))
colnames(statistics) <- trimws(colnames(statistics))
colnames(clinical) <- trimws(colnames(clinical))
colnames(gus_metadata) <- trimws(colnames(gus_metadata))

# ------------------------------------------------------------
# 3. Check required columns
# ------------------------------------------------------------

required_prev <- c(
  "GUS",
  "Prevalence_Pre_pct",
  "Prevalence_Post_pct",
  "CI95_Pre_lower",
  "CI95_Pre_upper",
  "CI95_Post_lower",
  "CI95_Post_upper",
  "p_value"
)

required_abundance <- "GUS"

required_stats <- c(
  "GUS",
  "p_value",
  "p_adjusted"
)

required_clinical <- c(
  "ID",
  "Menopause_or_not"
)

required_gus_metadata <- c(
  "new_id",
  "gus_type"
)

if (!all(required_prev %in% colnames(prevalence))) {
  stop("Missing required columns in prevalence statistics.")
}

if (!all(required_abundance %in% colnames(abundance))) {
  stop("Missing GUS column in relative abundance file.")
}

if (!all(required_stats %in% colnames(statistics))) {
  stop("Missing required columns in abundance statistics.")
}

if (!all(required_clinical %in% colnames(clinical))) {
  stop("Clinical metadata must contain ID and Menopause_or_not.")
}

if (!all(required_gus_metadata %in% colnames(gus_metadata))) {
  stop("GUS_metadata.tsv must contain new_ID and gus_type.")
}

# ------------------------------------------------------------
# 4. Prepare GUS class mapping
# ------------------------------------------------------------

gus_metadata$new_id <- trimws(
  as.character(gus_metadata$new_id)
)

gus_metadata$gus_type <- trimws(
  as.character(gus_metadata$gus_type)
)

if (anyDuplicated(gus_metadata$new_id)) {
  stop("Duplicated GUS identifiers in GUS_metadata.tsv.")
}

class_mapping <- setNames(
  gus_metadata$gus_type,
  gus_metadata$new_id
)

# ------------------------------------------------------------
# 5. Prepare clinical metadata
# ------------------------------------------------------------

clinical$ID <- trimws(
  as.character(clinical$ID)
)

if (anyDuplicated(clinical$ID)) {
  stop("Duplicated sample IDs in clinical metadata.")
}

clinical$Menopause_group <- normalize_menopause(
  clinical$Menopause_or_not
)

if (anyNA(clinical$Menopause_group)) {
  stop(
    "Some menopause values could not be classified. ",
    "Check Menopause_or_not in clinical_metadata.tsv."
  )
}

# ------------------------------------------------------------
# 6. Convert numerical columns
# ------------------------------------------------------------

numeric_prev <- c(
  "Prevalence_Pre_pct",
  "Prevalence_Post_pct",
  "CI95_Pre_lower",
  "CI95_Pre_upper",
  "CI95_Post_lower",
  "CI95_Post_upper",
  "p_value"
)

for (col in numeric_prev) {
  prevalence[[col]] <- as.numeric(prevalence[[col]])
}

statistics$p_value <- as.numeric(statistics$p_value)
statistics$p_adjusted <- as.numeric(statistics$p_adjusted)

# Standardize GUS identifiers
prevalence$GUS <- trimws(as.character(prevalence$GUS))
abundance$GUS <- trimws(as.character(abundance$GUS))
statistics$GUS <- trimws(as.character(statistics$GUS))

# Check duplicated GUS identifiers
if (anyDuplicated(abundance$GUS)) {
  stop("Duplicated GUS identifiers in abundance matrix.")
}

if (anyDuplicated(statistics$GUS)) {
  stop("Duplicated GUS identifiers in abundance statistics.")
}

# ------------------------------------------------------------
# 7. PREVALENCE FIGURE
# ------------------------------------------------------------

selected_prev <- prevalence[
  !is.na(prevalence$p_value) &
    prevalence$p_value < 0.1,
  ,
  drop = FALSE
]

selected_prev <- selected_prev[
  order(selected_prev$p_value),
  ,
  drop = FALSE
]

cat(
  "\nGUS selected for prevalence figure:",
  nrow(selected_prev),
  "\n"
)

prevalence_pdf <- file.path(
  output_dir,
  "GUS_prevalence_pvalue_filtered.pdf"
)

cairo_pdf(
  filename = prevalence_pdf,
  width = 11,
  height = max(6, nrow(selected_prev) * 0.40 + 2),
  family = "sans",
  bg = "white"
)

if (nrow(selected_prev) == 0) {

  plot.new()

  text(
    0.5, 0.5,
    "No GUS with unadjusted p-value < 0.1",
    cex = 1.2
  )

} else {

  gus_class_labels <- unname(
    class_mapping[selected_prev$GUS]
  )

  gus_class_labels[is.na(gus_class_labels)] <- "Unknown"

  class_labels <- paste0(
    selected_prev$GUS,
    " (",
    gus_class_labels,
    ")"
  )

  significant <- !is.na(selected_prev$p_value) &
    selected_prev$p_value < 0.05

  class_labels[significant] <- paste0(
    class_labels[significant],
    " *"
  )

  n_prev <- nrow(selected_prev)
  y_positions <- seq_len(n_prev)

  par(mar = c(5, 11, 4, 3))

  plot(
    NA,
    xlim = c(0, 100),
    ylim = c(0.5, n_prev + 0.5),
    xlab = "Prevalence (%)",
    ylab = "",
    yaxt = "n",
    main = "GUS prevalence: Pre vs Post",
    xaxs = "i",
    yaxs = "i"
  )

  abline(
    v = seq(0, 100, by = 20),
    col = "grey90"
  )

  axis(
    side = 2,
    at = y_positions,
    labels = class_labels,
    las = 1,
    cex.axis = 0.8
  )

  for (i in seq_len(n_prev)) {

    # Premenopausal
    y_pre <- y_positions[i] + 0.12

    segments(
      x0 = selected_prev$CI95_Pre_lower[i],
      y0 = y_pre,
      x1 = selected_prev$CI95_Pre_upper[i],
      y1 = y_pre,
      col = "steelblue",
      lwd = 2
    )

    points(
      selected_prev$Prevalence_Pre_pct[i],
      y_pre,
      pch = 16,
      col = "steelblue",
      cex = 1.2
    )

    # Postmenopausal
    y_post <- y_positions[i] - 0.12

    segments(
      x0 = selected_prev$CI95_Post_lower[i],
      y0 = y_post,
      x1 = selected_prev$CI95_Post_upper[i],
      y1 = y_post,
      col = "tomato",
      lwd = 2
    )

    points(
      selected_prev$Prevalence_Post_pct[i],
      y_post,
      pch = 16,
      col = "tomato",
      cex = 1.2
    )
  }

  legend(
    "bottomright",
    legend = c("Premenopausal", "Postmenopausal"),
    col = c("steelblue", "tomato"),
    pch = 16,
    lwd = 2,
    bty = "n",
    cex = 0.9
  )

  mtext(
    "* Unadjusted p-value < 0.05",
    side = 1,
    line = 3.5,
    adj = 1,
    cex = 0.8
  )
}

dev.off()

# ------------------------------------------------------------
# 8. Select GUS for boxplots
# ------------------------------------------------------------

selected_stats <- statistics[
  !is.na(statistics$p_value) &
    statistics$p_value < 0.05,
  ,
  drop = FALSE
]

selected_stats <- selected_stats[
  order(selected_stats$p_value),
  ,
  drop = FALSE
]

# Keep the order determined by raw p-value
selected_gus <- selected_stats$GUS[
  selected_stats$GUS %in% abundance$GUS
]

if (length(selected_gus) != nrow(selected_stats)) {
  warning(
    "Some significant GUS were not found in the abundance matrix."
  )
}

# Match statistics to the selected GUS
selected_stats <- selected_stats[
  match(selected_gus, selected_stats$GUS),
  ,
  drop = FALSE
]

cat(
  "\nGUS selected for boxplots:",
  length(selected_gus),
  "\n"
)

# ------------------------------------------------------------
# 9. Prepare abundance matrix
# ------------------------------------------------------------

sample_ids <- colnames(abundance)[-1]

if (anyDuplicated(sample_ids)) {
  stop("Duplicated sample names in abundance matrix.")
}

abundance_matrix <- as.matrix(
  abundance[
    match(selected_gus, abundance$GUS),
    sample_ids,
    drop = FALSE
  ]
)

storage.mode(abundance_matrix) <- "numeric"

rownames(abundance_matrix) <- selected_gus

if (anyNA(abundance_matrix) ||
    any(!is.finite(abundance_matrix))) {
  stop("Invalid values in the relative abundance matrix.")
}

# ------------------------------------------------------------
# 10. Assign menopause groups to samples
# ------------------------------------------------------------

group <- clinical$Menopause_group[
  match(sample_ids, clinical$ID)
]

if (anyNA(group)) {
  missing_samples <- sample_ids[is.na(group)]

  stop(
    "Some samples are missing from clinical metadata: ",
    paste(missing_samples, collapse = ", ")
  )
}

if (!all(c("Pre", "Post") %in% unique(group))) {
  stop("Both Pre and Post groups are required.")
}

# ------------------------------------------------------------
# 11. Generate multipage boxplots
# ------------------------------------------------------------

plots_per_row <- 3
plots_per_column <- 4
plots_per_page <- plots_per_row * plots_per_column

number_selected <- length(selected_gus)

boxplot_pdf <- file.path(
  output_dir,
  "GUS_boxplots_pvalue_filtered.pdf"
)

cairo_pdf(
  filename = boxplot_pdf,
  width = 12,
  height = 10,
  family = "sans",
  bg = "white"
)

if (number_selected == 0) {

  plot.new()

  text(
    0.5, 0.5,
    "No GUS with unadjusted p-value < 0.05",
    cex = 1.2
  )

} else {

  number_pages <- ceiling(
    number_selected / plots_per_page
  )

  for (page in seq_len(number_pages)) {

    start_index <- (page - 1) * plots_per_page + 1

    end_index <- min(
      page * plots_per_page,
      number_selected
    )

    current_indices <- start_index:end_index

    par(
      mfrow = c(plots_per_column, plots_per_row),
      mar = c(4, 4, 5, 1),
      oma = c(1, 1, 3, 1)
    )

    for (i in current_indices) {

      gus_name <- selected_gus[i]

      gus_class <- unname(class_mapping[gus_name])

      if (is.na(gus_class) || length(gus_class) == 0) {
        gus_class <- "Unknown"
      }

      values <- abundance_matrix[i, ]

      plot_data <- list(
        Pre = values[group == "Pre"],
        Post = values[group == "Post"]
      )

      raw_p <- selected_stats$p_value[i]
      fdr <- selected_stats$p_adjusted[i]

      boxplot(
        plot_data,
        main = paste0(
          gus_name,
          " (",
          gus_class,
          ")\n",
          "p = ",
          format.pval(raw_p, digits = 3, eps = 0.001),
          " | FDR = ",
          format.pval(fdr, digits = 3, eps = 0.001)
        ),
        ylab = "Relative abundance (%)",
        xlab = "",
        col = c("steelblue", "tomato"),
        outline = FALSE,
        cex.main = 0.70
      )

      stripchart(
        plot_data,
        vertical = TRUE,
        method = "jitter",
        add = TRUE,
        pch = 16,
        cex = 0.5
      )
    }

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
}

dev.off()

# ------------------------------------------------------------
# 12. Finish
# ------------------------------------------------------------

cat("\nRefined figures completed.\n")
cat("Prevalence PDF:", prevalence_pdf, "\n")
cat("Boxplots PDF:", boxplot_pdf, "\n")
