#!/usr/bin/env Rscript

# ============================================================
# Generar matrices de recuentos:
# 1. GUS x muestra
# 2. Clase GUS x muestra
#
# Reglas:
# - GUS:
#   * uniq hit -> contar por sseqid
#   * uncertain o uncertain+clase -> contar como Uncertain
#
# - Clase:
#   * uniq hit -> contar por gus_class
#   * uncertain+clase -> contar por gus_class
#   * uncertain -> contar como Uncertain
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        "Uso: Rscript count_GUS_matrices.R ",
        "samples.tsv results_directory output_directory"
    )
}

metadata_file <- args[1]
results_dir <- args[2]
output_dir <- args[3]

if (!file.exists(metadata_file)) {
    stop("No existe el archivo metadata: ", metadata_file)
}

if (!dir.exists(results_dir)) {
    stop("No existe el directorio de resultados: ", results_dir)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Leer identificadores de muestras
# ------------------------------------------------------------

metadata <- read.delim(
    metadata_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
)

if (ncol(metadata) < 2) {
    stop("samples.tsv debe contener al menos dos columnas.")
}

sample_ids <- as.character(metadata[[1]])

if (anyDuplicated(sample_ids)) {
    stop("Hay identificadores de muestra duplicados en samples.tsv.")
}

# Listas para almacenar los recuentos
gus_counts <- list()
class_counts <- list()

# ------------------------------------------------------------
# 2. Procesar cada muestra
# ------------------------------------------------------------

for (sample in sample_ids) {

    input_file <- file.path(
        results_dir,
        paste0(sample, ".GUS_filtered_best.tsv")
    )

    cat("\nProcesando:", sample, "\n")

    if (!file.exists(input_file)) {
        stop("No existe el archivo: ", input_file)
    }

    df <- read.delim(
        input_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    required <- c("sseqid", "gus_class", "hit_status")

    missing <- setdiff(required, colnames(df))

    if (length(missing) > 0) {
        stop(
            "Faltan columnas en ", input_file, ": ",
            paste(missing, collapse = ", ")
        )
    }

    # --------------------------------------------------------
    # 2A. Recuento por GUS
    # --------------------------------------------------------

    if (nrow(df) == 0) {

        gus_sample <- integer(0)

    } else {

        # Todas las asignaciones inciertas se agrupan
        # en una única categoría llamada Uncertain.
        gus_labels <- ifelse(
            grepl("^uncertain", df$hit_status, ignore.case = TRUE),
            "Uncertain",
            df$sseqid
        )

        gus_sample <- table(gus_labels)
    }

    gus_counts[[sample]] <- gus_sample

    # --------------------------------------------------------
    # 2B. Recuento por clase GUS
    # --------------------------------------------------------

    if (nrow(df) == 0) {

        class_sample <- integer(0)

    } else {

        # uncertain+clase conserva su clase estructural.
        # uncertain sin clase se agrupa como Uncertain.
        class_labels <- ifelse(
            tolower(df$hit_status) == "uncertain",
            "Uncertain",
            df$gus_class
        )

        # Evitar valores vacíos o NA
        class_labels[
            is.na(class_labels) | class_labels == ""
        ] <- "Unknown"

        class_sample <- table(class_labels)
    }

    class_counts[[sample]] <- class_sample

    cat("Lecturas asignadas:", nrow(df), "\n")
}

# ------------------------------------------------------------
# 3. Función para construir matrices
# ------------------------------------------------------------

build_matrix <- function(count_list, sample_ids) {

    categories <- sort(unique(unlist(
        lapply(count_list, names)
    )))

    # Colocar Uncertain al final, si existe
    categories <- categories[categories != "Uncertain"]

    if (any(vapply(
        count_list,
        function(x) "Uncertain" %in% names(x),
        logical(1)
    ))) {
        categories <- c(categories, "Uncertain")
    }

    if (length(categories) == 0) {
        categories <- character(0)
    }

    result <- matrix(
        0,
        nrow = length(categories),
        ncol = length(sample_ids),
        dimnames = list(categories, sample_ids)
    )

    for (sample in sample_ids) {

        counts <- count_list[[sample]]

        if (length(counts) > 0) {
            result[names(counts), sample] <- as.numeric(counts)
        }
    }

    result
}

# ------------------------------------------------------------
# 4. Construir las dos matrices
# ------------------------------------------------------------

gus_matrix <- build_matrix(gus_counts, sample_ids)

class_matrix <- build_matrix(class_counts, sample_ids)

# ------------------------------------------------------------
# 5. Guardar resultados
# ------------------------------------------------------------

gus_output <- file.path(
    output_dir,
    "GUS_x_sample_counts.tsv"
)

class_output <- file.path(
    output_dir,
    "GUS_class_x_sample_counts.tsv"
)

gus_df <- data.frame(
    GUS = rownames(gus_matrix),
    gus_matrix,
    check.names = FALSE
)

class_df <- data.frame(
    GUS_class = rownames(class_matrix),
    class_matrix,
    check.names = FALSE
)

write.table(
    gus_df,
    file = gus_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

write.table(
    class_df,
    file = class_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

cat("\n====================================\n")
cat("MATRICES GENERADAS CORRECTAMENTE\n")
cat("====================================\n")
cat("Muestras:", length(sample_ids), "\n")
cat("GUS:", nrow(gus_matrix), "\n")
cat("Clases:", nrow(class_matrix), "\n")
cat("Matriz GUS:", gus_output, "\n")
cat("Matriz clases:", class_output, "\n")
