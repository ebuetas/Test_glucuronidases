
#!/usr/bin/env Rscript

# ============================================================
# Filtrado de hits DIAMOND y selección de una GUS por lectura
#
# Filtros:
#   Query coverage >= 90%
#   Identity >= 50%
#   E-value <= 1e-5
#
# Selección:
#   Diferencia de bit score >= 5%: uniq hit
#   Hits ambiguos de la misma clase: selección aleatoria
#   Hits ambiguos de distintas clases: conservar el mejor
#
# Salidas:
#   1. *_filtered_hits.tsv: todos los hits que pasan filtros
#   2. *_filtered_best.tsv: una asignación por qseqid
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Uso: Rscript filter_select_GUS.R archivo_GUS.tsv")
}

input_file <- args[1]

if (!file.exists(input_file)) {
  stop("No existe el archivo: ", input_file)
}

# Semilla fija para reproducibilidad de las selecciones aleatorias
set.seed(20261117)

# Leer resultados previamente anotados
df <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Comprobar columnas necesarias
required_columns <- c(
  "qseqid",
  "sseqid",
  "pident",
  "evalue",
  "bitscore",
  "query_coverage",
  "gus_class"
)

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(
    "Faltan estas columnas en el archivo: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Asegurar que las variables numéricas son numéricas
numeric_columns <- c(
  "pident",
  "evalue",
  "bitscore",
  "query_coverage"
)

for (column in numeric_columns) {
  df[[column]] <- as.numeric(df[[column]])
}

# ------------------------------------------------------------
# 1. FILTRADO
# ------------------------------------------------------------

filtered <- df[
  !is.na(df$qseqid) &
  !is.na(df$pident) &
  !is.na(df$evalue) &
  !is.na(df$bitscore) &
  !is.na(df$query_coverage) &
  df$query_coverage >= 90 &
  df$pident >= 50 &
  df$evalue <= 1e-5,
]

# Guardar todos los hits que pasan los filtros
filtered_file <- sub(
  "\\.tsv$",
  "_filtered_hits.tsv",
  input_file,
  ignore.case = TRUE
)

if (identical(filtered_file, input_file)) {
  filtered_file <- paste0(input_file, "_filtered_hits.tsv")
}

write.table(
  filtered,
  file = filtered_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Hits iniciales:", nrow(df), "\n")
cat("Hits después del filtrado:", nrow(filtered), "\n")

# ------------------------------------------------------------
# 2. SELECCIÓN DE UN HIT POR qseqid
# ------------------------------------------------------------

if (nrow(filtered) == 0) {

  best_hits <- filtered
  best_hits$hit_status <- character(0)

} else {

  # Ordenar los hits por:
  # 1. Mayor bit score
  # 2. Menor E-value
  # 3. Mayor identidad
  # 4. ID de GUS para desempate reproducible

  filtered <- filtered[
    order(
      filtered$qseqid,
      -filtered$bitscore,
      filtered$evalue,
      -filtered$pident,
      filtered$sseqid
    ),
  ]

  # Agrupar las filas por lectura
  groups <- split(filtered, filtered$qseqid)

  selected_list <- lapply(groups, function(group) {

    # El primer hit es el de mayor bit score
    group <- group[
      order(
        -group$bitscore,
        group$evalue,
        -group$pident,
        group$sseqid
      ),
    ]

    best <- group[1, , drop = FALSE]

    # Si solo hay un hit que pasa los filtros
    if (nrow(group) == 1) {
      best$hit_status <- "uniq hit"
      return(best)
    }

    # Segundo mejor hit
    second_score <- group$bitscore[2]
    best_score <- group$bitscore[1]

    # Diferencia relativa entre primer y segundo hit
    score_difference <- (best_score - second_score) / second_score

    # Si el mejor supera al segundo en al menos un 5%
    if (score_difference >= 0.05) {

      best$hit_status <- "uniq hit"
      return(best)

    } else {

      # Hits dentro del margen del 5% respecto al mejor score
      near_best <- group[
        group$bitscore >= best_score / 1.05,
        ,
        drop = FALSE
      ]

      classes <- unique(near_best$gus_class)

      # Comprobar si todos pertenecen a la misma clase conocida
      same_known_class <-
        length(classes) == 1 &&
        !is.na(classes[1]) &&
        classes[1] != "" &&
        classes[1] != "Unknown"

      if (same_known_class) {

        # Elegir aleatoriamente entre los hits de esa clase
        chosen_index <- sample(
          seq_len(nrow(near_best)),
          size = 1
        )

        best <- near_best[chosen_index, , drop = FALSE]

        best$hit_status <- paste0(
          "uncertain+",
          as.character(best$gus_class)
        )

      } else {

        # Clases diferentes o clase desconocida:
        # conservar el hit con mayor bit score
        best$hit_status <- "uncertain"
      }

      return(best)
    }
  })

  best_hits <- do.call(rbind, selected_list)
  rownames(best_hits) <- NULL
}

# ------------------------------------------------------------
# 3. GUARDAR RESULTADOS FINALES
# ------------------------------------------------------------

best_file <- sub(
  "\\.tsv$",
  "_filtered_best.tsv",
  input_file,
  ignore.case = TRUE
)

if (identical(best_file, input_file)) {
  best_file <- paste0(input_file, "_filtered_best.tsv")
}

write.table(
  best_hits,
  file = best_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ------------------------------------------------------------
# 4. RESUMEN
# ------------------------------------------------------------

cat("\nANÁLISIS COMPLETADO\n")
cat("----------------------------\n")
cat("Lecturas con asignación:", nrow(best_hits), "\n")

if (nrow(best_hits) > 0) {
  print(table(best_hits$hit_status))
}

cat("\nArchivos generados:\n")
cat(filtered_file, "\n")
cat(best_file, "\n")
