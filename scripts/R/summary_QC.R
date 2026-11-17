# ============================================================
# RESUMEN DE QC: TRIMMOMATIC + BOWTIE2
# ============================================================

library(dplyr)
library(stringr)
library(readr)
library(tibble)

# ============================================================================
# 1. CONFIGURACION
# ============================================================================

PROJECT_DIR <- "/data/home/molivares/elena/glucuronidases"

METADATA <- file.path(
  PROJECT_DIR,
  "data",
  "samples.tsv"
)

LOG_DIR <- file.path(
  PROJECT_DIR,
  "logs"
)

RESULTS_DIR <- file.path(
  PROJECT_DIR,"data",
  "processed"
)

QC_RAW <- file.path(RESULTS_DIR, "fastqc_raw")
QC_TRIMMED <- file.path(RESULTS_DIR, "fastqc_trimmed")
QC_FINAL <- file.path(RESULTS_DIR, "fastqc_nonhuman")

OUTPUT_FILE <- file.path(
 PROJECT_DIR,"results",
  "QC_summary.tsv"
)


# ============================================================================
# 2. COMPROBAR DIRECTORIOS Y METADATOS
# ============================================================================

if (!file.exists(METADATA)) {
  stop("No se encuentra el archivo samples.tsv")
}

if (!dir.exists(LOG_DIR)) {
  stop("No se encuentra el directorio de logs")
}

# Leer archivo de muestras: dos columnas, sin encabezado
samples <- read.delim(
  METADATA,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

colnames(samples) <- c("Sample", "FASTQ")

samples$Sample <- trimws(samples$Sample)
samples$FASTQ <- trimws(samples$FASTQ)

if (anyDuplicated(samples$Sample)) {
  stop("Hay identificadores de muestra duplicados en samples.tsv")
}

cat("Muestras encontradas:", nrow(samples), "\n")


# ============================================================================
# 3. FUNCION PARA EXTRAER RECUENTOS DE TRIMMOMATIC
# ============================================================================

read_trimmomatic <- function(log_file) {

  if (!file.exists(log_file)) {
    warning("Falta el log de Trimmomatic: ", log_file)
    return(NULL)
  }

  lines <- readLines(log_file, warn = FALSE)

  extract_count <- function(pattern) {

    line <- grep(
      pattern,
      lines,
      value = TRUE
    )

    if (length(line) == 0) {
      return(NA_real_)
    }

    match <- regexec(
      "([0-9,]+)",
      line[1]
    )

    value <- regmatches(
      line[1],
      match
    )[[1]]

    if (length(value) < 2) {
      return(NA_real_)
    }

    as.numeric(gsub(",", "", value[2]))
  }

  data.frame(
    Raw_reads = extract_count("^Input Reads:"),
    Reads_after_Trimmomatic = extract_count("^Surviving:"),
    Reads_removed_Trimmomatic = extract_count("^Dropped:")
  )
}


# ============================================================================
# 4. FUNCION PARA EXTRAER RECUENTOS DE BOWTIE2
# ============================================================================

read_bowtie2 <- function(log_file) {

  if (!file.exists(log_file)) {
    warning("Falta el log de Bowtie2: ", log_file)
    return(NULL)
  }

  lines <- readLines(log_file, warn = FALSE)

  # Total de lecturas analizadas
  total_line <- grep(
    "reads; of these:",
    lines,
    value = TRUE
  )

  # Lecturas que no se alinearon: NO humanas
  nonhuman_line <- grep(
    "aligned 0 times",
    lines,
    value = TRUE
  )

  if (length(total_line) == 0 ||
      length(nonhuman_line) == 0) {

    warning("No se pudieron interpretar los datos de: ", log_file)
    return(NULL)
  }

  extract_number <- function(line) {

    match <- regexec(
      "^\\s*([0-9,]+)",
      line
    )

    value <- regmatches(
      line,
      match
    )[[1]]

    if (length(value) < 2) {
      return(NA_real_)
    }

    as.numeric(gsub(",", "", value[2]))
  }

  total <- extract_number(total_line[1])

  nonhuman <- extract_number(nonhuman_line[1])

  human <- total - nonhuman

  data.frame(
    Bowtie2_reads = total,
    Human_reads = human,
    Nonhuman_reads = nonhuman
  )
}


# ============================================================================
# 5. FUNCION PARA COMPROBAR INFORMES FASTQC
# ============================================================================

fastqc_report_exists <- function(directory, sample) {

  if (!dir.exists(directory)) {
    return(FALSE)
  }

  files <- list.files(
    directory,
    pattern = "_fastqc\\.html$"
  )

  any(
    startsWith(files, sample)
  )
}


# ============================================================================
# 6. PROCESAR TODAS LAS MUESTRAS
# ============================================================================

results <- vector("list", nrow(samples))

for (i in seq_len(nrow(samples))) {

  sample <- samples$Sample[i]

  cat("Procesando:", sample, "\n")

  # Archivos de log
  trim_log <- file.path(
    LOG_DIR,
    paste0(sample, "_trimmomatic.log")
  )

  bowtie_log <- file.path(
    LOG_DIR,
    paste0(sample, "_bowtie2.log")
  )

  # Leer resultados
  trim_data <- read_trimmomatic(trim_log)

  bowtie_data <- read_bowtie2(bowtie_log)

  # Crear fila aunque falte algún log
  if (is.null(trim_data)) {
    trim_data <- data.frame(
      Raw_reads = NA_real_,
      Reads_after_Trimmomatic = NA_real_,
      Reads_removed_Trimmomatic = NA_real_
    )
  }

  if (is.null(bowtie_data)) {
    bowtie_data <- data.frame(
      Bowtie2_reads = NA_real_,
      Human_reads = NA_real_,
      Nonhuman_reads = NA_real_
    )
  }

  # Comprobar informes HTML
  raw_report <- fastqc_report_exists(
    QC_RAW,
    sample
  )

  trimmed_report <- fastqc_report_exists(
    QC_TRIMMED,
    sample
  )

  final_report <- fastqc_report_exists(
    QC_FINAL,
    sample
  )

  # Combinar información
  results[[i]] <- cbind(
    data.frame(Sample = sample),
    trim_data,
    bowtie_data,
    data.frame(
      FastQC_raw = raw_report,
      FastQC_trimmed = trimmed_report,
      FastQC_nonhuman = final_report
    )
  )
}


# ============================================================================
# 7. UNIR RESULTADOS
# ============================================================================

QC <- do.call(rbind, results)


# ============================================================================
# 8. CALCULAR PORCENTAJES
# ============================================================================

QC$Reads_removed_total <- with(
  QC,
  Raw_reads - Nonhuman_reads
)

QC$Trimmomatic_retention_pct <- with(
  QC,
  ifelse(
    !is.na(Raw_reads) & Raw_reads > 0,
    100 * Reads_after_Trimmomatic / Raw_reads,
    NA_real_
  )
)

QC$Trimmomatic_removed_pct <- with(
  QC,
  ifelse(
    !is.na(Raw_reads) & Raw_reads > 0,
    100 * Reads_removed_Trimmomatic / Raw_reads,
    NA_real_
  )
)

QC$Human_reads_pct <- with(
  QC,
  ifelse(
    !is.na(Bowtie2_reads) & Bowtie2_reads > 0,
    100 * Human_reads / Bowtie2_reads,
    NA_real_
  )
)

QC$Nonhuman_reads_pct <- with(
  QC,
  ifelse(
    !is.na(Bowtie2_reads) & Bowtie2_reads > 0,
    100 * Nonhuman_reads / Bowtie2_reads,
    NA_real_
  )
)

QC$Nonhuman_retention_pct <- with(
  QC,
  ifelse(
    !is.na(Raw_reads) & Raw_reads > 0,
    100 * Nonhuman_reads / Raw_reads,
    NA_real_
  )
)


# ============================================================================
# 9. COMPROBACIONES DE CONSISTENCIA
# ============================================================================

QC$Trimmomatic_check <- with(
  QC,
  ifelse(
    !is.na(Raw_reads) &
    !is.na(Reads_after_Trimmomatic) &
    !is.na(Reads_removed_Trimmomatic),
    Raw_reads ==
      Reads_after_Trimmomatic +
      Reads_removed_Trimmomatic,
    NA
  )
)

QC$Bowtie2_check <- with(
  QC,
  ifelse(
    !is.na(Bowtie2_reads) &
    !is.na(Reads_after_Trimmomatic),
    Bowtie2_reads == Reads_after_Trimmomatic,
    NA
  )
)


# ============================================================================
# 10. ORDENAR COLUMNAS
# ============================================================================

QC <- QC[, c(
  "Sample",
  "Raw_reads",
  "Reads_after_Trimmomatic",
  "Reads_removed_Trimmomatic",
  "Trimmomatic_retention_pct",
  "Trimmomatic_removed_pct",
  "Bowtie2_reads",
  "Human_reads",
  "Human_reads_pct",
  "Nonhuman_reads",
  "Nonhuman_reads_pct",
  "Nonhuman_retention_pct",
  "Reads_removed_total",
  "Trimmomatic_check",
  "Bowtie2_check",
  "FastQC_raw",
  "FastQC_trimmed",
  "FastQC_nonhuman"
)]


# ============================================================================
# 11. GUARDAR TABLA
# ============================================================================

write.table(
  QC,
  file = OUTPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("\n============================================\n")
cat("RESUMEN DE QC COMPLETADO\n")
cat("============================================\n")

cat("Muestras:", nrow(QC), "\n")
cat("Archivo generado:", OUTPUT_FILE, "\n\n")

print(QC)
