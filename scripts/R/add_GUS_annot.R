# Usage:
# Rscript add_GUS_annotations.R input_GUS.tsv GUS_classes.tsv

args <- commandArgs(trailingOnly = TRUE)

if (length(args)!= 2) {
  stop("Usage: Rscript add_GUS_annotations.R input_GUS.tsv GUS_classes.tsv")
}

diamond_file <- args[1]
classes_file <- args[2]

# Read DIAMOND output (12 original columns, no header)
df <- read.delim(
  diamond_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(df) <- c(
  "qseqid", "sseqid", "pident", "length",
  "qstart", "qend", "sstart", "send",
  "evalue", "bitscore", "qlen", "slen"
)

# Read GUS structural classes
classes <- read.delim(
  classes_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Query coverage (%)
df$query_coverage <- (
  abs(df$qend - df$qstart) + 1
) / df$qlen * 100

# Subject coverage (%)
df$subject_coverage <- (
  abs(df$send - df$sstart) + 1
) / df$slen * 100

# Add GUS structural class
class_dict <- setNames(classes$gus_type, classes$new_id)

df$gus_class <- unname(class_dict[df$sseqid])
df$gus_class[is.na(df$gus_class)] <- "Unknown"

# Count distinct GUS proteins per read
gus_counts <- aggregate(
  sseqid ~ qseqid,
  data = df,
  FUN = function(x) length(unique(x))
)

names(gus_counts)[2] <- "n_distinct_GUS"

df$n_distinct_GUS <- gus_counts$n_distinct_GUS[
  match(df$qseqid, gus_counts$qseqid)
]

# Unique or multiple GUS matches
df$match_type <- ifelse(
  df$n_distinct_GUS == 1,
  "Unique_GUS",
  "Multi_GUS"
)

# Round coverage
df$query_coverage <- round(df$query_coverage, 2)
df$subject_coverage <- round(df$subject_coverage, 2)

# Backup original file
backup <- paste0(diamond_file, ".backup")

if (!file.exists(backup)) {
  file.copy(diamond_file, backup)
}

# Overwrite original file with annotations
write.table(
  df,
  file = diamond_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Summary
cat("\nGUS annotation completed\n")
cat("------------------------\n")
cat("Total alignments:", nrow(df), "\n")
cat("Reads with GUS hits:", length(unique(df$qseqid)), "\n")
cat("Unique GUS reads:",
    sum(gus_counts$n_distinct_GUS == 1), "\n")
cat("Multiple GUS reads:",
    sum(gus_counts$n_distinct_GUS > 1), "\n")
cat("Updated file:", diamond_file, "\n")
