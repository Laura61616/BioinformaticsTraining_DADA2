# Tele02 12S eDNA metabarcoding training pipeline
# Run this script from top to bottom in RStudio.

# ===============================
# Preparation: set paths and load data
# ===============================

# CHANGE ME: use the complete path to the edna-training-repo folder.
project_dir <- "SET_YOUR_PROJECT_DIRECTORY/edna-training-repo"

# CHANGE ME: use the complete path to cutadapt.exe (Windows) or cutadapt (macOS/Linux).
cutadapt <- "SET_YOUR_CUTADAPT_INSTALLATION_DIRECTORY/Scripts/cutadapt.exe"

data_dir <- file.path(project_dir, "data")
reference_dir <- file.path(project_dir, "reference")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
trimmed_dir <- file.path(results_dir, "trimmed")
filtered_dir <- file.path(results_dir, "filtered")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(trimmed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(filtered_dir, showWarnings = FALSE, recursive = TRUE)

library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)

manifest <- read.delim(file.path(data_dir, "sample_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
primers <- read.delim(file.path(data_dir, "primer.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

sample_names <- manifest$sample_id
fnFs <- file.path(data_dir, manifest$forward_file)
fnRs <- file.path(data_dir, manifest$reverse_file)
names(fnFs) <- sample_names
names(fnRs) <- sample_names

FWD <- primers$forward_primer_5to3[1]
REV <- primers$reverse_primer_5to3[1]
taxonomy_reference <- file.path(reference_dir, "tele02_dada2_taxonomy.fa.gz")

stopifnot(all(file.exists(fnFs)), all(file.exists(fnRs)))
stopifnot(file.exists(taxonomy_reference), file.exists(cutadapt))

system2(cutadapt, "--version")
manifest

# ===============================
# Step 1: remove primers
# ===============================

all_orientations <- function(primer) {
  dna <- DNAString(primer)
  c(
    Forward = as.character(dna),
    Complement = as.character(complement(dna)),
    Reverse = as.character(reverse(dna)),
    RevComp = as.character(reverseComplement(dna))
  )
}

primer_hits <- function(primer, fastq_file) {
  reads <- sread(readFastq(fastq_file))
  sum(vcountPattern(primer, reads, fixed = TRUE) > 0)
}

FWD_orientations <- all_orientations(FWD)
REV_orientations <- all_orientations(REV)

primer_hits_before <- rbind(
  FWD_in_R1 = sapply(FWD_orientations, primer_hits, fastq_file = fnFs[1]),
  FWD_in_R2 = sapply(FWD_orientations, primer_hits, fastq_file = fnRs[1]),
  REV_in_R1 = sapply(REV_orientations, primer_hits, fastq_file = fnFs[1]),
  REV_in_R2 = sapply(REV_orientations, primer_hits, fastq_file = fnRs[1])
)
primer_hits_before

raw_quality_F <- plotQualityProfile(fnFs)
raw_quality_R <- plotQualityProfile(fnRs)
ggsave(file.path(figures_dir, "step1_raw_quality_forward.png"), raw_quality_F, width = 8, height = 5, dpi = 300)
ggsave(file.path(figures_dir, "step1_raw_quality_forward.pdf"), raw_quality_F, width = 8, height = 5)
ggsave(file.path(figures_dir, "step1_raw_quality_reverse.png"), raw_quality_R, width = 8, height = 5, dpi = 300)
ggsave(file.path(figures_dir, "step1_raw_quality_reverse.pdf"), raw_quality_R, width = 8, height = 5)

FWD_RC <- as.character(reverseComplement(DNAString(FWD)))
REV_RC <- as.character(reverseComplement(DNAString(REV)))
adapter_R1 <- paste0(FWD, ";required...", REV_RC, ";optional")
adapter_R2 <- paste0(REV, ";required...", FWD_RC, ";optional")

trimFs <- file.path(trimmed_dir, basename(fnFs))
trimRs <- file.path(trimmed_dir, basename(fnRs))
names(trimFs) <- sample_names
names(trimRs) <- sample_names

for (i in seq_along(fnFs)) {
  cutadapt_args <- c(
    "-a", shQuote(adapter_R1),
    "-A", shQuote(adapter_R2),
    "-e", "0.1",
    "-O", "3",
    "--discard-untrimmed",
    "--pair-filter=any",
    "--minimum-length=50",
    "-o", shQuote(trimFs[i]),
    "-p", shQuote(trimRs[i]),
    shQuote(fnFs[i]),
    shQuote(fnRs[i])
  )
  if (system2(cutadapt, cutadapt_args) != 0) stop("Cutadapt failed for ", sample_names[i])
}

raw_reads <- sapply(fnFs, function(x) countFastq(x)$records)
trimmed_reads <- sapply(trimFs, function(x) countFastq(x)$records)

primer_removal <- data.frame(
  sample_id = sample_names,
  input_pairs = as.numeric(raw_reads),
  primer_trimmed_pairs = as.numeric(trimmed_reads),
  retained_percent = round(100 * trimmed_reads / raw_reads, 2)
)
write.table(primer_removal, file.path(results_dir, "primer_removal_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
primer_removal

primer_hits_after <- rbind(
  FWD_in_R1 = sapply(FWD_orientations, primer_hits, fastq_file = trimFs[1]),
  FWD_in_R2 = sapply(FWD_orientations, primer_hits, fastq_file = trimRs[1]),
  REV_in_R1 = sapply(REV_orientations, primer_hits, fastq_file = trimFs[1]),
  REV_in_R2 = sapply(REV_orientations, primer_hits, fastq_file = trimRs[1])
)
primer_hits_after

trimmed_quality_F <- plotQualityProfile(trimFs)
trimmed_quality_R <- plotQualityProfile(trimRs)
ggsave(file.path(figures_dir, "step1_trimmed_quality_forward.png"), trimmed_quality_F, width = 8, height = 5, dpi = 300)
ggsave(file.path(figures_dir, "step1_trimmed_quality_forward.pdf"), trimmed_quality_F, width = 8, height = 5)
ggsave(file.path(figures_dir, "step1_trimmed_quality_reverse.png"), trimmed_quality_R, width = 8, height = 5, dpi = 300)
ggsave(file.path(figures_dir, "step1_trimmed_quality_reverse.pdf"), trimmed_quality_R, width = 8, height = 5)

# ===============================
# Step 2: filter reads and learn errors
# ===============================

filtFs <- file.path(filtered_dir, basename(trimFs))
filtRs <- file.path(filtered_dir, basename(trimRs))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

filter_stats <- filterAndTrim(
  trimFs, filtFs,
  trimRs, filtRs,
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  truncLen = c(0, 0),
  minLen = 50,
  rm.phix = TRUE,
  compress = TRUE,
  matchIDs = TRUE,
  multithread = FALSE
)
filter_stats

set.seed(20261108)
errF <- learnErrors(filtFs, nbases = 50000000, randomize = TRUE, multithread = FALSE)
set.seed(20261108)
errR <- learnErrors(filtRs, nbases = 50000000, randomize = TRUE, multithread = FALSE)

error_plot_F <- plotErrors(errF, nominalQ = TRUE)
error_plot_R <- plotErrors(errR, nominalQ = TRUE)
ggsave(file.path(figures_dir, "step2_error_rates_forward.png"), error_plot_F, width = 8, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "step2_error_rates_forward.pdf"), error_plot_F, width = 8, height = 6)
ggsave(file.path(figures_dir, "step2_error_rates_reverse.png"), error_plot_R, width = 8, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "step2_error_rates_reverse.pdf"), error_plot_R, width = 8, height = 6)

# ===============================
# Step 3: infer ASVs, merge pairs and remove chimeras
# ===============================

dadaFs <- dada(filtFs, err = errF, pool = FALSE, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, pool = FALSE, multithread = FALSE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 20, maxMismatch = 0, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)

# Keep only Tele02 ASVs between 100 and 200 bp.
keep_asvs <- nchar(colnames(seqtab_nochim)) >= 100 & nchar(colnames(seqtab_nochim)) <= 200
seqtab_target <- seqtab_nochim[, keep_asvs, drop = FALSE]
table(nchar(getSequences(seqtab_target)))

getN <- function(x) sum(getUniques(x))

read_tracking <- data.frame(
  sample_id = sample_names,
  raw = as.numeric(raw_reads),
  primer_removed = as.numeric(trimmed_reads),
  filtered = as.numeric(filter_stats[, "reads.out"]),
  denoised_forward = sapply(dadaFs, getN),
  denoised_reverse = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  non_chimeric = rowSums(seqtab_nochim),
  target_length = rowSums(seqtab_target),
  row.names = NULL
)
write.table(read_tracking, file.path(results_dir, "read_tracking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
read_tracking

asv_order <- order(-colSums(seqtab_target), colnames(seqtab_target))
seqtab_final <- seqtab_target[, asv_order, drop = FALSE]
asv_sequences <- colnames(seqtab_final)
asv_ids <- paste0("ASV", seq_along(asv_sequences))

asv_table <- data.frame(
  ASV = asv_ids,
  Sequence = asv_sequences,
  as.data.frame(t(seqtab_final), check.names = FALSE),
  check.names = FALSE
)
write.table(asv_table, file.path(results_dir, "asv_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(as.vector(rbind(paste0(">", asv_ids), asv_sequences)), file.path(results_dir, "asvs.fa"))

asv_lengths <- data.frame(length = nchar(asv_sequences))
length_plot <- ggplot(asv_lengths, aes(x = length)) +
geom_histogram(binwidth = 1, boundary = 0, fill = "#3366A6", color = "white") +
geom_vline(xintercept = c(100, 200), color = "#B2182B", linetype = "dashed") +
labs(title = "Length distribution of retained ASVs", x = "ASV length (bp)", y = "Number of ASVs") +
theme_bw() +
theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "bold")
)
ggsave(file.path(figures_dir, "step3_asv_length_distribution.png"), length_plot, width = 8, height = 5, dpi = 300)
ggsave(file.path(figures_dir, "step3_asv_length_distribution.pdf"), length_plot, width = 8, height = 5)

# ===============================
# Step 4: assign taxonomy and export results
# ===============================

set.seed(20261108)
taxonomy_result <- assignTaxonomy(
  asv_sequences,
  taxonomy_reference,
  minBoot = 50,
  tryRC = TRUE,
  outputBootstraps = TRUE,
  multithread = FALSE,
  verbose = TRUE
)
taxonomy <- as.data.frame(taxonomy_result$tax, stringsAsFactors = FALSE)
bootstraps <- as.data.frame(taxonomy_result$boot, stringsAsFactors = FALSE)
names(bootstraps) <- paste0(names(bootstraps), "_bootstrap")

asv_taxonomy <- data.frame(
  ASV = asv_ids,
  Sequence = asv_sequences,
  Length = nchar(asv_sequences),
  Total_reads = colSums(seqtab_final),
  taxonomy,
  bootstraps,
  check.names = FALSE
)
write.table(asv_taxonomy, file.path(results_dir, "asv_taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Read proportions are molecular signals, not organism abundance.
order_labels <- ifelse(is.na(asv_taxonomy$Order) | asv_taxonomy$Order == "", "Unassigned", asv_taxonomy$Order)
composition_list <- vector("list", length(sample_names))

for (i in seq_along(sample_names)) {
  order_reads <- tapply(seqtab_final[i, ], order_labels, sum)
  composition_list[[i]] <- data.frame(sample_id = sample_names[i], Order = names(order_reads), Reads = as.numeric(order_reads))
}

taxonomic_composition <- do.call(rbind, composition_list)
taxonomic_composition$Read_proportion <- ave(taxonomic_composition$Reads, taxonomic_composition$sample_id, FUN = function(x) x / sum(x))
order_levels <- c(sort(setdiff(unique(taxonomic_composition$Order), "Unassigned")), "Unassigned")
taxonomic_composition$Order <- factor(taxonomic_composition$Order, levels = order_levels)
write.table(taxonomic_composition, file.path(results_dir, "taxonomic_composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

order_colors <- setNames(ggsci::pal_d3("category20")(length(order_levels)), order_levels)
composition_plot <- ggplot(taxonomic_composition, aes(x = sample_id, y = Read_proportion, fill = Order)) +
geom_col(width = 0.75) +
scale_fill_manual(values = order_colors, drop = FALSE) +
scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
labs(
  title = "Order-level composition of retained reads",
  x = "Sample",
  y = "Read proportion",
  caption = "Read proportions are molecular signals, not organism abundance."
) +
theme_bw() +
theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "bold"),
  legend.position = "right"
)
ggsave(file.path(figures_dir, "step4_taxonomic_composition.png"), composition_plot, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(figures_dir, "step4_taxonomic_composition.pdf"), composition_plot, width = 9, height = 5.5)

cat("\nPipeline complete.\n")
cat("Results:", results_dir, "\n")
cat("Figures:", figures_dir, "\n")
