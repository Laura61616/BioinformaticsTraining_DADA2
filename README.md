# Tele02 12S eDNA Metabarcoding Training

This beginner-friendly tutorial processes paired-end Tele02 12S eDNA reads with Cutadapt and DADA2. You will remove primers, filter reads, infer amplicon sequence variants (ASVs), remove chimeras, retain 100–200 bp ASVs, and assign taxonomy from Kingdom to Species.

The project contains two deterministic training samples with 200,000 read pairs per sample. Run the single R script section by section in RStudio.

## Project files

```text
edna-training-repo/
|-- README.md
|-- edna-training-repo.Rproj
|-- data/                         # Training FASTQ files and sample information
|-- reference/                    # Standardized DADA2 reference FASTA files for vertebrates based on the Tele02 primer set; mammalian sequences are excluded to reduce file size
|-- scripts/
|   `-- tele02_dada2_training.R  # The only R script
|-- results/                      # Result tables and intermediate FASTQ files
|-- figures/                      # PNG and PDF figures
`-- validation/
    `-- validation_report.md     # Expected results and figure-reading guide
```

## Installation

### 1. Prepare R

#### Windows

1. Install R from https://cran.r-project.org/bin/windows/base/
2. Install RStudio Desktop from https://posit.co/download/rstudio-desktop/

#### macOS

1. Install R from https://cran.r-project.org/bin/macosx/
2. Choose the installer matching your Mac processor.
3. Install RStudio Desktop from https://posit.co/download/rstudio-desktop/

#### Linux

Follow the instructions for your distribution at https://cran.r-project.org/bin/linux/

For Ubuntu, the basic command is:

```bash
sudo apt update
sudo apt install r-base r-base-dev
```

Then install RStudio Desktop from https://posit.co/download/rstudio-desktop/

#### Install the R packages

Open the RStudio Console and run:

```r
install.packages(c("BiocManager", "ggplot2", "ggsci", "scales"))
BiocManager::install(c("dada2", "ShortRead", "Biostrings"))
```

Check the installation:

```r
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
```

If these commands finish without an error, the R packages are ready.

### 2. Install Cutadapt

The official installation page is https://cutadapt.readthedocs.io/en/stable/installation.html

The examples below install Cutadapt in a separate Python environment. Replace `SET_YOUR_CUTADAPT_INSTALLATION_DIRECTORY` with a directory that you choose and can find again.

#### Windows

Install Python from https://www.python.org/downloads/ and run in PowerShell:

```powershell
$cutadapt_env = "SET_YOUR_CUTADAPT_INSTALLATION_DIRECTORY"
py -m venv $cutadapt_env
& "$cutadapt_env\Scripts\python.exe" -m pip install --upgrade pip
& "$cutadapt_env\Scripts\python.exe" -m pip install cutadapt
& "$cutadapt_env\Scripts\cutadapt.exe" --version
```

Record the complete path to this file:

```text
SET_YOUR_CUTADAPT_INSTALLATION_DIRECTORY/Scripts/cutadapt.exe
```

#### macOS and Linux

Make sure Python 3 is installed and run in Terminal:

```bash
CUTADAPT_ENV="SET_YOUR_CUTADAPT_INSTALLATION_DIRECTORY"
python3 -m venv "$CUTADAPT_ENV"
"$CUTADAPT_ENV/bin/python" -m pip install --upgrade pip
"$CUTADAPT_ENV/bin/python" -m pip install cutadapt
"$CUTADAPT_ENV/bin/cutadapt" --version
```

Record the complete path to this file:

```text
SET_YOUR_CUTADAPT_INSTALLATION_DIRECTORY/bin/cutadapt
```

Keep this path. The R script needs the path to the executable itself, not only the installation directory.

## Open the project in RStudio

1. Download and unzip the complete `edna-training-repo` folder.
2. Double-click `edna-training-repo.Rproj`.
3. Open `scripts/tele02_dada2_training.R`.
4. Run code with **Ctrl+Enter** on Windows/Linux or **Command+Enter** on macOS.

Complete Preparation and then Steps 1–4 in order.

## Preparation: set paths and load data

Change the first two path values in the R script.

`project_dir` must point to the downloaded `edna-training-repo` folder:

```r
project_dir <- "SET_YOUR_PROJECT_DIRECTORY/edna-training-repo"
```

`cutadapt` must point to the executable file. Use the line for your operating system:

```r
# Windows
cutadapt <- "SET_YOUR_CUTADAPT_INSTALLATION_DIRECTORY/Scripts/cutadapt.exe"

# macOS or Linux
cutadapt <- "SET_YOUR_CUTADAPT_INSTALLATION_DIRECTORY/bin/cutadapt"
```

Use forward slashes in R paths. Do not leave `SET_YOUR_...` in either value.

Run the complete Preparation section. It creates output folders, loads packages, reads `sample_manifest.tsv` and `primer.tsv`, and finds the FASTQ and taxonomy reference files. The Console should print the Cutadapt version and a two-row sample manifest.

## Step 1: remove primers

Run the section labelled:

```r
# Step 1: remove primers
```

This step checks Tele02 primer orientations, plots raw quality, removes primers with Cutadapt, checks the trimmed reads, and plots trimmed quality.

Check `primer_hits_before` and `primer_hits_after` in the Console. The expected orientations should be common before trimming and absent or nearly absent afterwards.

Open:

- `results/primer_removal_summary.tsv`
- `figures/step1_raw_quality_forward.png`
- `figures/step1_raw_quality_reverse.png`
- `figures/step1_trimmed_quality_forward.png`
- `figures/step1_trimmed_quality_reverse.png`

Trimmed FASTQ files are saved in `results/trimmed/`.

## Step 2: filter reads and learn errors

Run the section labelled:

```r
# Step 2: filter reads and learn errors
```

DADA2 removes reads with ambiguous bases or excessive expected errors, then learns separate forward and reverse error models.

Check `filter_stats` in the Console. An unexpectedly large loss suggests that the quality profiles and filtering parameters should be reviewed.

Open:

- `figures/step2_error_rates_forward.png`
- `figures/step2_error_rates_reverse.png`

Filtered FASTQ files are saved in `results/filtered/`.

## Step 3: infer ASVs, merge pairs and remove chimeras

Run the section labelled:

```r
# Step 3: infer ASVs, merge pairs and remove chimeras
```

This step infers sequence variants, merges paired reads with at least 20 bp overlap and no mismatches, removes chimeras, and retains only ASVs from 100 to 200 bp inclusive.

Check `read_tracking` in the Console and open:

- `results/read_tracking.tsv`
- `results/asv_table.tsv`
- `results/asvs.fa`
- `figures/step3_asv_length_distribution.png`

`asv_table.tsv` is the final sample-by-ASV count table. Each row contains the ASV ID, sequence, and read count in every sample.

## Step 4: assign taxonomy and export results

Run the section labelled:

```r
# Step 4: assign taxonomy and export results
```

`assignTaxonomy()` classifies each retained ASV from Kingdom to Species with `minBoot = 50`. The script uses `tryRC = TRUE`, so DADA2 also tests the reverse-complement orientation when assigning taxonomy.

This tutorial intentionally does not run `assignSpecies()`. The `Species` column is produced by `assignTaxonomy()` and should be interpreted together with `Species_bootstrap`; it is not a separate exact-match confirmation.

Open:

- `results/asv_taxonomy.tsv`
- `results/taxonomic_composition.tsv`
- `figures/step4_taxonomic_composition.png`

`asv_taxonomy.tsv` has one row for every retained ASV and contains:

- ASV ID, sequence, length, and total reads;
- Kingdom, Phylum, Class, Order, Family, Genus, and Species;
- a bootstrap value for every taxonomic rank.

An empty lower-rank field means that the sequence was not confidently resolved at that rank. Do not replace it with a guessed taxon.

## Main results

| File | Description |
|---|---|
| `primer_removal_summary.tsv` | Read pairs retained after Cutadapt |
| `read_tracking.tsv` | Read counts through all major steps |
| `asv_table.tsv` | ASV sequences and counts for each sample |
| `asvs.fa` | Retained ASV sequences in FASTA format |
| `asv_taxonomy.tsv` | One row per ASV with Kingdom–Species assignments and bootstraps |
| `taxonomic_composition.tsv` | Order-level read proportions used in the final plot |

### Expected output for the training data

| Stage | SAM001 | SAM002 | Total |
|---|---:|---:|---:|
| Raw read pairs | 200,000 | 200,000 | 400,000 |
| After primer removal | 169,384 | 142,393 | 311,777 |
| After quality filtering | 154,944 | 128,115 | 283,059 |
| Merged reads | 151,992 | 126,779 | 278,771 |
| Non-chimeric reads | 151,591 | 126,736 | 278,327 |
| Reads in 100–200 bp ASVs | 151,369 | 126,624 | 277,993 |

The validated run produced 81 non-chimeric ASVs before length filtering and retained 63 ASVs from 100 to 200 bp. With `assignTaxonomy(minBoot = 50, tryRC = TRUE)`, 19 of the 63 ASVs received a non-empty Species assignment. Small differences can occur with different DADA2, Cutadapt, R, or reference versions.

See `validation/validation_report.md` for a step-by-step guide to the expected tables, quality profiles, error-learning plots, ASV length plot, and taxonomy composition plot.

## Important interpretation notes

- An ASV is an inferred DNA sequence variant, not automatically a species.
- A taxonomic assignment depends on the reference database and bootstrap support.
- Read counts and proportions are molecular signals, not direct measurements of organism abundance, biomass, or density.
- These two training samples have no negative controls, so they cannot assess contamination or ecological absence.
- The reference source and redistribution license still require confirmation before external redistribution.

## References

- DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
- DADA2 taxonomy assignment: https://benjjneb.github.io/dada2/assign.html
- Cutadapt user guide: https://cutadapt.readthedocs.io/en/stable/guide.html
- Tele02 marker paper: https://doi.org/10.1111/2041-210X.13485
- Tutorial style reference: https://github.com/jingliu92/Amplicon_QIIME2
