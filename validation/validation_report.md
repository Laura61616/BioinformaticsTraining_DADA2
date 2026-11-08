# Validation and output-reading guide

**Project:** Tele02 12S eDNA Metabarcoding Training  
**Validation date:** 2026-11-08  
**Training data:** two samples, 200,000 paired reads per sample  
**Overall status:** passed

## What was validated

The Cutadapt-to-DADA2 pipeline was run on the supplied training FASTQ files. It completed primer removal, quality filtering, error learning, ASV inference, pair merging, chimera removal, the 100–200 bp ASV filter, taxonomy assignment, result export, and figure generation.

The final taxonomy stage was rerun on the validated 63-ASV table with minBoot = 50, tryRC = TRUE, and outputBootstraps = TRUE. assignSpecies() is not used in the final tutorial. Species values in results/asv_taxonomy.tsv come from assignTaxonomy() and have a corresponding Species_bootstrap.

The R script contains path placeholders. A learner must set project_dir and cutadapt before running it.

## Expected read tracking

| Stage | SAM001 | SAM002 | Total |
|---|---:|---:|---:|
| Raw read pairs | 200,000 | 200,000 | 400,000 |
| Primer removed | 169,384 | 142,393 | 311,777 |
| Quality filtered | 154,944 | 128,115 | 283,059 |
| Denoised forward | 154,932 | 128,108 | 283,040 |
| Denoised reverse | 154,930 | 128,106 | 283,036 |
| Merged | 151,992 | 126,779 | 278,771 |
| Non-chimeric | 151,591 | 126,736 | 278,327 |
| Retained in 100–200 bp ASVs | 151,369 | 126,624 | 277,993 |

Expected ASV totals:

- 81 non-chimeric ASVs before length filtering;
- 63 ASVs after retaining lengths from 100 through 200 bp;
- 63 rows in both asv_table.tsv and asv_taxonomy.tsv;
- 19 ASVs with a non-empty Species assignment at minBoot = 50 and tryRC = TRUE.

Minor numerical differences can occur if R, DADA2, Cutadapt, or the reference FASTA changes.

## Preparation: what to check

After running Preparation, the Console should show:

1. a Cutadapt version rather than an error;
2. a two-row manifest containing SAM001 and SAM002;
3. no failed stopifnot() check.

If stopifnot() fails, check that project_dir points to the complete edna-training-repo folder, cutadapt points to cutadapt.exe on Windows or cutadapt on macOS/Linux, the four training FASTQ files remain in data/, and reference/tele02_dada2_taxonomy.fa.gz exists.

## Step 1: remove primers

### Primer checks

primer_hits_before checks the first sample in four primer orientations. In the validated data, the strongest expected hits were:

| Check | Forward orientation | Reverse-complement orientation |
|---|---:|---:|
| Forward primer in R1 | 190,191 | 0 |
| Reverse primer in R2 | 172,408 | 0 |
| Forward primer in R2 | 0 | 24,567 |
| Reverse primer in R1 | 0 | 25,238 |

This pattern is consistent with paired-end reads covering the amplicon in opposite directions. After Cutadapt, primer_hits_after should be zero or very small. Persistent large values indicate that the adapter definitions, primer sequences, or read orientation should be checked.

### primer_removal_summary.tsv

| Sample | Input pairs | Primer-trimmed pairs | Retained |
|---|---:|---:|---:|
| SAM001 | 200,000 | 169,384 | 84.69% |
| SAM002 | 200,000 | 142,393 | 71.20% |

The samples do not need identical retention. A sample-specific loss can reflect incomplete primer detection, short inserts, sequencing quality, or biological/library differences. A near-zero result for every sample would indicate a likely primer or orientation problem.

### How to read the quality-profile figures

Files:

- step1_raw_quality_forward.png and .pdf;
- step1_raw_quality_reverse.png and .pdf;
- step1_trimmed_quality_forward.png and .pdf;
- step1_trimmed_quality_reverse.png and .pdf.

In a DADA2 quality profile:

- the x-axis is the base position in the read;
- the y-axis is the Phred quality score;
- gray shading shows how frequently each quality score occurs at each position;
- the green line is the mean quality score;
- orange lines show the quality-score quartiles;
- the red line is the scaled fraction of reads extending to that position.

The raw figures describe reads before primer removal. The trimmed figures begin after the required primer and can have variable terminal positions because Cutadapt removes adapters and discards unsuitable pairs.

Look for a broad or sudden quality collapse across many cycles. A gradual decline near read ends is common, and reverse reads are often poorer than forward reads. Enough high-quality sequence must remain for paired reads to overlap later.

The official DADA2 tutorial explains the same gray heat map and colored lines: https://benjjneb.github.io/dada2/tutorial.html

## Step 2: filter reads and learn errors

### Filtering

The validated output retained:

- 154,944 of 169,384 primer-trimmed SAM001 pairs;
- 128,115 of 142,393 primer-trimmed SAM002 pairs.

The tutorial uses maxN = 0, maxEE = c(2, 2), truncQ = 2, no fixed truncLen, and minLen = 50. Large losses should be interpreted together with the trimmed quality profiles before changing parameters.

### How to read the learnErrors figures

Files:

- step2_error_rates_forward.png and .pdf;
- step2_error_rates_reverse.png and .pdf.

Each panel represents a nucleotide transition such as A→C. For substitution panels:

- black points are observed error frequencies at each consensus quality score;
- the solid black line is the error rate fitted by DADA2;
- the red line is the theoretical error rate implied by the nominal Phred score.

A reasonable model has a black fitted line that broadly follows the observed points and generally decreases as quality increases. It does not need to overlap the red line: DADA2 is estimating run-specific errors instead of assuming perfect Phred calibration.

Sparse points at very high quality can make a fitted line look less smooth. The A→A, C→C, G→G, and T→T panels represent correct calls rather than substitution errors and should not be interpreted like off-diagonal panels.

learnErrors() normally reports the number of bases used and the self-consistency rounds. If convergence is reported, proceed. If it stops before convergence, inspect whether the black lines still follow the observations. The DADA2 FAQ notes that a nearly converged model is often usable when the visual fit is reasonable: https://benjjneb.github.io/dada2/faq.html

Warnings about zero, missing, or infinite values on a logarithmic scale can occur when a transition has no observation at a quality score. They are not automatically an analysis failure. A missing error plot, an invalid error object, or fitted lines that consistently miss dense observed points would be a real concern.

## Step 3: infer ASVs, merge pairs and remove chimeras

### read_tracking.tsv

Read this table from left to right. Counts usually decline because each column adds a quality or biological-consistency requirement.

- raw: input read pairs;
- primer_removed: pairs accepted by Cutadapt;
- filtered: pairs passing DADA2 quality filters;
- denoised_forward and denoised_reverse: reads assigned after error correction;
- merged: pairs with sufficient overlap and no mismatches;
- non_chimeric: reads remaining after chimera removal;
- target_length: reads in ASVs from 100 to 200 bp.

A large loss at one transition helps locate a problem. A large filtered-to-merged loss suggests insufficient overlap, orientation problems, or discordant pairs. A large merged-to-non_chimeric loss suggests many inferred chimeras or a problematic library. A large non_chimeric-to-target_length loss suggests off-target sequence lengths.

### asv_table.tsv

This is the final ASV count table. It has one row per retained ASV. ASV is the identifier, Sequence is the inferred nucleotide sequence, and each sample column contains integer read counts.

For each sample, the sum of ASV counts must equal target_length in read_tracking.tsv.

### ASV length figure

step3_asv_length_distribution.png and .pdf show the number of retained ASVs at each length. Dashed red lines mark 100 and 200 bp. Because filtering occurs before plotting, all bars must lie within these limits, including the boundaries.

The validated run retained 63 ASVs, and every exported sequence is 100–200 bp.

## Step 4: assign taxonomy and export results

### asv_taxonomy.tsv

This table contains one row for every row in asv_table.tsv. It includes:

1. ASV, Sequence, Length, and Total_reads;
2. Kingdom, Phylum, Class, Order, Family, Genus, and Species;
3. a bootstrap column for every taxonomic rank.

minBoot = 50 means a rank is retained only when it reaches the configured bootstrap threshold. Blank lower ranks mean unresolved classification, not confirmed biological absence.

tryRC = TRUE allows DADA2 to classify either the supplied ASV orientation or its reverse complement. The validated table has:

- 63 ASV rows;
- 25 non-empty Genus assignments;
- 19 non-empty Species assignments;
- no separate exact-species-match column.

Do not treat the Species field as proof that the marker uniquely identifies that organism. It is a classifier result supported by the supplied reference and bootstrap value.

### Taxonomic-composition figure

step4_taxonomic_composition.png and .pdf are stacked bar charts:

- each bar is one sample;
- colors are assigned Orders;
- segment height is the fraction of retained reads;
- each complete bar should sum to 100%;
- Unassigned means the ASV was not confidently classified at Order.

This figure summarizes read proportions, not organism abundance. Primer bias, mitochondrial copy number, DNA shedding, degradation, extraction, amplification, and reference coverage can all affect segment sizes.

## Final file checks

The validated project contains:

- results/primer_removal_summary.tsv;
- results/read_tracking.tsv;
- results/asv_table.tsv;
- results/asvs.fa;
- results/asv_taxonomy.tsv;
- results/taxonomic_composition.tsv;
- eight non-empty PNG figures and matching PDF figures.

Consistency checks passed:

- all 63 ASVs occur in both final tables;
- ASV sequences and IDs agree between the tables;
- every ASV length is from 100 to 200 bp;
- sample count sums agree with target_length;
- asv_taxonomy.tsv includes all seven ranks from Kingdom through Species;
- the obsolete Exact_species_match field is absent.

## Non-fatal warnings and limits

Windows can print locale warnings for C.UTF-8. DADA2 error plots can also report zero or missing values on a logarithmic scale. In the validated run these warnings did not stop the analysis or alter the expected counts.

This validation shows that the teaching workflow is internally consistent on the supplied data. It does not validate organism abundance, contamination, ecological absence, or reference-database completeness. The reference source and redistribution license still require confirmation before external redistribution.
