# methylation_array

Methylation_array is a Snakemake workflow for processing Illumina methylation array data. The pipeline performs quality control, creates a SummarizedExperiment object, generates PCA plots, and performs differential methylation analysis using the Bioconductor package [SeSAME](https://bioconductor.org/packages/release/bioc/html/sesame.html).

The workflow is optimized for the Van Andel Institute HPC environment but can be adapted to other Unix-based systems.

## Installation

Clone the repository:
```
git clone https://github.com/vari-bbc/methylation_array.git
```
Example data from [GEO GSE198222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198222) are included for testing.

## Configuration

To run the workflow on your own data:

### 1. Copy idat raw files

Copy raw idat files from the array into idat_files directory. 

---

### 2. Add the sample sheet

Place `SampleSheet.csv` in `samples/`.  
Count the number of rows before the `Sample_Name` header — this value will be used in the config file.

---

### 2. Edit `samples/meta.tsv`

This file links IDAT files to samples.

Required fields:

- **Group** – grouping variable for differential methylation analysis  
- **Sample_Name** – must match the first column of `SampleSheet.csv`

For multi-factor experiments (e.g. genotype × diet), combine factors into a single group label (e.g. `KO_hf`, `KO_sc`, `WT_hf`, `WT_sc`).

---

### 3. Edit `samples/comparisons.tsv`

Define comparisons for differential methylation.

- `group1` – test group  
- `group2` – reference group  

The analysis computes **group1 − group2**, so positive values indicate higher methylation in `group1`.

---

### 4. Edit `bin/config.yaml`

Key parameters:

- **samplesheet_header_rows** – number of header rows in the sample sheet
- **array_type**
  - `EPICv2`, `EPIC`, `HM450` (human)
  - `MM285` (mouse)
- **prep_method**
  - recommended: `QCDPB` (human arrays)
  - recommended: `TQCDPB` (MM285)

See the [SeSAME vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/sesame.html) for details.

## Running the workflow

From the repository root directory:

```bash
sbatch bin/run_snake.sh
```

Dry run:

```bash
module load bbc2/snakemake
snakemake -np
```

## Output

The workflow generates QC metrics, methylation summaries, differential methylation results, and visualization plots.

Key outputs include:

- **QC results** – per-chip quality metrics (`qc.rds`)
- **SummarizedExperiment (`se.rds`)** – beta values and sample metadata
- **Differential methylation results**
  - CpG-level statistics
  - Differentially methylated regions (DMRs)

Plots (in `analysis/plots`) include:

- **PCA plots** – sample clustering based on beta values
- **Volcano plots** – CpG-level methylation differences for each comparison

Note: p-values shown in volcano plots are **not adjusted for multiple testing** and are intended for visualization only.
