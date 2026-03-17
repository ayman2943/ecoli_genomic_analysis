# Usage Guide

## Table of Contents
1. [Quick Start](#quick-start)
2. [Data Preparation](#data-preparation)
3. [Running Analyses](#running-analyses)
4. [Customization](#customization)
5. [Troubleshooting](#troubleshooting)

## Quick Start

### Minimal Working Example

```bash
# 1. Clone repository
git clone https://github.com/yourusername/ecoli-genomic-analysis.git
cd ecoli-genomic-analysis

# 2. Install dependencies
Rscript -e "source('requirements.txt')"

# 3. Prepare your data (see Data Preparation below)

# 4. Run analyses
Rscript scripts/01_card_analysis.R
Rscript scripts/02_vfdb_analysis.R
Rscript scripts/03_plasmid_analysis.R
```

## Data Preparation

### Directory Structure

Your data must follow this structure:

```
ecoli-genomic-analysis/
├── data/
│   ├── card_summary/
│   │   ├── ST10/
│   │   │   ├── Escherichia_coli_genome1_card.tsv
│   │   │   ├── Escherichia_coli_genome2_card.tsv
│   │   │   └── ...
│   │   ├── ST131/
│   │   ├── ST69/
│   │   ├── ST73/
│   │   └── ST95/
│   ├── vfdb_summary/
│   │   ├── ST10/
│   │   │   ├── Escherichia_coli_genome1_vfdb.tsv
│   │   │   └── ...
│   │   ├── ST131/
│   │   ├── ST69/
│   │   ├── ST73/
│   │   └── ST95/
│   ├── metadata/
│   │   ├── ST10_filtered.xlsx
│   │   ├── ST131_filtered.xlsx
│   │   ├── ST69_filtered.xlsx
│   │   ├── ST73_filtered.xlsx
│   │   └── ST95_filtered.xlsx
│   ├── plasmid_CARD_merged.tsv
│   └── plasmid_VFDB_merged.tsv
```

### File Formats

#### 1. CARD Summary Files (.tsv)

Each file should contain:
- Gene presence/absence or counts
- One row per genome
- Columns: Gene names (e.g., `aac(3)-IIa`, `blaCTX-M-15`)

Example:
```
aac(3)-IIa	aac(6')-Ib-cr	blaCTX-M-15	sul1	tet(A)
1	0	1	1	0
0	1	1	0	1
```

#### 2. VFDB Summary Files (.tsv)

Similar format to CARD, but with virulence genes:
- Columns: VF gene names (e.g., `fimH`, `iucD`, `ompA`)

#### 3. Metadata Files (.xlsx)

Required columns:
- `Name` or `Genome_ID` - Must match genome IDs in summary files
- `Collection Year` or `Collection_Year` - Numeric year
- `Country` - Country of isolation
- `Continent` - Geographic continent
- `Source Niche` or `Source_Niche` - Source category

Optional columns:
- `Source_Category` - Broader categorization
- `Host` - Host organism
- `Clinical_Status` - "Clinical" or "Non-clinical"

Example:
```
Name                    Collection_Year  Country        Continent      Source_Niche
Sample_001              2018            USA            North America  Human_Clinical
Sample_002              2019            UK             Europe         Animal
Sample_003              2020            India          Asia           Environmental
```

#### 4. Plasmid Files (.tsv)

Merged plasmid data for all STs:
- `plasmid_CARD_merged.tsv` - Plasmid ARG data
- `plasmid_VFDB_merged.tsv` - Plasmid VF data

Required columns:
- `Name` - Genome identifier
- `ST` - Sequence type (ST10, ST131, etc.)
- `Collection Year` - Numeric year
- `Country`, `Continent` - Geographic info
- Additional columns: Gene presence/absence

## Running Analyses

### 1. AMR Analysis

```bash
Rscript scripts/01_card_analysis.R
```

**Outputs:**
- `outputs/R_analysis_outputs/plots/` - 50+ publication-quality figures
- `outputs/R_analysis_outputs/tables/AMR_COMPLETE_SUMMARY.xlsx` - 21 worksheets
- `outputs/R_analysis_outputs/tables/Publication_Tables_CARD.docx` - Formatted tables

**Key Analyses:**
- Temporal trends (Mann-Kendall)
- ST-specific enrichment (log2 fold-change)
- Drug class analysis
- Clinical vs. non-clinical comparison
- Diversity metrics (Shannon, NMDS)
- Co-occurrence networks

**Runtime:** ~15-30 minutes (depends on dataset size)

### 2. Virulence Factor Analysis

```bash
Rscript scripts/02_vfdb_analysis.R
```

**Outputs:**
- `outputs/VFDB_analysis_outputs/plots/` - VF visualizations
- `outputs/VFDB_analysis_outputs/tables/VFDB_COMPLETE_SUMMARY.xlsx`

**Key Analyses:**
- VF class composition
- Gene prevalence profiles
- Temporal dynamics
- Clinical enrichment per ST
- Shannon diversity

**Runtime:** ~10-20 minutes

### 3. Plasmid Analysis

```bash
Rscript scripts/03_plasmid_analysis.R
```

**Outputs:**
- `outputs/plasmid_outputs/` - 5 key visualizations
- `plasmid_summary.xlsx` - Statistical summaries

**Key Analyses:**
- Plasmid ARG/VF burden
- MDR plasmid prevalence
- High-risk gene detection (ESBL, carbapenemases, MCR)
- ARG-VF correlation

**Runtime:** ~5-10 minutes

## Customization

### Modifying Sequence Types

Edit the `STs` variable in each main script:

```r
# Default
STs <- c("ST10", "ST131", "ST69", "ST73", "ST95")

# Custom
STs <- c("ST101", "ST127", "ST405")
```

### Changing Output Directories

Modify directory paths in main scripts:

```r
OUT_DIR <- "outputs/custom_output_directory"
```

### Adjusting Statistical Thresholds

Edit in respective module files:

```r
# Example: significance level
alpha <- 0.01  # Default is 0.05

# Example: minimum prevalence for gene inclusion
min_prevalence <- 0.1  # 10% instead of default 5%
```

### Custom Color Palettes

Modify in `scripts/utils/common_utils.R`:

```r
ST_COLORS <- c(
  "ST10"  = "#YOUR_COLOR_1",
  "ST131" = "#YOUR_COLOR_2",
  # ... etc
)
```

### Plot Dimensions

Adjust in individual module files or globally:

```r
# In save_plot() calls
save_plot(plot, "name", width = 14, height = 10)  # Default: 12 x 7
```

## Troubleshooting

### Common Issues

#### 1. Metadata Matching Failures

**Symptom:** Low match rates (<90%)

**Solutions:**
- Check genome ID format in summary files vs. metadata
- Verify column names in metadata (`Name` vs `Genome_ID`)
- Review `batch_match()` strategies in `scripts/utils/metadata_matching.R`
- Add custom matching patterns if needed

#### 2. Missing Packages

**Symptom:** `Error: package 'X' is not installed`

**Solution:**
```r
install.packages("package_name")
# or for all packages:
source("requirements.txt")
```

#### 3. Memory Issues

**Symptom:** `Error: cannot allocate vector of size X`

**Solutions:**
- Reduce parallel cores: `n_cores <- 2`
- Process subsets of data
- Increase system RAM allocation
- Use `data.table` instead of `data.frame` for large datasets

#### 4. Plot Rendering Issues

**Symptom:** Plots not displaying or garbled text

**Solutions:**
- Ensure graphics device is available
- Update graphics packages: `update.packages(c("ggplot2", "grid"))`
- Check font availability on system
- Use `cairo_pdf` device for better font support

#### 5. File Path Issues (Windows)

**Symptom:** `Error: cannot open the connection`

**Solutions:**
- Use forward slashes: `data/card_summary/` not `data\card_summary\`
- Or use `file.path()`: `file.path("data", "card_summary")`
- Set working directory: `setwd("path/to/project")`

### Getting Help

1. Check existing [GitHub Issues](https://github.com/yourusername/ecoli-genomic-analysis/issues)
2. Review error messages carefully
3. Run with verbose output: `options(verbose = TRUE)`
4. Open a new issue with:
   - Error message
   - Session info: `sessionInfo()`
   - Minimal reproducible example

## Advanced Usage

### Running Specific Modules Only

```r
# Load utilities
source("scripts/utils/common_utils.R")
source("scripts/utils/metadata_matching.R")

# Run only temporal analysis
source("scripts/modules/amr/05_temporal_analysis.R")
temporal_results <- run_temporal_analysis(data, STs, OUT_DIR)
```

### Batch Processing Multiple Datasets

```bash
#!/bin/bash
# batch_analysis.sh

for dataset in dataset1 dataset2 dataset3; do
  echo "Processing $dataset..."
  cd $dataset
  Rscript ../scripts/01_card_analysis.R
  Rscript ../scripts/02_vfdb_analysis.R
  Rscript ../scripts/03_plasmid_analysis.R
  cd ..
done
```

### Parallel Execution

```r
# Use all available cores
n_cores <- parallel::detectCores()

# Or specify manually
n_cores <- 8
```

## Performance Optimization

### For Large Datasets (>10,000 genomes)

1. **Use data.table for loading:**
```r
library(data.table)
data <- fread("large_file.tsv")
```

2. **Subset by ST first:**
```r
data_ST10 <- data %>% filter(ST == "ST10")
# Process ST10 separately
```

3. **Reduce visualization complexity:**
```r
# Sample for visualization
data_sample <- data %>% sample_n(5000)
```

4. **Use parallel processing:**
```r
library(parallel)
results <- mclapply(gene_list, analyze_gene, mc.cores = n_cores)
```

## Output Interpretation

See [docs/OUTPUT_GUIDE.md](OUTPUT_GUIDE.md) for detailed explanation of:
- Statistical test interpretations
- Plot descriptions
- Table contents
- How to use results in publications
