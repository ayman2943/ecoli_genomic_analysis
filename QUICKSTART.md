# Quick Start Guide

Get up and running with the E. coli Genomic Analysis Pipeline in 5 minutes.

## Prerequisites

- R version ≥ 4.0.0
- 8 GB RAM minimum
- 10 GB free disk space

## Installation

### Step 1: Clone Repository

```bash
git clone https://github.com/yourusername/ecoli-genomic-analysis.git
cd ecoli-genomic-analysis
```

### Step 2: Install R Packages

```bash
Rscript install.R
```

This will install all 30+ required packages automatically.

## Prepare Your Data

### Step 3: Organize Data Files

Place your files in the `data/` directory:

```
data/
├── card_summary/
│   ├── ST10/
│   ├── ST131/
│   ├── ST69/
│   ├── ST73/
│   └── ST95/
├── vfdb_summary/
│   ├── ST10/
│   ├── ST131/
│   ├── ST69/
│   ├── ST73/
│   └── ST95/
├── metadata/
│   ├── ST10_filtered.xlsx
│   ├── ST131_filtered.xlsx
│   ├── ST69_filtered.xlsx
│   ├── ST73_filtered.xlsx
│   └── ST95_filtered.xlsx
├── plasmid_CARD_merged.tsv
└── plasmid_VFDB_merged.tsv
```

**Don't have data yet?** Use the example data:

```bash
cp -r examples/data/* data/
```

## Run Analysis

### Step 4: Execute Pipeline

Run all three analyses:

```bash
# AMR Analysis (15-30 min)
Rscript scripts/01_card_analysis.R

# Virulence Factor Analysis (10-20 min)
Rscript scripts/02_vfdb_analysis.R

# Plasmid Analysis (5-10 min)
Rscript scripts/03_plasmid_analysis.R
```

Or run them all at once:

```bash
for script in scripts/0{1,2,3}_*.R; do
  Rscript "$script"
done
```

## View Results

### Step 5: Check Outputs

Your results are in the `outputs/` directory:

```
outputs/
├── R_analysis_outputs/
│   ├── plots/              # 50+ PNG visualizations
│   └── tables/             # Excel & CSV summaries
├── VFDB_analysis_outputs/
│   ├── plots/              # VF visualizations
│   └── tables/             # VF summaries
└── plasmid_outputs/
    ├── P1-P5_*.png        # Plasmid plots
    └── plasmid_summary.xlsx
```

## Key Output Files

### Publication-Ready Materials

1. **Main Summary**
   - `outputs/R_analysis_outputs/tables/AMR_COMPLETE_SUMMARY.xlsx`
   - 21 worksheets with all statistical results

2. **Publication Tables**
   - `outputs/R_analysis_outputs/tables/Publication_Tables_CARD.docx`
   - Formatted Table 1 (study population) and Table 2 (gene prevalence)

3. **Combined Figures**
   - `outputs/R_analysis_outputs/plots/COMBINED_publication_figure.png`
   - Multi-panel figure ready for manuscripts

### Key Visualizations

- **M01-M05**: Metadata and sample distributions
- **01-10**: Core AMR burden and prevalence plots
- **A1-A3**: ST-specific gene enrichment
- **D1-D6b**: Clinical vs. non-clinical comparisons
- **E1-E2**: Diversity analysis (Shannon, NMDS)
- **F1-F2**: Gene co-occurrence networks
- **V01-V15**: Virulence factor analyses
- **P1-P5**: Plasmid-associated genes

## Troubleshooting

### Common Issues

**Problem:** Package installation fails

**Solution:**
```r
# Install packages manually
install.packages(c("tidyverse", "ggplot2", "vegan"))
```

**Problem:** Low metadata matching rate

**Solution:** Check that genome IDs in your CARD/VFDB files match the `Name` column in your metadata Excel files.

**Problem:** Out of memory

**Solution:** Reduce parallel cores in the scripts:
```r
n_cores <- 2  # Instead of max cores
```

## Next Steps

### Customize Your Analysis

1. **Modify sequence types** (if not using ST10/131/69/73/95):
   ```r
   # In each main script
   STs <- c("ST101", "ST127", "ST405")
   ```

2. **Adjust plot aesthetics**:
   - Edit `scripts/utils/common_utils.R`
   - Modify color palettes, themes, dimensions

3. **Run specific modules only**:
   ```r
   source("scripts/utils/common_utils.R")
   source("scripts/modules/amr/05_temporal_analysis.R")
   results <- run_temporal_analysis(data, STs, OUT_DIR)
   ```

### Learn More

- [Detailed Usage Guide](docs/USAGE.md) - Complete documentation
- [Project Structure](docs/STRUCTURE.md) - Directory organization
- [Contributing Guidelines](CONTRIBUTING.md) - How to contribute

## Getting Help

1. Check [documentation](docs/)
2. Search [existing issues](https://github.com/yourusername/ecoli-genomic-analysis/issues)
3. Open a [new issue](https://github.com/yourusername/ecoli-genomic-analysis/issues/new)

## Example Workflow

Complete example from start to finish:

```bash
# 1. Setup
git clone https://github.com/yourusername/ecoli-genomic-analysis.git
cd ecoli-genomic-analysis
Rscript install.R

# 2. Add your data
# ... copy your files to data/ directory ...

# 3. Run analysis
Rscript scripts/01_card_analysis.R
Rscript scripts/02_vfdb_analysis.R
Rscript scripts/03_plasmid_analysis.R

# 4. Review results
ls outputs/R_analysis_outputs/plots/
open outputs/R_analysis_outputs/tables/AMR_COMPLETE_SUMMARY.xlsx

# 5. Generate manuscript figures
# Figures are already publication-ready at 300 DPI!
```

## What You Get

After running all analyses, you'll have:

- ✅ 70+ publication-quality figures (300 DPI PNG)
- ✅ Comprehensive statistical summaries (Excel)
- ✅ Formatted publication tables (Word)
- ✅ Individual CSV exports for each analysis
- ✅ Complete reproducible workflow

## Time Investment

- **Setup:** 10 minutes
- **Data preparation:** 30-60 minutes (first time)
- **Analysis runtime:** 30-60 minutes
- **Results review:** 30+ minutes

**Total:** ~2-3 hours for complete analysis pipeline

## Success Criteria

You'll know it worked when:

1. ✅ All scripts run without errors
2. ✅ Output directories contain files
3. ✅ Plots render correctly
4. ✅ Excel files open without issues
5. ✅ Metadata matching rate >90%

---

**Questions?** See [USAGE.md](docs/USAGE.md) or open an issue on GitHub!
