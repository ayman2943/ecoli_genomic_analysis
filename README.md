# E. coli Genomic Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue)](https://www.r-project.org/)

A comprehensive R-based pipeline for analyzing antimicrobial resistance (AMR), virulence factors (VF), and plasmid-associated genes in *Escherichia coli* genome sequences across multiple sequence types (ST10, ST131, ST69, ST73, ST95).

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Components](#pipeline-components)
- [Output](#output)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [License](#license)
- [Contributing](#contributing)

## 🔬 Overview

This pipeline integrates three major analytical components:

1. **AMR Analysis (CARD)** - Comprehensive Antibiotic Resistance Database analysis
2. **Virulence Factor Analysis (VFDB)** - Virulence Factor Database profiling
3. **Plasmid Analysis** - Plasmid-borne resistance and virulence characterization

The pipeline performs publication-quality statistical analyses including temporal trends, ST-specific enrichment, clinical vs. non-clinical comparisons, and multi-dimensional visualizations.

## ✨ Features

### AMR Analysis
- Multi-strategy genome ID reconciliation with metadata
- Temporal trend analysis (Mann-Kendall tau, linear models)
- ST-specific gene enrichment (log2 fold-change, Fisher's exact test)
- Drug class-specific analyses (aminoglycosides, β-lactams)
- Clinical vs. non-clinical source comparisons
- MDR (Multi-Drug Resistance) profiling
- Co-occurrence network analysis
- Diversity metrics (Shannon, NMDS, PERMANOVA)

### Virulence Factor Analysis
- VF class composition and temporal dynamics
- Gene prevalence profiling across STs
- Shannon diversity of VF profiles
- Clinical enrichment analysis per ST
- Co-occurrence correlation matrices
- Statistical frameworks (Kruskal-Wallis, Fisher's exact, logistic regression)

### Plasmid Analysis
- Plasmid-borne ARG burden tracking
- MDR plasmid prevalence
- High-risk gene detection (ESBL, carbapenemases, MCR)
- VF-ARG correlation on plasmids
- Temporal trend visualization

## 🚀 Installation

### Prerequisites

- R (≥ 4.0.0)
- RStudio (recommended)

### Clone the Repository

```bash
git clone https://github.com/yourusername/ecoli-genomic-analysis.git
cd ecoli-genomic-analysis
```

### Install Dependencies

```r
# Install CRAN packages
install.packages(c(
  "tidyverse", "readxl", "writexl", "scales", "RColorBrewer",
  "viridis", "ggpubr", "gridExtra", "patchwork", "ggsci",
  "data.table", "parallel", "stringi", "ggridges", "ggtext",
  "broom", "ggrepel", "vegan", "MASS", "rstatix", "igraph",
  "ggraph", "tidygraph", "corrplot", "flextable", "officer",
  "epitools", "tidytext", "colorspace", "pheatmap"
))
```

## 📊 Usage

### Quick Start

1. **Prepare your data structure:**

```
project_root/
├── card_summary/          # CARD analysis TSV files
│   ├── ST10/
│   ├── ST131/
│   ├── ST69/
│   ├── ST73/
│   └── ST95/
├── vfdb_summary/          # VFDB analysis TSV files
│   ├── ST10/
│   ├── ST131/
│   ├── ST69/
│   ├── ST73/
│   └── ST95/
├── metadata/              # Metadata Excel files
│   ├── ST10_filtered.xlsx
│   ├── ST131_filtered.xlsx
│   ├── ST69_filtered.xlsx
│   ├── ST73_filtered.xlsx
│   └── ST95_filtered.xlsx
├── plasmid_CARD_merged.tsv
└── plasmid_VFDB_merged.tsv
```

2. **Run the analyses:**

```bash
# AMR Analysis
Rscript scripts/01_card_analysis.R

# Virulence Factor Analysis
Rscript scripts/02_vfdb_analysis.R

# Plasmid Analysis
Rscript scripts/03_plasmid_analysis.R
```

### Advanced Usage

See [docs/USAGE.md](docs/USAGE.md) for detailed parameter configuration and customization options.

## 🔧 Pipeline Components

### 1. AMR Analysis (`scripts/01_card_analysis.R`)

**Modules:**
- `01_metadata_matching.R` - Multi-strategy ID reconciliation
- `02_data_loading.R` - CARD data import and validation
- `03_base_visualizations.R` - Core AMR plots
- `04_enrichment_analysis.R` - ST-specific gene enrichment
- `05_temporal_analysis.R` - Time-series statistical models
- `06_clinical_comparison.R` - Clinical vs. non-clinical analysis
- `07_diversity_analysis.R` - Shannon, NMDS, PERMANOVA
- `08_cooccurrence_network.R` - Gene co-occurrence patterns
- `09_statistical_models.R` - GLM, logistic regression
- `10_export_tables.R` - Summary table generation

### 2. Virulence Factor Analysis (`scripts/02_vfdb_analysis.R`)

**Modules:**
- `01_vfdb_loading.R` - VFDB data import
- `02_vf_class_analysis.R` - Class composition profiling
- `03_prevalence_analysis.R` - Gene prevalence calculations
- `04_temporal_vf_trends.R` - VF temporal dynamics
- `05_clinical_vf_comparison.R` - Clinical enrichment per ST
- `06_diversity_metrics.R` - VF profile diversity
- `07_statistical_tests.R` - Comprehensive statistical framework
- `08_visualizations.R` - Publication-quality plots

### 3. Plasmid Analysis (`scripts/03_plasmid_analysis.R`)

**Modules:**
- `01_plasmid_loading.R` - Plasmid data processing
- `02_arg_burden.R` - ARG burden quantification
- `03_highrisk_detection.R` - ESBL/carbapenemase/MCR detection
- `04_vf_plasmid.R` - Plasmid-borne virulence factors
- `05_correlation_analysis.R` - ARG-VF co-occurrence

## 📈 Output

### Directory Structure

```
outputs/
├── R_analysis_outputs/
│   ├── plots/                    # PNG visualizations (300 DPI)
│   │   ├── M01-M05_*.png        # Metadata plots
│   │   ├── 01-10_*.png          # Base AMR visualizations
│   │   ├── A1-A3_*.png          # ST enrichment
│   │   ├── B1-B3_*.png          # Aminoglycoside analysis
│   │   ├── C1-C3_*.png          # β-lactam analysis
│   │   ├── D1-D6b_*.png         # Clinical comparisons
│   │   ├── E1-E2_*.png          # Diversity metrics
│   │   ├── F1-F2_*.png          # Co-occurrence networks
│   │   └── G1-G2_*.png          # Statistical models
│   └── tables/                   # Statistical outputs
│       ├── AMR_COMPLETE_SUMMARY.xlsx  # 21 worksheets
│       ├── Publication_Tables_CARD.docx
│       └── *.csv                # Individual CSV tables
├── VFDB_analysis_outputs/
│   ├── plots/
│   │   ├── V01-V15_*.png        # VF visualizations
│   │   └── combined_*.png       # Multi-panel figures
│   ├── tables/
│   │   └── VFDB_COMPLETE_SUMMARY.xlsx
│   └── stats/                    # Statistical test results
└── plasmid_outputs/
    ├── P1-P5_*.png              # Plasmid visualizations
    └── plasmid_summary.xlsx
```

### Key Outputs

- **Publication Tables**: Formatted Word/Excel documents with Table 1 (study population) and Table 2 (gene prevalence)
- **Statistical Summaries**: Comprehensive Excel workbooks with 20+ analysis sheets
- **High-Resolution Plots**: 300 DPI PNG images for publication
- **CSV Exports**: Individual tables for external analysis

## 📦 Dependencies

### Core Packages
- `tidyverse` (≥ 2.0.0) - Data manipulation and visualization
- `data.table` (≥ 1.14.0) - High-performance data processing
- `ggplot2` (≥ 3.4.0) - Advanced plotting

### Statistical Analysis
- `vegan` - Ecological diversity metrics
- `MASS` - Statistical modeling
- `rstatix` - Statistical tests
- `broom` - Model tidying
- `epitools` - Epidemiological tools

### Visualization
- `patchwork` - Multi-panel plot composition
- `ggpubr` - Publication-ready themes
- `viridis` / `RColorBrewer` - Color palettes
- `ggraph` / `igraph` - Network visualization
- `corrplot` - Correlation matrices
- `pheatmap` - Heatmap generation

### Reporting
- `flextable` / `officer` - Word document generation
- `writexl` / `readxl` - Excel I/O

See [requirements.txt](requirements.txt) for complete version specifications.

## 📚 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{ecoli_genomic_analysis_2025,
  author = {Your Name},
  title = {E. coli Genomic Analysis Pipeline},
  year = {2025},
  url = {https://github.com/yourusername/ecoli-genomic-analysis},
  version = {1.0.0}
}
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🤝 Contributing

Contributions are welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

### Development Guidelines

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 🐛 Issues

Found a bug or have a feature request? Please open an issue on [GitHub Issues](https://github.com/yourusername/ecoli-genomic-analysis/issues).

## 👥 Authors

- **Your Name** - *Initial work* - [YourGitHub](https://github.com/yourusername)

## 🙏 Acknowledgments

- CARD Database - Comprehensive Antibiotic Resistance Database
- VFDB - Virulence Factor Database
- The R community for excellent statistical packages

## 📞 Contact

For questions or collaboration inquiries, please contact: your.email@example.com

---

**Last Updated:** March 2025
