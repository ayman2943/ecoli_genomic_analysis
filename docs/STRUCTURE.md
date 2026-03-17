# Project Structure

Complete directory structure and file organization for the E. coli Genomic Analysis Pipeline.

## Overview

```
ecoli-genomic-analysis/
├── data/                          # Input data (gitignored)
├── docs/                          # Documentation
├── examples/                      # Example datasets
├── outputs/                       # Analysis outputs (gitignored)
├── scripts/                       # Analysis scripts
│   ├── modules/                   # Modular analysis components
│   │   ├── amr/                   # AMR-specific modules
│   │   ├── vfdb/                  # VFDB-specific modules
│   │   └── plasmid/               # Plasmid-specific modules
│   └── utils/                     # Shared utility functions
├── tests/                         # Unit tests (future)
├── .gitignore                     # Git ignore rules
├── CHANGELOG.md                   # Version history
├── CONTRIBUTING.md                # Contribution guidelines
├── install.R                      # Package installation script
├── LICENSE                        # MIT License
├── README.md                      # Main documentation
└── requirements.txt               # R package dependencies
```

## Detailed Structure

### `/data/` - Input Data Directory

**Purpose:** Stores all input data files (excluded from git)

```
data/
├── card_summary/                  # CARD analysis results
│   ├── ST10/
│   │   ├── Escherichia_coli_genome1_card.tsv
│   │   ├── Escherichia_coli_genome2_card.tsv
│   │   └── ...
│   ├── ST131/
│   ├── ST69/
│   ├── ST73/
│   └── ST95/
├── vfdb_summary/                  # VFDB analysis results
│   ├── ST10/
│   ├── ST131/
│   ├── ST69/
│   ├── ST73/
│   └── ST95/
├── metadata/                      # Genome metadata
│   ├── ST10_filtered.xlsx
│   ├── ST131_filtered.xlsx
│   ├── ST69_filtered.xlsx
│   ├── ST73_filtered.xlsx
│   └── ST95_filtered.xlsx
├── plasmid_CARD_merged.tsv       # Merged plasmid ARG data
└── plasmid_VFDB_merged.tsv       # Merged plasmid VF data
```

### `/docs/` - Documentation

**Purpose:** Additional documentation and guides

```
docs/
├── USAGE.md                       # Detailed usage guide
├── OUTPUT_GUIDE.md               # Output interpretation (future)
├── API.md                        # Function documentation (future)
└── STRUCTURE.md                  # This file
```

### `/examples/` - Example Data

**Purpose:** Sample datasets for testing and demonstration

```
examples/
├── data/                         # Example input files
│   ├── card_summary/
│   ├── vfdb_summary/
│   ├── metadata/
│   ├── plasmid_CARD_merged.tsv
│   └── plasmid_VFDB_merged.tsv
└── README.md                     # Example data guide
```

### `/outputs/` - Analysis Outputs

**Purpose:** Generated results (excluded from git)

```
outputs/
├── R_analysis_outputs/           # AMR analysis results
│   ├── plots/                    # PNG visualizations
│   │   ├── M01-M05_*.png        # Metadata plots
│   │   ├── 01-10_*.png          # Base AMR plots
│   │   ├── A1-A3_*.png          # Enrichment analysis
│   │   ├── B1-B3_*.png          # Aminoglycoside analysis
│   │   ├── C1-C3_*.png          # β-lactam analysis
│   │   ├── D1-D6b_*.png         # Clinical comparisons
│   │   ├── E1-E2_*.png          # Diversity analysis
│   │   ├── F1-F2_*.png          # Network analysis
│   │   └── G1-G2_*.png          # Statistical models
│   └── tables/                   # Data tables
│       ├── AMR_COMPLETE_SUMMARY.xlsx
│       ├── Publication_Tables_CARD.docx
│       └── *.csv                # Individual tables
├── VFDB_analysis_outputs/        # VF analysis results
│   ├── plots/
│   │   ├── V01-V15_*.png        # VF visualizations
│   │   └── combined_*.png       # Multi-panel figures
│   ├── tables/
│   │   └── VFDB_COMPLETE_SUMMARY.xlsx
│   └── stats/                    # Statistical test results
└── plasmid_outputs/              # Plasmid analysis results
    ├── P1-P5_*.png              # Plasmid visualizations
    └── plasmid_summary.xlsx
```

### `/scripts/` - Analysis Scripts

**Purpose:** Main analysis pipeline and modules

```
scripts/
├── 01_card_analysis.R            # Main AMR pipeline
├── 02_vfdb_analysis.R            # Main VFDB pipeline
├── 03_plasmid_analysis.R         # Main plasmid pipeline
├── modules/                      # Modular components
│   ├── amr/                      # AMR-specific modules
│   │   ├── 01_metadata_matching.R
│   │   ├── 02_data_loading.R
│   │   ├── 03_base_visualizations.R
│   │   ├── 04_enrichment_analysis.R
│   │   ├── 05_temporal_analysis.R
│   │   ├── 06_clinical_comparison.R
│   │   ├── 07_diversity_analysis.R
│   │   ├── 08_cooccurrence_network.R
│   │   ├── 09_statistical_models.R
│   │   └── 10_export_tables.R
│   ├── vfdb/                     # VFDB-specific modules
│   │   ├── 01_vfdb_loading.R
│   │   ├── 02_vf_class_analysis.R
│   │   ├── 03_prevalence_analysis.R
│   │   ├── 04_temporal_vf_trends.R
│   │   ├── 05_clinical_vf_comparison.R
│   │   ├── 06_diversity_metrics.R
│   │   ├── 07_statistical_tests.R
│   │   └── 08_export_tables.R
│   └── plasmid/                  # Plasmid-specific modules
│       ├── 01_plasmid_loading.R
│       ├── 02_arg_burden.R
│       ├── 03_highrisk_detection.R
│       ├── 04_vf_plasmid.R
│       └── 05_correlation_analysis.R
└── utils/                        # Shared utilities
    ├── common_utils.R            # Themes, colors, helpers
    └── metadata_matching.R       # ID matching functions
```

### `/tests/` - Unit Tests

**Purpose:** Automated testing (future implementation)

```
tests/
├── test_metadata_matching.R      # Test matching functions
├── test_data_loading.R           # Test data loading
├── test_statistical_functions.R # Test stats functions
└── test_integration.R            # End-to-end tests
```

## File Naming Conventions

### Scripts
- **Main pipelines:** `01_card_analysis.R`, `02_vfdb_analysis.R`, etc.
- **Modules:** Numbered by execution order (e.g., `01_metadata_matching.R`)
- **Utilities:** Descriptive names (e.g., `common_utils.R`)

### Plots
- **Metadata plots:** `M01`, `M02`, etc.
- **Base analysis:** `01`, `02`, etc.
- **Extended analysis:** Letter + number (e.g., `A1`, `B2`, `D6b`)
- **VF plots:** `V01`, `V02`, etc.
- **Plasmid plots:** `P1`, `P2`, etc.

### Tables
- **Summary files:** `*_SUMMARY.xlsx`
- **Publication tables:** `Publication_Tables_*.docx`
- **CSV exports:** Descriptive names matching content

## Module Descriptions

### AMR Modules

| Module | File | Purpose |
|--------|------|---------|
| 1 | `01_metadata_matching.R` | Multi-strategy ID reconciliation |
| 2 | `02_data_loading.R` | Load and validate CARD data |
| 3 | `03_base_visualizations.R` | Core AMR plots |
| 4 | `04_enrichment_analysis.R` | ST-specific gene enrichment |
| 5 | `05_temporal_analysis.R` | Time-series analysis |
| 6 | `06_clinical_comparison.R` | Clinical vs. non-clinical |
| 7 | `07_diversity_analysis.R` | Shannon, NMDS, PERMANOVA |
| 8 | `08_cooccurrence_network.R` | Gene co-occurrence |
| 9 | `09_statistical_models.R` | GLM, logistic regression |
| 10 | `10_export_tables.R` | Generate summary tables |

### VFDB Modules

| Module | File | Purpose |
|--------|------|---------|
| 1 | `01_vfdb_loading.R` | Load VF data with metadata |
| 2 | `02_vf_class_analysis.R` | VF class composition |
| 3 | `03_prevalence_analysis.R` | Gene prevalence profiling |
| 4 | `04_temporal_vf_trends.R` | Temporal dynamics |
| 5 | `05_clinical_vf_comparison.R` | Clinical enrichment |
| 6 | `06_diversity_metrics.R` | VF profile diversity |
| 7 | `07_statistical_tests.R` | Comprehensive testing |
| 8 | `08_export_tables.R` | Summary table export |

### Plasmid Modules

| Module | File | Purpose |
|--------|------|---------|
| 1 | `01_plasmid_loading.R` | Load merged plasmid data |
| 2 | `02_arg_burden.R` | ARG burden analysis |
| 3 | `03_highrisk_detection.R` | ESBL/carb/MCR detection |
| 4 | `04_vf_plasmid.R` | Plasmid VF analysis |
| 5 | `05_correlation_analysis.R` | ARG-VF correlation |

## Data Flow

```
Input Data
    ↓
Metadata Matching (Module 01)
    ↓
Data Loading & Validation (Module 02)
    ↓
    ├─→ Base Visualizations (Module 03)
    ├─→ Enrichment Analysis (Module 04)
    ├─→ Temporal Analysis (Module 05)
    ├─→ Clinical Comparison (Module 06)
    ├─→ Diversity Analysis (Module 07)
    ├─→ Network Analysis (Module 08)
    └─→ Statistical Models (Module 09)
    ↓
Export Tables & Reports (Module 10)
    ↓
Final Outputs (plots/, tables/)
```

## Configuration Files

| File | Purpose |
|------|---------|
| `.gitignore` | Exclude data/outputs from git |
| `requirements.txt` | R package dependencies |
| `install.R` | Automated installation |
| `LICENSE` | MIT License text |
| `README.md` | Main project documentation |
| `CHANGELOG.md` | Version history |
| `CONTRIBUTING.md` | Contribution guidelines |

## Best Practices

### Adding New Modules

1. Create module file in appropriate directory
2. Follow naming convention: `##_descriptive_name.R`
3. Include function documentation
4. Return results as named list
5. Update main pipeline script
6. Document in CHANGELOG.md

### Modifying Existing Modules

1. Test changes with example data
2. Update documentation if needed
3. Maintain backward compatibility
4. Update version in CHANGELOG.md

### Output Organization

- Keep plots and tables separate
- Use consistent file naming
- Include metadata in filenames when relevant
- Generate index/summary files for large output sets

## Size Estimates

Typical disk usage (approximate):

```
data/              500 MB - 5 GB (depends on genome count)
outputs/           100 MB - 500 MB
scripts/           1 MB
docs/              < 1 MB
examples/          < 10 MB
Total:             600 MB - 6 GB
```

## Dependencies

### System Requirements
- R ≥ 4.0.0
- RAM: 8 GB minimum, 16 GB recommended
- Disk space: 10 GB free (for data + outputs)
- OS: Linux, macOS, or Windows

### R Package Count
- Core: 30+ packages
- All dependencies: 100+ packages (auto-installed)

## Version Control

### What to Commit
- All scripts and modules
- Documentation
- Example data (small files only)
- Configuration files

### What to Ignore
- Actual data files (too large)
- Generated outputs
- Temporary files
- User-specific configurations

## Future Structure Additions

Planned for future versions:

```
config/                           # Configuration files
├── default_config.yaml          # Default parameters
└── user_config.yaml             # User overrides

docker/                          # Docker containerization
├── Dockerfile
└── docker-compose.yml

workflows/                       # Workflow definitions
├── snakemake/
└── nextflow/

reports/                         # Generated reports
├── html/
└── pdf/
```
