# Changelog

All notable changes to the E. coli Genomic Analysis Pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-03-17

### Added
- Initial release of modular analysis pipeline
- AMR (CARD) analysis module with 10 sub-modules
- VFDB analysis module with 8 sub-modules
- Plasmid analysis module with 5 sub-modules
- Comprehensive metadata matching system (13 strategies)
- Publication-quality visualization themes
- Statistical testing framework (15+ tests)
- Multi-format export (Excel, CSV, DOCX, PNG)
- Parallel processing support
- Detailed documentation and usage guides

### Features

#### AMR Analysis
- Multi-strategy genome ID reconciliation
- Temporal trend analysis (Mann-Kendall tau)
- ST-specific gene enrichment (log2 fold-change, Fisher's exact)
- Drug class-specific analyses (aminoglycosides, β-lactams)
- Clinical vs. non-clinical comparisons (per-ST and overall)
- MDR profiling and classification
- Shannon diversity and NMDS ordination
- PERMANOVA for beta diversity
- Gene co-occurrence networks (Phi correlation, igraph)
- Negative binomial GLM for incidence rate ratios
- Logistic regression for high-burden predictors
- 50+ publication-quality visualizations
- 21-sheet Excel summary workbook

#### VFDB Analysis
- VF class composition and profiling
- Gene prevalence analysis (top 40 genes)
- Temporal dynamics of VF classes
- Clinical enrichment per ST
- Shannon diversity of VF profiles
- Co-occurrence correlation matrices
- Comprehensive statistical testing
- 15+ specialized visualizations
- Multi-panel combined figures

#### Plasmid Analysis
- Plasmid-borne ARG burden quantification
- MDR plasmid prevalence tracking
- High-risk gene detection (ESBL, carbapenemases, MCR)
- Plasmid VF profiling
- ARG-VF correlation analysis
- Temporal trend visualization
- 5 key analytical plots

#### Statistical Methods
- Kruskal-Wallis with Bonferroni-Dunn post-hoc
- Mann-Kendall trend tests per ST
- Linear models with Year × ST interaction
- Jaccard PERMANOVA / adonis2
- NMDS ordination with stress annotation
- Fisher's exact test with OR and 95% CI
- Log2FC enrichment with BH-FDR correction
- Spearman correlation analyses
- Negative binomial GLM
- Logistic regression
- Chi-squared heterogeneity tests
- Wilcoxon tests with BH-FDR correction

### Documentation
- Comprehensive README with installation guide
- Detailed USAGE.md with examples
- CONTRIBUTING.md for developers
- MIT License
- Code of Conduct (implicit in CONTRIBUTING.md)
- Example data structure diagrams

### Infrastructure
- Modular script organization
- Shared utility functions
- Consistent error handling
- Progress reporting
- Automated output directory creation
- Session info logging

---

## [Unreleased]

### Planned Features
- Support for additional sequence types
- Integration with phylogenetic trees
- Machine learning predictions (ARG/VF presence)
- Interactive HTML reports (Shiny dashboard)
- Docker containerization
- Snakemake workflow integration
- Additional databases (ResFinder, PlasmidFinder)
- Automated quality control checks
- Batch processing scripts
- Cloud deployment options (AWS, Google Cloud)

### Under Consideration
- GUI interface for non-programmers
- Real-time monitoring dashboard
- Integration with public databases (NCBI, ENA)
- Automated report generation
- Multi-language support
- Additional export formats (LaTeX, Markdown)

---

## Version History

### Version Numbering
- **Major version (X.0.0)**: Incompatible API changes
- **Minor version (0.X.0)**: New features, backwards compatible
- **Patch version (0.0.X)**: Bug fixes, backwards compatible

### Support Policy
- **Current version**: Full support with updates
- **Previous major version**: Security fixes only
- **Older versions**: No support (upgrade recommended)

---

## Migration Guides

### From Monolithic Scripts (Pre-1.0.0)

If migrating from the original combined scripts:

1. **Data Structure**: Ensure data follows new directory structure (see USAGE.md)
2. **Dependencies**: Install all packages from requirements.txt
3. **Configuration**: Update file paths in main scripts
4. **Execution**: Run individual pipeline scripts instead of monolithic code

**Breaking Changes:**
- Output directory structure reorganized
- Function names standardized
- Module-based execution required
- Configuration via script variables (no config file yet)

---

## Contributors

See [CONTRIBUTORS.md](CONTRIBUTORS.md) for list of contributors.

---

## Links

- [Repository](https://github.com/yourusername/ecoli-genomic-analysis)
- [Issues](https://github.com/yourusername/ecoli-genomic-analysis/issues)
- [Releases](https://github.com/yourusername/ecoli-genomic-analysis/releases)

---

**Note:** This changelog follows the format recommended by [Keep a Changelog](https://keepachangelog.com/).
