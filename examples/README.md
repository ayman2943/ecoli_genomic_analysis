# Example Data

This directory contains example datasets to help you understand the required data format and test the pipeline.

## Directory Structure

```
examples/
├── data/
│   ├── card_summary/
│   │   ├── ST10/
│   │   │   └── Escherichia_coli_example1_card.tsv
│   │   └── ST131/
│   │       └── Escherichia_coli_example2_card.tsv
│   ├── vfdb_summary/
│   │   ├── ST10/
│   │   │   └── Escherichia_coli_example1_vfdb.tsv
│   │   └── ST131/
│   │       └── Escherichia_coli_example2_vfdb.tsv
│   ├── metadata/
│   │   ├── ST10_filtered.xlsx
│   │   └── ST131_filtered.xlsx
│   ├── plasmid_CARD_merged.tsv
│   └── plasmid_VFDB_merged.tsv
└── README.md (this file)
```

## File Formats

### CARD Summary Files

Tab-separated files with ARG presence/absence or counts.

**File naming:** `Escherichia_coli_{GENOME_ID}_card.tsv`

**Format:**
```tsv
aac(3)-IIa	aac(6')-Ib-cr	blaCTX-M-15	sul1	tet(A)
1	0	1	1	0
```

### VFDB Summary Files

Tab-separated files with VF gene presence/absence or counts.

**File naming:** `Escherichia_coli_{GENOME_ID}_vfdb.tsv`

**Format:**
```tsv
fimH	iucD	ompA	papC	traT
1	1	1	0	1
```

### Metadata Files

Excel files with genome metadata.

**Required columns:**
- `Name` - Must match genome IDs in summary files
- `Collection Year` - Numeric year (e.g., 2018)
- `Country` - Country of isolation
- `Continent` - Geographic continent
- `Source Niche` - Source category (e.g., Human_Clinical, Animal, Environmental)

**Example:**
| Name | Collection Year | Country | Continent | Source Niche |
|------|----------------|---------|-----------|--------------|
| example1 | 2018 | USA | North America | Human_Clinical |
| example2 | 2019 | UK | Europe | Animal |

### Plasmid Files

Merged plasmid data for all STs in tab-separated format.

**Required columns:**
- `Name` - Genome identifier
- `ST` - Sequence type
- `Collection Year` - Numeric year
- `Country`, `Continent` - Geographic info
- Gene columns - ARG or VF presence/absence

## Running Example Analysis

### Quick Test

```bash
# Navigate to project root
cd ecoli-genomic-analysis

# Copy example data to main data directory
cp -r examples/data/* data/

# Run analyses
Rscript scripts/01_card_analysis.R
Rscript scripts/02_vfdb_analysis.R
Rscript scripts/03_plasmid_analysis.R
```

### Expected Outputs

After running the example, you should see:

```
outputs/
├── R_analysis_outputs/
│   ├── plots/
│   │   ├── M01_sample_distribution.png
│   │   ├── 01_ARG_burden_by_ST.png
│   │   └── ... (50+ plots)
│   └── tables/
│       ├── AMR_COMPLETE_SUMMARY.xlsx
│       └── ... (CSV tables)
├── VFDB_analysis_outputs/
│   └── ... (VF plots and tables)
└── plasmid_outputs/
    └── ... (plasmid plots)
```

## Modifying for Your Data

1. **Replace example files** with your actual data
2. **Ensure file naming** follows the pattern shown above
3. **Match genome IDs** between summary files and metadata
4. **Update sequence types** if using different STs:
   ```r
   # In main scripts
   STs <- c("ST10", "ST131", "ST69", "ST73", "ST95")
   ```

## Data Quality Checklist

Before running analysis, verify:

- [ ] All genome IDs match between CARD/VFDB files and metadata
- [ ] Metadata has required columns with correct names
- [ ] No missing values in critical fields (ST, Year, Country)
- [ ] File encodings are UTF-8
- [ ] Tab-separated files use actual tabs (not spaces)
- [ ] Excel files are not corrupted
- [ ] Directory structure matches expected format

## Troubleshooting

**Issue:** Low metadata matching rate

**Solution:** Check genome ID formats. The pipeline tries 13 different matching strategies, but you may need to standardize naming.

**Issue:** Missing gene columns

**Solution:** Ensure gene names are consistent across all files. The pipeline expects exact matches.

**Issue:** Excel file errors

**Solution:** Verify Excel files open correctly and contain no merged cells or special formatting.

## Need Help?

- Review [docs/USAGE.md](../docs/USAGE.md) for detailed instructions
- Check [GitHub Issues](https://github.com/yourusername/ecoli-genomic-analysis/issues)
- Contact maintainers

## Notes

The example data provided is synthetic and for demonstration purposes only. It does not represent real biological data and should not be used for research conclusions.
