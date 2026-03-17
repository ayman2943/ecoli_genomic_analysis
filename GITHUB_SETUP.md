# GitHub Setup Instructions

Follow these steps to upload your repository to GitHub.

## Prerequisites

- Git installed on your computer
- GitHub account created
- Repository created on GitHub (empty, no README)

## Step-by-Step Upload

### 1. Initialize Git Repository (Local)

```bash
cd ecoli-genomic-analysis

# Initialize git
git init

# Add all files
git add .

# Make first commit
git commit -m "Initial commit: E. coli Genomic Analysis Pipeline v1.0.0"
```

### 2. Connect to GitHub

```bash
# Add your GitHub repository as remote
# Replace 'yourusername' with your actual GitHub username
git remote add origin https://github.com/yourusername/ecoli-genomic-analysis.git

# Verify remote
git remote -v
```

### 3. Push to GitHub

```bash
# Push to main branch
git branch -M main
git push -u origin main
```

## Alternative: Using GitHub Desktop

1. Open GitHub Desktop
2. File → Add Local Repository
3. Select the `ecoli-genomic-analysis` folder
4. Click "Publish repository"
5. Choose repository name and visibility
6. Click "Publish repository"

## Post-Upload Tasks

### Add Topics/Tags

On GitHub, add these topics to help others find your repository:

- `r`
- `bioinformatics`
- `genomics`
- `e-coli`
- `antimicrobial-resistance`
- `amr`
- `virulence-factors`
- `card`
- `vfdb`
- `data-analysis`
- `pipeline`

### Enable GitHub Pages (Optional)

For hosting documentation:

1. Settings → Pages
2. Source: Deploy from branch
3. Branch: main / docs folder
4. Save

### Set Up Issues

1. Go to Settings → General
2. Enable "Issues"
3. Consider creating issue templates

### Add Branch Protection (Recommended)

1. Settings → Branches
2. Add branch protection rule for `main`
3. Require pull request reviews
4. Require status checks

## Repository Settings

### About Section

Add to your repository's About section:

**Description:**
```
Comprehensive R pipeline for E. coli genomic analysis: AMR (CARD), virulence factors (VFDB), and plasmid-associated genes
```

**Website:**
```
https://yourusername.github.io/ecoli-genomic-analysis
```

**Topics:**
```
r, bioinformatics, genomics, e-coli, antimicrobial-resistance, amr, virulence-factors
```

### README Badges

Your README already includes these badges:
- License badge
- R version badge

Consider adding:
- ![GitHub issues](https://img.shields.io/github/issues/yourusername/ecoli-genomic-analysis)
- ![GitHub stars](https://img.shields.io/github/stars/yourusername/ecoli-genomic-analysis)
- ![GitHub forks](https://img.shields.io/github/forks/yourusername/ecoli-genomic-analysis)

## Create First Release

### 1. Create a Git Tag

```bash
git tag -a v1.0.0 -m "Release version 1.0.0"
git push origin v1.0.0
```

### 2. Create GitHub Release

1. Go to Releases → Draft a new release
2. Tag version: `v1.0.0`
3. Release title: `v1.0.0 - Initial Release`
4. Description: Copy from CHANGELOG.md
5. Attach any binary files (if applicable)
6. Publish release

## Continuous Integration (Future)

Consider setting up GitHub Actions for:

- Automated testing
- Code style checking
- Documentation building

Example `.github/workflows/test.yml`:

```yaml
name: R Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: r-lib/actions/setup-r@v2
      - name: Install dependencies
        run: Rscript install.R
      - name: Run tests
        run: Rscript -e 'source("tests/test_all.R")'
```

## Maintenance Schedule

- **Weekly**: Review and respond to issues
- **Monthly**: Update dependencies
- **Quarterly**: Major feature releases
- **Annually**: Major version updates

## Troubleshooting

### Large Files

If you have files >100MB:

```bash
# Use Git LFS
git lfs install
git lfs track "*.xlsx"
git lfs track "*.png"
git add .gitattributes
git commit -m "Add Git LFS"
```

### Authentication Issues

Use personal access token instead of password:

1. GitHub → Settings → Developer settings
2. Personal access tokens → Generate new token
3. Select scopes: repo, workflow
4. Use token as password when pushing

### Rename Repository

If you need to rename later:

1. GitHub → Settings → Repository name
2. Update local remote:
   ```bash
   git remote set-url origin https://github.com/yourusername/new-name.git
   ```

## Best Practices

1. ✅ Write clear commit messages
2. ✅ Use branches for new features
3. ✅ Keep main branch stable
4. ✅ Review pull requests thoroughly
5. ✅ Update CHANGELOG.md for each release
6. ✅ Tag releases with semantic versioning
7. ✅ Respond to issues promptly
8. ✅ Keep documentation up-to-date

## Getting Help

- [GitHub Docs](https://docs.github.com)
- [Git Documentation](https://git-scm.com/doc)
- [GitHub Community](https://github.community)

---

**After upload, your repository will be live at:**
`https://github.com/yourusername/ecoli-genomic-analysis`

Update the URLs in README.md and other files accordingly!
