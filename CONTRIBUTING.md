# Contributing to E. coli Genomic Analysis Pipeline

First off, thank you for considering contributing to this project! 🎉

## Code of Conduct

This project and everyone participating in it is governed by our Code of Conduct. By participating, you are expected to uphold this code.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check existing issues to avoid duplicates. When you create a bug report, include as many details as possible:

- **Use a clear and descriptive title**
- **Describe the exact steps to reproduce the problem**
- **Provide specific examples** (sample data, code snippets)
- **Describe the behavior you observed** and what you expected
- **Include R session info** (`sessionInfo()`)
- **Include error messages** (complete stack traces)

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion:

- **Use a clear and descriptive title**
- **Provide a detailed description** of the suggested enhancement
- **Explain why this enhancement would be useful**
- **List any alternative solutions** you've considered

### Pull Requests

1. Fork the repository
2. Create a new branch from `main`
3. Make your changes
4. Add or update tests as appropriate
5. Update documentation
6. Ensure all tests pass
7. Submit a pull request

#### Pull Request Guidelines

- Follow the existing code style
- Write clear, descriptive commit messages
- Update the README.md if needed
- Add comments for complex logic
- Include examples if adding new features

## Development Setup

### Prerequisites

```bash
# Clone your fork
git clone https://github.com/yourusername/ecoli-genomic-analysis.git
cd ecoli-genomic-analysis

# Install dependencies
Rscript -e "install.packages(c('tidyverse', 'data.table', ...))"
```

### Code Style

- Use `tidyverse` style guide
- 2-space indentation
- Maximum line length: 80 characters
- Use meaningful variable names
- Add comments for complex operations

### Testing

Before submitting:

```r
# Run test data through pipeline
source("scripts/tests/test_pipeline.R")

# Check for warnings
devtools::check()
```

## Documentation

- Update README.md for major changes
- Add inline comments for complex functions
- Update CHANGELOG.md following [Keep a Changelog](https://keepachangelog.com/)

## Commit Messages

- Use present tense ("Add feature" not "Added feature")
- Use imperative mood ("Move cursor to..." not "Moves cursor to...")
- Limit first line to 72 characters
- Reference issues and pull requests when relevant

Examples:
```
Add temporal analysis for plasmid genes

- Implement Mann-Kendall trend test
- Add visualization for trends
- Update documentation

Closes #123
```

## Review Process

1. Maintainers will review your PR
2. Changes may be requested
3. Once approved, your PR will be merged
4. Your contribution will be acknowledged in CHANGELOG.md

## Recognition

Contributors will be:
- Listed in CONTRIBUTORS.md
- Acknowledged in release notes
- Credited in any publications using the enhanced pipeline

## Questions?

Feel free to open an issue with the `question` label or contact the maintainers.

Thank you for your contributions! 🙏
