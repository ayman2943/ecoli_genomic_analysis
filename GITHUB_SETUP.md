# GitHub setup

Instructions for publishing this repository on GitHub.

## Prerequisites

- A GitHub account
- Git installed locally

## 1. Local repo

This directory is already a git repository with the full pipeline committed
(stage files with `git add -A`, then `git commit`).

## 2. Create the remote repository

Create a new (empty) GitHub repository named `ecoli_genomic_analysis` — either
at https://github.com/new or with:

```bash
# optional: install gh CLI
gh repo create ecoli_genomic_analysis --public --source . --remote origin --push
```

Do **not** initialise it with a README, `.gitignore` or licence (this repo
already contains them).

## 3. Link and push

If the repo was created on the web (not via `gh`), link and push manually:

```bash
git remote add origin https://github.com/<username>/ecoli_genomic_analysis.git
git branch -M main
git push -u origin main
```

The current remote already points at `origin
https://github.com/ayman2943/ecoli_genomic_analysis.git`. If that repository
was deleted/re-created, first remove the stale remote:

```bash
git remote remove origin
git remote add origin https://github.com/<username>/ecoli_genomic_analysis.git
```

## 4. After pushing

- Check **Settings → Secrets** (or use `.env` locally) so `NCBI_API_KEY` is
  never committed — the downloader reads it from the environment.
- Add a description and topic tags (e.g. `ecoli`, `pangenomics`, `amr`,
  `virulence`) on the repository page.
