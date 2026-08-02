# GitHub upload steps

## Option A: Web upload

1. Create a new GitHub repository.
2. Recommended name: `caz-avi-evidence-composite-pkpd`.
3. Keep it **private** until journal submission policy and file contents are confirmed.
4. Upload the files and folders from this repository package.
5. Add a short description: `Provenance-aware evidence-composite workflow for ceftazidime-avibactam PK/PD simulation`.

## Option B: Command line

```bash
git init
git add .
git commit -m "Initial CAZ-AVI evidence-composite reproducibility package"
git branch -M main
git remote add origin https://github.com/<USER>/caz-avi-evidence-composite-pkpd.git
git push -u origin main
```

## Before making public

- Confirm whether the journal permits posting the submitted manuscript PDF.
- Remove any copyrighted article PDFs. This package intentionally does not redistribute source article PDFs.
- Confirm whether the supplementary workbook should be public at submission or after acceptance.
- Add DOI after Zenodo archiving, if desired.
