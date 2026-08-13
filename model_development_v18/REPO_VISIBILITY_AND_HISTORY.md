# Why this repository's history begins on 12 August 2026

The git history here starts abruptly. That is deliberate and worth explaining, because a
reproducibility repository with no past invites the question.

## What happened

During development, two files containing **unpublished manuscript prose** were pushed to
this repository: a draft of a separate manuscript, and a working note that quoted passages
of the submitted article verbatim. Neither belongs in a public code-and-data repository
before publication.

Removing them from the current commit was not enough. Git keeps every version, and on
GitHub the commits referenced by a pull request stay reachable by their SHA even after the
branch is rewritten — so the text remained retrievable through the pull-request refs, which
cannot be deleted through the interface.

The repository was therefore **rebuilt from scratch on 12 August 2026**: deleted and
recreated under the same name, with a fresh history containing no manuscript prose in any
commit. A full mirror of the previous repository, including all pull requests and their
discussion, is retained offline by the author.

## What did not change

- **The URL is unchanged**, so the Data availability and Code availability statements in
  the article resolve exactly as published.
- **The file content is unchanged.** The republished tree is byte-identical to the previous
  one — the same 340 files, and `SHA256SUMS.txt` still verifies.
- Nothing was removed from the scientific record. What was removed was manuscript text that
  had never been part of it.

## What this repository does and does not contain

**Contains:** processed simulation inputs, output workbooks, provenance registers, runtime
metadata, checksums, analysis code and figure-generation scripts, the external evidence base
with per-source provenance, and the model-development work in `model_development_v18/`.

**Does not contain:** manuscript text, and PDFs of source articles, which are not
redistributed for copyright reasons. Each entry under `data_external/` records the DOI so
the source can be obtained through the reader's own institutional access.
