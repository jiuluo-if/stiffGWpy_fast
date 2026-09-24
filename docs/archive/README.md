# Historical documentation and local experiment data

[中文说明](README_zh.md) | English

This directory holds material retained for provenance after it stopped serving
as a current report or active navigation target.

- `baselines/baseline_54d65e3.md` is the 2026-09-11 two-thread snapshot. The
  latest detailed accuracy audit is a dated 2026-09-12 record in
  `../fast_v02_audit_report.md`; the `../validation/validation_manifest.json`
  is dated 2026-09-17. Neither is a current-HEAD full-space certification.
- `oracle_prufer/raw_logs/` is a local-only holding area for three 2026-09-11
  diagnostic captures moved out of the active `docs/` root. The repository's
  `*.log` ignore rule keeps these raw captures out of Git; later JSON records
  and the assessment remain the shared, navigable evidence. A fresh clone will
  not contain the ignored raw logs.
- `reference/head_20260910/` contains an older benchmark reference snapshot.

Archived material is historical evidence. It must not be read as a current
production claim without checking its recorded commit and protocol.
