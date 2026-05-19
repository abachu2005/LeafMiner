# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.1] - 2026-05-19

### Added
- `.zenodo.json` for automatic Zenodo archival + DOI minting on every release.

## [1.0.0] - 2026-05-19

First public release suitable for general OSS distribution.

### Added
- 4-step in-browser **Setup Wizard** overlay (mode → config → input source → finish),
  auto-shown on first visit, re-openable from a button at the top.
- `bin/leafcutter-setup` CLI wizard (venv, deps, Quest profile to
  `~/.leafcutter/config.json` with optional SSH test, smoke test).
- `tests/fixtures/sample.SJ.out.tab` for the wizard smoke test.
- `LICENSE` (MIT), `CITATION.cff`, `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`.
- GitHub Actions CI: pytest + ruff on Python 3.9–3.12.
- Dockerfile + GHCR publish workflow.
- pytest tests for the LC2 pipeline helpers and webapp surface.

### Changed
- Default execution mode is now **Local** (Quest Slurm still available).
- Webapp and CLI no longer carry hardcoded NetID / Slurm account / Quest paths;
  the wizard collects them per-user and persists locally.
- `recount3_to_bed.py` / `recount3_to_star_sj.py` read `$LEAFCUTTER_ROOT`
  instead of hardcoded `/Users/...` paths.

### Removed
- Personal/collaborator artifacts: `email_work_summary.*`,
  `scripts/aggregate_for_saba.py`, `scripts/generate_summary_for_saba.py`,
  `scripts/plot.py`, `scripts/plot2.py`, `scripts/build_file_list.py`,
  `scripts/recount3_to_bed_TEST.py`.

[Unreleased]: https://github.com/abachu2005/Leaf_Cutter/compare/v1.0.1...HEAD
[1.0.1]: https://github.com/abachu2005/Leaf_Cutter/releases/tag/v1.0.1
[1.0.0]: https://github.com/abachu2005/Leaf_Cutter/releases/tag/v1.0.0
