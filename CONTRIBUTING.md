# Contributing to Leaf_Cutter

Thanks for your interest! Pipeline improvements, bug fixes, new tissue presets, and documentation updates are all welcome.

## Setup

```bash
git clone https://github.com/abachu2005/Leaf_Cutter.git
cd Leaf_Cutter
python3 bin/leafcutter-setup   # creates webapp/.venv, installs deps, smoke-tests
```

For development:

```bash
./webapp/.venv/bin/pip install -e ".[dev]"
./webapp/.venv/bin/pre-commit install
```

## Tests

```bash
./webapp/.venv/bin/pytest -q
```

CI runs the same tests on Python 3.9–3.12.

## Code style

Formatting + linting via [ruff](https://docs.astral.sh/ruff/) (configured in `pyproject.toml`). Pre-commit runs it automatically.

## Pull requests

1. Fork and branch off `main`.
2. Add tests (under `tests/`) and update `CHANGELOG.md` "Unreleased" section.
3. Make sure `pytest` and `ruff check .` pass.
4. Open a PR. Use the template.

## A note on Quest-specific code

This repo supports both general local users and Northwestern Quest HPC users. Please don't commit hardcoded NetIDs, account names, or `/projects/...` paths. Use the wizard / form fields / environment variables.

## Code of Conduct

Participation governed by [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).
