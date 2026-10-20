# molecular-property-gnn

Research foundation for aqueous solubility prediction. The eventual goal is to compare
traditional molecular-descriptor machine learning with graph neural networks (GNNs),
then investigate scaffold-based evaluation, explainability, and uncertainty estimation.

The current implementation covers dataset validation and an initial data audit.
No dataset splits, model features, or models are implemented.

## Setup and use

Run from the repository root in WSL/Linux with Python 3.12 and
[uv](https://docs.astral.sh/uv/):

```bash
uv sync --extra dev --locked
.venv/bin/python -m src.data.audit data/raw/ESOL_delaney-processed.csv
.venv/bin/python -m pytest
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
```

The existing folders are preserved. Reusable code lives in `src/data/`; tests use
small synthetic fixtures in `tests/`. Raw data is excluded from Git by the existing
`.gitignore`. Place the supplied CSV at the path above; the audit never downloads,
rewrites, repairs, or filters it.

The command prints a JSON report to standard output. Exit codes are `0` for passing
content checks, `1` for reported content problems, and `2` for unreadable input or
schema/CSV errors. Whitespace and duplicates are reported findings, not automatic
failures. A pass indicates only that the documented validation checks succeeded.

To save a new audit snapshot outside the raw-data directory:

```bash
.venv/bin/python -m src.data.audit data/raw/ESOL_delaney-processed.csv > docs/esol_audit.json
```

See [dataset documentation](docs/dataset.md) for exact columns, target definition,
provenance, validation rules, initial findings, and limitations. A local
machine-readable audit snapshot is retained at `docs/esol_audit.json`, with
the input checksum and runtime versions. It is excluded from Git pending
separate publication approval.
