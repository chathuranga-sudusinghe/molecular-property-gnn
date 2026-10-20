# molecular-property-gnn

Research foundation for aqueous solubility prediction. The eventual goal is to compare
traditional molecular-descriptor machine learning with graph neural networks (GNNs),
then investigate scaffold-based evaluation, explainability, and uncertainty estimation.

The implementation covers **dataset validation, structural diagnostics, and the approved
CPU-based hybrid split generator**. The split generator is checked using synthetic
fixtures; real ESOL split generation and artifact acceptance require separate approval.
No real split, model features, or models have been created by this milestone.

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

## Primary structural split

See [the accepted primary split contract](docs/primary_split.md) for identity rules,
typed Murcko groups, acyclic similarity-connected components, the deterministic
allocator, artifact schemas, and the procedural test seal.

The generator uses fixed Morgan radius-2 sparse count fingerprints and the
project-specific inclusive Tanimoto threshold 0.70. It does not use targets or
historical predictions. It writes an unsealed candidate requiring owner acceptance.

The documented ESOL command requires separate real-data execution approval; do not
run it as part of synthetic implementation checks. GPU/CUDA handling is reserved for
a later model-training milestone. Existing audit evidence remains unchanged.
