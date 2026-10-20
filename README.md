# Evaluating Graph Neural Networks for Structure-Aware Aqueous Solubility Prediction Under Molecular Distribution Shift

## Project overview

`molecular-property-gnn` is a research project investigating aqueous solubility
prediction, beginning with the Delaney ESOL dataset. Its objective is to compare
traditional molecular-descriptor/fingerprint machine-learning models with Graph
Neural Networks (GNNs) under a shared, structure-aware evaluation protocol.

The research will examine how predictive performance changes with molecular
structural familiarity, followed by explainability and uncertainty estimation in
later milestones. GNN superiority is a research hypothesis, not an established
result. The primary evaluation uses a hybrid structure-aware split combining
Bemis-Murcko scaffolds and acyclic similarity-connected components.

See [research questions and hypotheses](docs/research_questions.md) for the accepted
research framing and unresolved statistical-design decisions.

## Current implementation status

The implementation covers **dataset validation, structural diagnostics, and the approved
CPU-based hybrid split generator**. The split generator is checked using synthetic
fixtures; real ESOL split generation and artifact acceptance require separate approval.
No real split, model features, or models have been created by this milestone.

## Dataset source

This project currently uses the Delaney ESOL aqueous-solubility dataset.

- **Dataset file:** `data/raw/ESOL_delaney-processed.csv`.
- **Molecular records:** 1,128; this is a record count, not a unique-structure count.
- **Structure field:** `smiles`.
- **Target:** measured log aqueous solubility, `logS = log10(S / (mol/L))`, where
  `S` is aqueous solubility. The exact target column is
  `measured log solubility in mols per litre`.
- **Download source (owner-provided):**
  [GLambard/Molecules_Dataset_Collection — ESOL_delaney-processed.csv](https://github.com/GLambard/Molecules_Dataset_Collection/blob/master/originals/ESOL_delaney-processed.csv).
- **Original reference:** Delaney, J. S. (2004). *ESOL: Estimating Aqueous
  Solubility Directly from Molecular Structure*. Journal of Chemical Information
  and Computer Sciences, 44(3), 1000–1005.
  [DOI: 10.1021/ci034243x](https://doi.org/10.1021/ci034243x).

The download source above was supplied by the owner; its upstream revision and
byte-for-byte correspondence with the local CSV have not been independently
verified. See [dataset documentation](docs/dataset.md) for the recorded local
checksum, validation findings, provenance limitations, and redistribution caveats.

## Repository structure

```text
molecular-property-gnn/
├── AGENTS.md                      # Governance, approval, and review rules
├── README.md                      # Project overview and usage
├── pyproject.toml                 # Dependencies and tool configuration
├── uv.lock                        # Locked dependency resolution
├── .gitignore                     # Local data, environment, and artifact exclusions
├── docs/
│   ├── dataset.md                 # Dataset contract, provenance, and limitations
│   ├── primary_split.md           # Accepted primary split specification
│   └── research_questions.md      # Accepted research questions and hypotheses
├── data/                          # Local datasets; excluded from Git
│   ├── raw/                       # Immutable source data
│   ├── interim/                   # Existing structural audit evidence
│   └── processed/                 # Reserved for approved derived artifacts
├── src/
│   ├── data/
│   │   ├── audit.py               # Dataset validation and initial audit
│   │   ├── structural_audit.py    # Label-blind molecular structure audit
│   │   ├── acyclic_diagnostic.py  # Acyclic similarity diagnostics
│   │   └── hybrid_split.py        # CPU-based candidate split generator
│   ├── features/                  # Reserved for future feature implementation
│   ├── models/                    # Reserved for future models
│   ├── training/                  # Reserved for future training code
│   ├── evaluation/                # Reserved for future evaluation code
│   └── utils/                     # Reserved for future shared utilities
├── tests/                         # Synthetic validation and structural tests
├── configs/                       # Reserved for future configurations
├── notebooks/                     # Reserved for future exploratory notebooks
└── outputs/                       # Reserved for future experiment outputs
    ├── figures/
    ├── metrics/
    └── models/
```

The tree summarizes the existing local layout and omits environment/cache folders,
Python package markers, and individual data artifacts. Reserved folders are not
implemented capabilities; empty directories may not appear in a fresh Git clone.

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
