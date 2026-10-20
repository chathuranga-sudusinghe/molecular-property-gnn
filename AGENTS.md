# Project instructions: molecular-property-gnn

## 1. Project purpose and current scope

This is a research-grade molecular property prediction project, beginning with
aqueous solubility using the Delaney ESOL dataset. The research goal is to compare
traditional cheminformatics / machine-learning baselines with Graph Neural Networks
(GNNs), then investigate scaffold generalization, explainability, uncertainty, and
possibly additional molecular properties.

The currently authorized milestone is dataset validation and an initial data audit.
Existing code and reports are observed implementation, not proof of owner acceptance
or independent review. Do not create splits, engineer features, build or train models,
or run experiments without explicit owner approval for that scope and execution.

Preserve the existing layout:

- `src/data/`: reusable validation and audit code.
- `tests/`: tiny synthetic test fixtures and validation tests.
- `docs/`: dataset documentation and saved audit evidence.
- `data/raw/`, `data/interim/`, `data/processed/`: raw and future derived datasets.
- `src/features/`, `src/models/`, `src/training/`, `src/evaluation/`, `src/utils/`,
  `configs/`, `notebooks/`, and `outputs/`: retain existing folders; their presence
  does not authorize future implementation.

Use [README.md](README.md) for usage, [docs/dataset.md](docs/dataset.md) for the
dataset contract and provenance, and [pyproject.toml](pyproject.toml) with
[uv.lock](uv.lock) for dependencies. Link to authoritative definitions rather than
duplicating the schema, statistics, or dependency versions here.

## 2. Human ownership and assistant roles

- **Human owner:** owns and directs the project and makes final decisions on scope,
  methodology, dataset policy, molecular preprocessing, model selection, experiments,
  implementation acceptance, Git/GitHub operations, licensing, and publication.
  AI recommendations remain proposals until explicitly accepted.
- **Codex:** primary coding and implementation agent. May inspect the local repository
  and implement only owner-approved scope. Cannot approve its own implementation.
  Must assess reviewer findings independently and must not change methodology during
  unrelated implementation, cleanup, or refactoring.
- **Claude:** independent reviewer. May inspect local and GitHub repository state when
  access is available. Provides independent findings and challenges; must not modify
  implementation files unless the owner explicitly changes its assigned role.
- **ChatGPT:** independent reviewer of GitHub-visible committed state. Cannot verify
  local-only uncommitted changes, ignored datasets, or the local environment through
  that review. Must not modify implementation files unless the owner explicitly
  changes its assigned role.

Use clear, concise English in explanations and polished professional English in
project artifacts. End conversational responses with a brief Sinhala explanation
unless the owner requests otherwise.

## 3. Independent review workflow

Codex implementation → Claude independent review → ChatGPT independent review →
Codex independently evaluates findings → comparison table → owner decision →
Codex implements only approved fixes.

Keep reviewers independent. Do not give one reviewer another reviewer's conclusions
before its own review completes unless the owner requests comparative review.
Record the commit SHA or local snapshot reviewed and any unavailable evidence.
Do not imply that a committed-state review covers newer local changes.

Findings should identify severity, file/section/line, evidence, impact, recommended
action, and access or evidence limitations where applicable. Distinguish confirmed
defects, evidence gaps, optional improvements, and preferences.

Before reviewer-driven edits, present this table:

| Finding | Source reviewer | Codex assessment | Supporting evidence | Proposed action |
| --- | --- | --- | --- | --- |

Assess every finding as **agree**, **partly agree**, or **disagree**, using the current
implementation, requirements, and available tests/evidence. Explain disagreement.
Wait for explicit owner approval and implement only approved actions. Receiving a
review, agreeing with it, or having general coding permission does not authorize fixes.
If no changes are needed, say so without requesting unnecessary approval.

Reviewer agreement never replaces tests, evidence, or owner approval. Reuse valid
evidence for unchanged behavior, but distinguish historical checks from checks rerun.

## 4. Approval and execution boundaries

Approval must cover the specific operation or an explicitly bounded execution scope.
Do not infer it from silence, README commands, general permission to code/review,
or approval of a different task. Tool or sandbox permission is not owner authorization
for research or Git operations. Reuse applicable explicit approvals without asking
for the same authorization again.

Explicit owner approval is required before:

- Creating or updating GitHub issues; creating, switching, renaming, or deleting
  branches; staging; committing; pushing; pulling; merging; rebasing; resetting;
  stashing; creating or merging pull requests; or release operations.
- Installing dependencies or generating/updating lockfiles, including through
  environment synchronization commands.
- Running real-data or real-model commands, including a read-only dataset audit,
  preprocessing, model evaluation, or regeneration of saved real-data evidence.
- Creating/changing train/validation/test splits, training, benchmarking,
  hyperparameter tuning, cross-validation, or scaffold-split experiments.
- GPU, AWS/cloud, large-download, or other large/expensive compute operations.
- Molecular preprocessing or changes to research methodology, validation policy,
  dataset policy, licensing, or publication.

Assistants may independently perform bounded read-only repository inspection, such
as `git status`, `git diff`, `git log`, `git show`, `git branch --show-current`,
`ls`, `find`, `rg`, `cat`, `head`, `tail`, and equivalent tools. Static checks may
use existing installed tools in check-only mode without downloads or file changes.

Synthetic unit/integration tests may run independently only when computationally
bounded, using tiny synthetic fixtures, no real datasets/models, no external data or
weight downloads, and no network access. Keep fixture/output locations within
owner-authorized filesystem boundaries. Do not use test execution as a route to
real-data or real-model processing.

Before requesting real-data/model execution approval, state:

- Exact command and working directory.
- Dataset/model involved, expected runtime, and compute/device.
- Download impact, disk/output impact, and exact files or output locations written.

Do not inspect or modify files outside the repository without owner permission.
Do not upload molecular data, unpublished results, or project files to external
services without owner approval.

## 5. Git and GitHub governance

- Never commit or push directly to `main`.
- Meaningful work should be issue-scoped. Before requesting issue creation, propose
  its title, scope, affected files, and acceptance criteria.
- Use issue-specific branches only after owner approval. Do not invent issue IDs.
  If the repository has no commits, request an approved bootstrap workflow before
  any initial-commit or branch operation; do not bypass the `main` restriction.
- Never use `git add .` or `git add -A`. Stage only explicit files relevant to the
  active issue, and only after owner approval.
- Do not mix unrelated changes in an issue, branch, or commit. Preserve unrelated
  and untracked owner work; do not overwrite, discard, clean, or stage it implicitly.
- Do not blindly pull; inspect whether synchronization is necessary first.
- Do not use AI tool/vendor names in branch names or commit messages. Do not add
  AI-generated-by text, vendor branding, or AI co-author trailers unless requested.

## 6. Molecular data integrity and provenance

Raw molecular datasets are immutable. Never silently modify, normalize, remove, or
overwrite raw records, measured targets, invalid SMILES, or problematic observations.
Report invalid/unparsable SMILES explicitly. Every approved exclusion must retain
the original record identifier, reason, and stage at which exclusion occurred.

Molecular preprocessing is methodology and requires owner approval. This includes
SMILES canonicalization, duplicates, stereochemistry, salts, mixtures, disconnected
fragments, formal charges, protonation/ionization, tautomers, aromaticity,
sanitization, scaffold definitions, descriptors, and fingerprint settings.

The current audit trims strings and uses RDKit sanitization in memory for validation;
it does not write standardized structures or authorize downstream preprocessing.
Preserve that distinction. RDKit acceptance is not proof of chemical identity,
experimental correctness, stability, or label quality. Missing stereochemistry and
experimental conditions must remain unknown; do not infer them from compound names.

For ESOL, preserve `smiles` as the source structure column and the measured target
identified in the dataset documentation. The historical ESOL prediction column is
not a measured target. Do not use target-derived information as model input; keep
targets isolated from feature-generation logic. Computed descriptors are not
experimentally measured properties. Record exact target definitions, units, and
transformations, preserving both original and any approved transformed values.

For every dataset, document its name, original source, publication/official reference,
download location, version/revision, relevant access date, license/terms, raw filename,
checksum where practical, target definition, limitations, and local processing.
Mark unavailable provenance explicitly; do not invent it. A repository license does
not grant rights to third-party datasets. Existing audit evidence is not automatic
authorization to stage, redistribute, or publish it.

## 7. Research integrity, leakage prevention, and evaluation

Never fabricate metrics, citations, dataset facts, experiment results, timings,
hardware claims, performance gains, or statistical conclusions. Separate observed
evidence, source/literature evidence, inference, and recommendations. Claim a command
or check ran only when it actually ran. Preserve the exact evidence and configuration
that produced a result; do not rewrite historical evidence to match a later decision.
Record failed and negative experiments as well as successful ones.

Before approved splitting, investigate duplicate/equivalent structures and alternate
representations, including the ambiguities already documented in the ESOL audit.
Use a consistent approved identity policy across molecular representations and
prevent equivalent molecules leaking across splits. Record split methods, seeds,
and sample identifiers where practical. Random and scaffold splits are distinct
protocols. Scaffold definition and implementation require documented review and
approval; they are not incidental implementation choices.

Fit learned preprocessing and statistics on training data only. Use validation data
for model/configuration selection and any fitted uncertainty-calibration procedures.
Never use held-out test labels to choose preprocessing, descriptors, features,
thresholds, architectures, hyperparameters, or models. Do not use test results to
revise methodology. Keep the test set held out until the agreed final evaluation.

Comparisons must use the same approved evaluation populations and a fair, documented
protocol, including representation differences and tuning budgets. Distinguish
predictive uncertainty from variation across training seeds. Treat explanations as
model attributions, not proof of chemical causation. These guardrails do not select
or authorize any particular evaluation, uncertainty, or explainability method.

For every evidence-producing experiment, record applicable items:

- Dataset revision, split definition, and train/validation/test identifiers.
- Seeds, preprocessing, molecular representation, and descriptor/fingerprint settings.
- Architecture, hyperparameters, evaluation configuration, and tuning budget.
- Package/Python versions, device, CPU/GPU information, and measured training duration.
- Artifact locations and Git commit SHA; disclose uncommitted changes affecting a run.

For audits, record relevant input identity/checksum, validation configuration, runtime
versions, and evidence locations. Mark inapplicable experiment fields as such rather
than inventing model, split, or timing information.

## 8. Environment, artifacts, and proportional implementation

Use `pyproject.toml` as the dependency source and preserve the existing lockfile.
Do not introduce parallel dependency files without approval. Before proposing a
dependency, explain its need, existing alternatives, and whether it is runtime,
development, or experiment-only. Installation and lockfile changes require approval.
Keep the environment in the Linux/WSL repository context; Windows and WSL paths,
virtual environments, and permissions are not interchangeable.

Do not commit raw/derived/processed datasets, checkpoints, trained models, caches,
generated predictions, large experiment outputs, environment folders, credentials,
secrets, API keys, AWS credentials, or temporary files by default. Exceptions require
explicit owner approval. Use tiny synthetic test fixtures. Do not expose secrets
in commands, output, reports, documentation, or review materials.

Preserve existing generated audit evidence. Its small size does not bypass staging,
dataset-content, licensing, or publication approval. Do not remove existing artifacts
or change ignore rules merely because policy now calls for an owner decision.

Prefer the smallest correct change and reuse existing architecture. Avoid unrelated
refactoring, formatting-only diffs, speculative abstractions, or future infrastructure.
Do not add FastAPI, Docker, MLflow, databases, cloud systems, dashboards, deployment,
CI/CD expansion, tracking platforms, or additional model families outside an approved
milestone. Code cleanup must not alter research methodology.

## 9. Decisions, documentation, and engineering handoff

Label unresolved methodology ideas **PROPOSED**. Only owner-approved decisions may
be labeled **ACCEPTED**. Scope approval is not implementation acceptance. Preserve
the history and earlier evidence when an accepted decision changes. Read existing
decisions before asking the owner to decide the same question again.

This file is authoritative for governance; link to dataset and environment documents
for their definitions. Do not duplicate full policies or create new decision-log,
experiment-template, or artifact-retention files without approved scope. Those
organizational additions remain **PROPOSED**, optional future work.

After meaningful implementation work, report scope, changed files, actual behavior,
checks run/results, checks not run, limitations, branch/state, and uncommitted changes.
Give the best next engineering step and identify Git actions still requiring owner
approval. Suggest targeted staging/commit commands and a professional commit message
only when appropriate; never execute them without approval. Recommend a pull request
only when repository state and approved workflow justify it, with base/head, title,
description, and testing details. Never imply an unperformed check passed.
