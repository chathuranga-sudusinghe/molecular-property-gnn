# Primary hybrid structure-aware split

## Decision and execution status

The owner accepted the primary methodology and allocation contract on 2026-10-20.
The protocol identifier is `esol-hybrid-structure-v1`.

Implementation and bounded synthetic checks are authorized. Real ESOL construction,
owner acceptance of generated artifacts, and procedural test sealing are separate
steps. This document and an implemented command do not authorize execution.

The implementation uses the local CPU. GPU/CUDA support belongs to a separately
approved model-training milestone and has no role in this structural split.

## Input and target isolation

The approved input is `data/raw/ESOL_delaney-processed.csv`, with SHA-256:

```text
8c06a76f0c6487d29ab0f903e6a7a7139f189ab3c1178f159c8be8964602f189
```

Expected populations are 1,128 records, including 317 empty-Murcko records and 811
nonempty-Murcko records. Raw bytes and source records remain immutable.

The CSV boundary verifies the checksum, unique required headers, field counts, and
CSV syntax. It projects each record immediately to `record_number`, `Compound ID`,
and `smiles`. Other fields are not interpreted or passed to structural functions.
Whole-file hashing necessarily includes all source bytes; that is not target analysis.

Measured targets, historical predictions, and existing descriptor columns must not
influence identities, fingerprints, groups, thresholds, allocation, or acceptance.
There is no target-distribution acceptance criterion.

Invalid/empty structures, unexpected disconnected molecules, population changes,
and checksum changes stop execution without automatic repair or exclusion.

## Molecular identity

Use RDKit 2024.09.6, pinned by the existing project environment:

1. Strip surrounding SMILES whitespace in memory.
2. Parse with sanitization, `parseName=False`, `allowCXSMILES=False`, and the
   remaining pinned RDKit defaults, including its normal hydrogen handling.
3. Clear atom maps on key copies.
4. Full identity: `MolToSmiles(canonical=True, isomericSmiles=True)`.
5. Conservative identity: remove stereochemistry on a copy, then use canonical
   isomeric SMILES; isotope distinctions remain.

Both identity matches are mandatory co-assignment constraints. Unknown stereochemistry
remains unknown. Do not desalt, neutralize, normalize tautomers, infer experimental
species, remove duplicates, or average measurements. Each original record retains
its own sample ID and record-count weight.

These are graph identities under the declared policy, not a claim that all chemically
equivalent protonation states, tautomers, or experimental species have been merged.

## Nonempty-Murcko route

Use the pinned RDKit implementation:

```python
scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
key = Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=False)
```

Route according to the scaffold atom count, not a hand-written ring heuristic.
Keep equal nonempty keys together. Preserve typed atoms and bonds; do not call
`MakeScaffoldGeneric`. RDKit's terminal/exocyclic decomposition behavior is part of
the definition. Scaffold serialization excludes stereo and isotope annotations;
molecular identity still distinguishes isotopes.

## Empty-Murcko route: Policy A

Create one node per conservative identity, retaining every associated raw record.

Use `GetMorganGenerator` and `GetSparseCountFingerprint`:

```text
radius=2
countSimulation=False
includeChirality=False
useBondTypes=True
onlyNonzeroInvariants=False
includeRingMembership=True
countBounds=None
fpSize=2048
atomInvariantsGenerator=None
bondInvariantsGenerator=None
includeRedundantEnvironments=False
```

`fpSize` is recorded for fidelity with the diagnostic; it does not fold this sparse
output into a 2,048-bit fingerprint. Sparse feature identifiers still have hashing
limitations. Counts retain multiplicity of local environments and are sensitive to
molecular size and repeated environments. These fingerprints are for structural
construction/audit only, not model features.

For nonnegative counts, define intersection `I = sum(min(x_k, y_k))` and union
`U = sum(max(x_k, y_k))`. Tanimoto is `I/U`. Reject undefined empty fingerprints.

Add an undirected edge iff `T >= 0.70`, evaluated exactly as `10*I >= 7*U`.
Take connected components, including isolated nodes.

Use **similarity-connected components**, not "chemical families." Transitive paths
can join distant endpoints: members need not be mutually similar. The guarantee is
that every acyclic pair crossing any partition boundary has `T < 0.70` under this
representation.

The existing diagnostic remains historical evidence:

| Threshold | Components | Largest record count | Minimum similarity in largest component |
| --- | ---: | ---: | ---: |
| 0.60 | 220 | 74 | 0.0580 |
| 0.70 | 255 | 18 | 0.2600 |
| 0.80 | 272 | 13 | 0.2157 |

At 0.70 the second-largest component has 16 records, minimum similarity 0.0988.
All three thresholds were count-feasible in the acyclic diagnostic; this did not
establish the result of the approved allocator. Preserve
`data/interim/esol_acyclic_similarity_diagnostic_v1.json` unchanged.

The chosen 0.70 threshold is an accepted project-specific operational boundary,
not a universal chemistry standard or a model-performance optimum.

## Combining constraints and allocation

Conservative nodes already co-assign duplicate records. Nonempty nodes group by
scaffold; empty nodes group by threshold components. Verify full and conservative
identity consistency with routes/scaffolds, and stop on unexpected conflicts.

Allocate routes independently, targeting train/validation/test fractions 4/5, 1/10,
1/10 by source-record count. Groups are indivisible.

For each route:

1. Sort by decreasing group record count, then complete stable key lexicographically.
2. Evaluate each candidate placement using the exact rational objective
   `sum((n_s - q_s*N)**2 / (q_s*N))`.
3. Exclude placements leaving fewer remaining groups than still-empty partitions.
4. Choose minimum cost; ties use train, validation, then test.
5. Make one pass; do not search seeds, thresholds, orders, or alternative candidates.

Scaffold keys are `murcko:` plus scaffold SMILES. Component keys are `acyclic:` plus
compact JSON of the sorted conservative identities. JSON uses UTF-8, sorted object
keys, no ASCII escaping, and separators `(",", ":")`.

Largest-first placement is for reproducibility and leakage control. It may favor
training placement for large groups and does not promise globally optimal counts.
It is not optimized for favorable model results.

Require both routes in every partition. Tolerance is inclusive +/-2 percentage
points around each target fraction, globally and within each route. Structural
integrity has priority. A tolerance failure stops after this allocation and leaves
an explicitly rejected, unsealed candidate for owner review. No automatic fallback.

## CLI: requires separate real-data execution approval

From the repository root:

```bash
.venv/bin/python -B -m src.data.hybrid_split \
  data/raw/ESOL_delaney-processed.csv \
  --expected-source-sha256 8c06a76f0c6487d29ab0f903e6a7a7139f189ab3c1178f159c8be8964602f189 \
  --expected-record-count 1128 \
  --expected-acyclic-count 317 \
  --output-prefix data/processed/esol_hybrid_structure_v1
```

The output parent must exist. There are no threshold, ratio, seed, order, sealing,
or overwrite flags. The source checksum is checked before construction and again
before publication. A different dataset revision requires a new owner-approved run
and artifact version.

Exit codes:

- `0`: validated candidate, pending owner artifact acceptance; not sealed.
- `1`: acceptance failure; no alternate allocation attempted.
- `2`: input, environment, configuration, or I/O failure.

## Artifacts and reproducibility

The command creates only these persistent outputs:

- `data/processed/esol_hybrid_structure_v1_manifest.csv`
- `data/processed/esol_hybrid_structure_v1_metadata.json`
- `data/processed/esol_hybrid_structure_v1_acceptance.json`

Manifest columns, in order:

```text
sample_id,source_record_number,compound_id,full_identity_id,
conservative_identity_id,structure_route,structural_group_id,split
```

The actual CSV header is one line. It uses UTF-8, LF line endings, and source-record
order. Sample IDs are `esol:<full-source-sha256>:<six-digit-record-number>`.
Compound IDs are preserved and need not be unique.

Key IDs are SHA-256 of compact JSON
`{"schema":1,"kind":<kind>,"key":<complete-key>}`, with the serialization rules above.
Kinds are `full_identity`, `conservative_identity`, and `structural_group`.
Metadata retains the complete key mappings and checks for digest collisions.

Metadata records configuration and its hash, source revision, group memberships,
Python/RDKit versions, OS/architecture, CPU execution, code/configuration/lockfile
hashes, read-only Git provenance where available, execution arguments, UTC times,
and measured construction/check duration excluding publication. An unborn repository
has a null commit; source hashes identify its local code.

The acceptance report contains check outcomes, counts/deviations, identity coverage,
group-size distributions, the five largest acyclic components' chaining diagnostics,
and per-held-out-identity nearest-training neighbors. Quantiles use linear
interpolation at `(n-1)*q`. Fixed bins are [0,.2), [.2,.4), [.4,.6), [.6,.7),
[.7,.8), [.8,1]. Nearest-neighbor ties use conservative identity lexical order.

Nearest-neighbor comparisons are unique conservative identities, not record weights.
They report all-training, same-route, and cross-route comparisons separately.
Only empty-to-empty threshold separation is enforced. Other high similarities are
descriptive findings, not an undeclared objective for split searching.

Metadata hashes the manifest and acceptance report; the report hashes the manifest.
All three hashes are printed. No file contains its own digest.

Outputs contain no labels or historical predictions, and no pairwise matrix is saved.
Temporary sibling files are removed on handled failures. Exclusive publication
refuses existing outputs, including racing writers. Acceptance is published last.
A usable bundle requires all three files, matching hashes, and passing checks;
a lone failure report or interrupted bundle is never usable. Separate files are
not a filesystem-wide atomic transaction.

Structural failures produce only a rejected acceptance report when construction
reaches validation. Input failures publish no artifacts. Tolerance-only failures
preserve all three artifacts, explicitly rejected and unsealed.

## Acceptance and procedural test seal

Required checks cover:

- Every approved record exactly once, unique sample IDs, correct manifest fields.
- Valid, mutually exclusive assignments and nonempty partitions.
- Full and conservative identity integrity.
- Nonempty scaffold integrity.
- Independently re-enumerated qualifying acyclic edges across all partition pairs.
- Complete group membership and group integrity.
- Both routes in every partition and global/per-route tolerance.
- Source checksum unchanged.
- Structural coverage and group-size/chaining diagnostics.
- Label-blind nearest-training-neighbor reporting.

Synthetic restricted-field and target-perturbation tests provide mechanical isolation
evidence. They cannot prove that a human has never considered public labels.

Passing construction does not freeze or seal artifacts. The owner must accept
the exact hashes and review any blockers before freezing one manifest for every
traditional ML and GNN comparison.

During later development, fit learned transformations and models on training only.
Use validation for model/configuration selection and approved calibration decisions.
Do not refit on train+validation without a later owner-approved protocol.
Development loaders/evaluation interfaces must enforce test access restrictions;
those model/evaluation interfaces are outside this split implementation.

Before final test release, freeze selected models, evaluation settings, seeds,
metrics, and reporting rules. Log release/evaluation; do not use test results to
reopen development and still call that test untouched. Public ESOL labels remain
accessible: the seal is procedural, not cryptographic.

For later evidence-producing experiments, additionally record the manifest hash,
representations, fitted preprocessing, model/hyperparameter configuration, tuning
budgets, training seeds, device, durations, predictions/checkpoints, metrics, and
failed runs. No models or experiment infrastructure are implemented here.

## Tests

Only bounded synthetic checks are authorized for this implementation milestone:

```bash
.venv/bin/python -B -m pytest -q -p no:cacheprovider \
  tests/test_hybrid_split.py \
  tests/test_acyclic_diagnostic.py \
  tests/test_structural_audit.py

.venv/bin/python -B -m ruff check \
  src/data/hybrid_split.py tests/test_hybrid_split.py

.venv/bin/python -B -m ruff format --check \
  src/data/hybrid_split.py tests/test_hybrid_split.py
```

These verify mechanisms with tiny synthetic inputs, not real ESOL allocation,
tolerance success, scientific validity, or actual test sealing.

## Literature basis and limits

The motivation is literature-supported; the exact hybrid, count configuration,
allocation, tolerance, and threshold are project-specific.

- Wu et al., *MoleculeNet: a benchmark for molecular machine learning*,
  Chemical Science 2018, 9, 513-530, DOI: 10.1039/C7SC02664A.
  Local PDF: `../research-materials/2018_moleculenet-benchmark.pdf`.
  Section 3.2, PDF p.5 supports multiple split strategies and 80/10/10 roles.
  Table 1, PDF p.4 documents historical ESOL random splitting and RMSE.
- Kovacs et al., *Quantitative interpretation explains machine learning models for
  chemical reaction prediction and uncovers bias*, Nature Communications 2021,
  12, 1695, DOI: 10.1038/s41467-021-21895-w.
  Local PDF: `../research-materials/2021_kovacs_ml-reaction-prediction-bias.pdf`.
  PDF pp.5-7 demonstrate reaction-prediction bias. Methods p.8 uses radius-3,
  1,024-bit fingerprints; the p.6 split is not our three-way count-component method.
- Kretschmer et al., *Coverage bias in small molecule machine learning*,
  Nature Communications 2025, 16, 554, DOI: 10.1038/s41467-024-55462-w.
  Local PDF: `../research-materials/2025_coverage-bias-small-molecule-ml.pdf`.
  PDF pp.2 and 5 limit generalization claims; pp.11-13 discuss fingerprint
  limitations and binary fingerprint examples, not direct testing of our counts.
- Li et al., *A Systematic Survey and Benchmark of Deep Learning for Molecular
  Property Prediction in the Foundation Model Era*, Journal of Chemical Theory
  and Computation 2026, 22, 4866-4887, DOI: 10.1021/acs.jctc.5c02081.
  Local PDF: `../research-materials/2026_molecular-property-prediction-survey.pdf`.
  Sections 6.2-6.3, PDF pp.11-13 discuss split realism and reproducibility;
  Sections 8.2-8.4, pp.15-16 discuss evaluation limitations and uncertainty.

The primary is one fixed structural evaluation, not evidence across split variation.
Training-seed variation is not split uncertainty. Small size, the large benzene
group, retained duplicates, experimental metadata limits, unnormalized chemistry,
achiral grouping, chaining, representation/size bias, and cross-route familiarity
limit interpretation. No deployment population or prospective realism is established.

Permitted after acceptance: precise claims about identity integrity, scaffold
integrity, acyclic cross-partition Tanimoto bounds, shared frozen assignments, and
the documented label-blind construction process.

Do not claim universal chemical generalization, completely novel chemistry, mutually
similar component members, a standard 0.70 boundary, perfect similarity, global
hybrid Tanimoto separation, inaccessible test labels, exact MoleculeNet reproduction,
guaranteed exact proportions/optimal allocation, or statistical superiority from
one fixed split.

## Secondary protocols

A duplicate-aware random baseline and a conventional MoleculeNet-style scaffold
benchmark are approved directions only. Each requires a separately pinned
implementation specification and execution approval. Neither is implemented here.

Freeze their plans before development; run their development/evaluation after the
primary evaluation is completed and locked. Disclose cross-protocol label exposure.
Secondary results must not reopen primary decisions.
