# Delaney ESOL dataset and initial audit

## Purpose and target

The local source is `data/raw/ESOL_delaney-processed.csv`. Its intended unit of
analysis is one compound measurement record; this audit does not assume that a
SMILES string uniquely identifies a standardized chemical entity.

- Structure column: **`smiles`**, containing Simplified Molecular Input Line Entry
  System (SMILES) representations.
- Regression target: **`measured log solubility in mols per litre`**. It represents
  aqueous `logS = log10(S / (mol/L))`, where `S` is measured solubility. A one-unit
  increase corresponds to tenfold higher molar solubility. Keep this target in its
  supplied logarithmic units.
- **`ESOL predicted log solubility in mols per litre` is a historical prediction,
  not a measured label.** The existing descriptor and prediction columns are
  profiled as source metadata only. They are not selected as model inputs, and no
  descriptors are calculated during this stage.

## Provenance

Original research: John S. Delaney (2004), *ESOL: Estimating Aqueous Solubility
Directly from Molecular Structure*, Journal of Chemical Information and Computer
Sciences, 44(3), 1000–1005. [DOI: 10.1021/ci034243x](https://doi.org/10.1021/ci034243x).
The publisher provides [supporting data](https://acs.figshare.com/articles/dataset/ESOL_Estimating_Aqueous_Solubility_Directly_from_Molecular_Structure/7944677).

DeepChem distributes a processed 1,128-record benchmark; its
[Delaney loader](https://github.com/deepchem/deepchem/blob/master/deepchem/molnet/load_function/delaney_datasets.py)
identifies `smiles` and the measured target and points to
[delaney-processed.csv](https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv).
These references establish the dataset family and target mapping. The local file's
exact download URL, download date, upstream revision, and processing history have
not been supplied or independently verified. Its Windows download-zone sidecar
does not record a source URL. Do not treat the filename as proof of byte identity
with a particular upstream release.

The local CSV SHA-256 is:

```text
8c06a76f0c6487d29ab0f903e6a7a7139f189ab3c1178f159c8be8964602f189
```

The audit hashes the same byte snapshot it parses and records library versions.
No source data is redistributed by this change. Dataset redistribution terms
must be checked against the actual source before any later publication.

## Validation contract

Implementation: `src/data/audit.py`. Execution instructions are in the README.

1. Read the CSV in memory only. No row removal, data export, canonicalization,
   chemical standardization, feature engineering, splitting, or training occurs.
2. Require the ten exact headers below in their supplied order. Reject missing,
   extra, repeated, or reordered headers, malformed quoting, and ragged or blank
   CSV records with an explicit input error. Reject NUL bytes to prevent
   parser-dependent truncation. A header-only file is reported as
   empty and fails content validation.
3. Report pandas-inferred data types; read compound identifiers and SMILES as
   text. Missing means null or empty/whitespace-only text. Literal strings such
   as `NA` are preserved rather than interpreted as missing-value sentinels.
4. Report missingness for every column. Report duplicate full rows, exact
   duplicate populated SMILES, and duplicate populated SMILES after stripping
   surrounding whitespace in an in-memory copy. `extra_rows` counts occurrences
   after the first; `affected_rows` counts all members of duplicate groups.
   Each duplicate-SMILES group lists its one-based data-record numbers.
   Missing SMILES are excluded from duplicate-key checks and reported separately.
5. Parse each trimmed SMILES with RDKit's default sanitization. Require at least
   one atom. Disable optional name-suffix and CXSMILES parsing so a valid SMILES
   prefix cannot hide trailing text. Empty/non-string values and parsing or
   sanitization failures are invalid. Multi-fragment structures and specified
   stereochemistry are accepted if RDKit accepts them. RDKit diagnostics remain
   visible on stderr; structured failures appear in the JSON report.
6. Require finite numeric measured targets. Report invalid/missing/infinite
   targets individually. Compute statistics and a ten-bin histogram only on
   finite targets, explicitly reporting the contributing count. No physical
   target bounds or outlier-removal rule is imposed.
7. Content validation passes only for nonempty data with no missing cells,
   invalid SMILES, or invalid targets. Whitespace and duplicate findings do not
   fail validation or authorize cleaning. A passing report is not model-readiness
   certification.

Invalid values are retained in the report with their compound identifier and
one-based data-record number (the header is excluded). Record numbers are based
on CSV record order, not DataFrame index labels or physical file lines; quoted
multiline fields can span multiple lines. The report contains no cleaned dataset.
Schema/structural failures stop the audit before content profiling.

## Initial snapshot (2026-10-20)

Complete local evidence is retained in `docs/esol_audit.json` and is not
included in Git. Observed schema and inferred types are listed below in original
order; `object` denotes text in this snapshot.

| Exact column name | pandas dtype | Missing cells |
| --- | --- | --- |
| `Compound ID` | object | 0 |
| `ESOL predicted log solubility in mols per litre` | float64 | 0 |
| `Minimum Degree` | int64 | 0 |
| `Molecular Weight` | float64 | 0 |
| `Number of H-Bond Donors` | int64 | 0 |
| `Number of Rings` | int64 | 0 |
| `Number of Rotatable Bonds` | int64 | 0 |
| `Polar Surface Area` | float64 | 0 |
| `measured log solubility in mols per litre` | float64 | 0 |
| `smiles` | object | 0 |

There are **1,128 rows and 10 columns**. There are no exact duplicate rows or
exact duplicate SMILES. **212 SMILES (18.79%) have surrounding whitespace**.
This can affect text-based comparisons; it is reported without rewriting data.

All **1,128 SMILES pass RDKit 2024.09.6 parsing and sanitization** after in-memory
whitespace trimming: **zero invalid molecules**. All 1,128 measured targets are
finite numeric values: **zero invalid targets**. The dataset passes the defined
content checks, subject to the duplicate and scientific limitations below.

Trimming reveals **five duplicate SMILES groups: five extra occurrences across
ten affected records (0.89% of all records)**. Three groups contain differing
measured targets. All ten records are retained:

| Data-record numbers | Compound identifiers | Measured logS values |
| --- | --- | --- |
| 214, 977 | eucalyptol; 1,8-Cineole | -1.64; -1.74 |
| 234, 277 | mannitol; Sorbitol | 0.06; 1.09 |
| 451, 1020 | DDD; P,P'-DDD | -7.20; -7.20 |
| 681, 1070 | Guaiacol; o-Methoxyphenol | -1.96; -1.96 |
| 704, 823 | Metronidazole; Metranidazole | -1.22; -1.26 |

This is a confirmed representation/measurement ambiguity, not evidence that
either measurement is incorrect. In particular, the supplied mannitol/sorbitol
SMILES omits stereochemistry; this audit does not infer it from compound names.
Before later evaluation, investigate these groups and define an explicit identity
and measurement-resolution policy to avoid overlapping representations across splits.

| Measured logS statistic | Value |
| --- | ---: |
| Count | 1,128 |
| Mean | -3.050102 |
| Sample standard deviation (`ddof=1`) | 2.096441 |
| Minimum | -11.600000 |
| 25th percentile | -4.317500 |
| Median | -2.860000 |
| 75th percentile | -1.600000 |
| Maximum | 1.580000 |
| Sample skewness | -0.492370 |

The middle 50% spans -4.3175 to -1.6, and the overall range spans 13.18 log units.
The mean is below the median, consistent with a longer low-solubility tail.
The observed negative skewness supports this descriptive summary.
The JSON histogram records all bin edges and counts rather than inferring a
parametric distribution or imposing arbitrary solubility classes.

## Current limitations and later decisions

- RDKit acceptance checks syntax and sanitization, not experimental correctness,
  chemical identity, chemical stability, or label quality.
- Exact and whitespace-trimmed text duplicates do not detect all equivalent
  molecular graphs, alternate SMILES encodings, salts, tautomers, or stereochemical
  equivalence. A later, explicitly approved standardization policy must investigate
  chemical duplicates and conflicting measurements before evaluation splits.
- The CSV has no per-record pH, temperature, measurement method, measurement
  uncertainty, or primary measurement citation. No common experimental conditions
  or intrinsic-solubility interpretation are assumed. Temporal/freshness analysis
  cannot be performed because there are no measurement timestamps.
- Auxiliary descriptor values and historical predictions are not recomputed or
  scientifically validated here. Their calculation methods and possible overlap
  with historical fitting data require investigation before use in comparisons.
- A small benchmark cannot establish generalization outside its chemical domain.
  Distribution summaries are descriptive; they are not evidence of calibrated
  uncertainty or future predictive performance.
- Tests use synthetic records so they work without the ignored raw dataset. They
  verify validation behavior, failure reporting, statistics, CLI exit codes, and
  immutability. They do not verify experimental measurements, upstream lineage,
  molecular standardization, dataset splits, or model behavior.
- No cleaning, duplicate resolution, train/validation/test splits, descriptor
  generation, graph construction, or model training is part of this stage.
