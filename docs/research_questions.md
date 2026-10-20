# Research questions and hypotheses

## Status, research title, and objective

**ACCEPTED research framing:** recorded under [Issue #8](https://github.com/chathuranga-sudusinghe/molecular-property-gnn/issues/8).
Acceptance of the framing is not evidence for the hypotheses, acceptance of an
implementation, or authorization to execute experiments.

**Research title:** Evaluating Graph Neural Networks for Structure-Aware Aqueous
Solubility Prediction Under Molecular Distribution Shift

**Objective:** Compare Graph Neural Networks (GNNs) with traditional
molecular-descriptor/fingerprint machine-learning models for Delaney ESOL aqueous
solubility prediction, using the same frozen primary manifest and a fair
structure-aware evaluation protocol. The shared manifest is a requirement for later
comparisons; its generation, artifact acceptance, freezing, and procedural test
sealing remain separate owner-approved steps.

See [the dataset contract](dataset.md) for the measured target definition,
provenance, validation rules, and dataset limitations.

## Research questions

### Main research question

**RQ1:** Do Graph Neural Networks provide superior aqueous solubility prediction
compared with traditional molecular-descriptor machine-learning models when
evaluated under structure-aware molecular generalization constraints?

### Supporting research questions

**RQ2:** How does predictive performance change between interpolation-oriented
random evaluation and the approved structure-aware evaluation protocol?

**RQ3:** How does prediction error change as molecular structural similarity to the
training set decreases?

**RQ4:** Which molecular structures or substructures contribute most strongly to
Graph Neural Network aqueous-solubility predictions?

**RQ5:** How reliable are model predictions for molecules that are structurally
distant from the training data?

## Hypotheses

**H0:** There is no meaningful difference in predictive performance between Graph
Neural Networks and traditional molecular-descriptor machine-learning models under
the approved structure-aware evaluation protocol.

**H1:** Graph Neural Networks achieve better predictive performance than traditional
molecular-descriptor machine-learning models under the approved structure-aware
evaluation protocol.

**H2:** Predictive error increases as molecular structural similarity to the training
set decreases.

These are accepted verbal research hypotheses, not observed findings or completed
statistical-test definitions. "Meaningful difference" has no numerical definition
yet and must not be equated with statistical significance. H0 and H1 do not yet
provide an exhaustive formal decision rule, including when GNNs perform worse.
Formal null/directional hypotheses and interpretation of all outcomes require a
later owner-approved statistical design; the accepted wording above is preserved.

## Scope and claim boundaries

The initial study concerns aqueous solubility in ESOL. Model configurations,
representations, tuning budgets, and statistical procedures require later pinned
specifications and owner approval. Neither GNN superiority nor any performance
difference is assumed.

The approved structural constraints define a limited evaluation of structural
generalization. The research title does not establish a prospective deployment
population, universal chemical novelty, or universal out-of-distribution
generalization. A single fixed manifest does not quantify variation across splits.
Variation across training seeds is distinct from uncertainty due to split choice.

Structural similarity is representation-dependent. The proposed error-versus-
similarity analysis investigates an association; it does not establish causation.
Likewise, model attributions must not be presented as proof of chemical causation.

This documentation milestone implements no features, models, split assignments,
statistical tests, explainability methods, uncertainty methods, or experiments.
It does not authorize dataset/target access or real-data/model execution.

## Relationship to the accepted primary protocol

The [primary hybrid structure-aware split contract](primary_split.md) remains the
authoritative specification for molecular identities, grouping, allocation,
artifact acceptance, and procedural test sealing. This document does not change
or duplicate that specification.

The primary protocol combines typed Bemis-Murcko scaffold grouping for nonempty
scaffolds with similarity-connected components for empty scaffolds. Describe it
as a hybrid structure-aware split, not a pure scaffold split or grouping into
"chemical families." Its guarantees and limitations remain those in the contract.

All traditional ML and GNN primary comparisons must use the same owner-accepted,
frozen primary manifest and approved evaluation populations. Fit learned
preprocessing and models on training data only. Use validation for model and
configuration selection and approved calibration procedures. Do not refit on
train plus validation without a later owner-approved protocol change.

The final test partition is procedurally sealed; public ESOL labels are not
cryptographically inaccessible. Test evaluation requires separate explicit owner
authorization after evaluation settings are frozen. Test results must not reopen
model selection or methodology development.

## Relationship to future milestones

| Question or hypothesis | Future work and boundary |
| --- | --- |
| RQ1, H0, H1 | Compare approved traditional models and GNNs on the shared primary manifest with documented tuning budgets and statistical design. |
| RQ2 | Use the secondary duplicate-aware random baseline after its implementation specification and execution receive separate approval. |
| RQ3, H2 | Analyze prediction error against nearest-training-neighbor structural similarity under a separately approved analysis specification. |
| RQ4 | Investigate explainability in a later milestone; interpretation remains model attribution, not chemical causation. |
| RQ5 | Investigate uncertainty and applicability-domain behavior in a later milestone; similarity alone is not calibrated predictive uncertainty. |

Secondary protocols remain separate from the primary evaluation and must follow
the sequencing and label-exposure rules in the primary contract. Secondary results
must not influence or reopen primary decisions. No secondary protocol is
implemented or executed by this document.

## Open statistical-design decisions

The following are **PROPOSED future design work**, not accepted specifications:

- Statistical test.
- Repeated-seed design.
- Confidence-interval procedure.
- Smallest Effect Size of Interest (SESOI).
- Superiority margin.
- Multiplicity handling.

No p-value cutoff, effect-size threshold, superiority margin, power assumption,
seed count, or statistical conclusion is established here. The formal hypotheses,
including treatment of outcomes unfavorable to GNNs, remain to be specified.

Accepted research framing, future methodology decisions, and future experimental
results must remain distinguishable. Later methodology and execution approvals
do not arise from recording these questions or hypotheses.
