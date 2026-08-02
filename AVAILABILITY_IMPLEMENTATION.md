# Sender Score Implementation

Every observed cell type participates in pseudobulk construction, cell-fraction
adjustment, P/C/E references, scores, events, and permutations. `min_cells` is
used only to annotate cell-count QC.

Reaction activity is adjusted before P/C/E construction:

```text
reaction_activity(r,t)
    = geometric_mean({expression(g,t) + 1 | g belongs to reaction r}) - 1

adjusted_reaction_activity(t)
    = reaction_activity(t) * cell_fraction(t) ** sender_abundance_exponent
```

Genes within a reaction contribute equally. The enzyme prior has no numerical
`weight`; legacy/custom enzyme `weight` columns are ignored and removed during
prior validation.
The `reaction` field is required and must be non-empty; rows without a known
reaction identity are rejected instead of being merged into an artificial
multi-gene complex.
Reaction grouping uses `canonical_hmdb_id + reaction + direction`; metabolite
name is display metadata. Exact gene symbols define each reaction's gene set.
Within the same HMDB ID and direction, complete gene sets are compared
independent of gene order. Identical sets contribute once, and strict subsets
are omitted in favor of their supersets. Incomparable sets, including overlaps
where neither contains the other, retain their full gene sets and all
contribute. These inclusion-maximal reaction gene sets are summed into one
HMDB-level capacity, and reaction-count metadata counts the contributing sets.

For each metabolite and each direction `X` in `P`, `C`, and `E`, the reference
is calculated from strictly positive adjusted capacities:

```text
X_ref(m) = mean({X(m,t) | X(m,t) > 0})              default
X_ref(m) = median({X(m,t) | X(m,t) > 0})            optional
X_score(m,t) = 0                                    if X(m,t) = 0
X_score(m,t) = X(m,t) / (X(m,t) + X_ref(m))         otherwise
```

Because both `X(m,t)` and `X_ref(m)` are strictly positive in the division
branch, their sum is positive and no epsilon parameter is required.

The base sender score and bounded exporter modulation are:

```text
base_sender_score(m,t) = P_score(m,t)^2 / (P_score(m,t) + C_score(m,t))
E_factor(m,t) = (1 - export_weight) + export_weight * E_effective(m,t)
sender_score(m,t) = base_sender_score(m,t) * E_factor(m,t)
```

The base result is zero when the denominator is zero. `E_effective` equals
`E_score` for supported exporter evidence, 0.5 for missing or gene-unavailable
exporter priors, and 0 for measured all-zero exporter capacity. The default
`export_weight=0.2` bounds `E_factor` to `[0.8, 1.0]`; a zero weight recovers
the base result exactly.

Intermediates include raw `P`, `C`, and `E`, their positive-reference scores
and references, `base_availability`, `E_effective`, `E_factor`, pseudobulk
expression, expression fractions, cell counts, cell fractions, abundance
weights, and QC annotations. Obsolete signed contrasts and duplicate score
aliases are not returned.
