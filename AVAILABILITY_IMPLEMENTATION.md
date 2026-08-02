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
name is display metadata. Exact gene symbols are deduplicated within each group,
and multiple reactions in the same direction are summed into one HMDB-level
capacity.

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

The formal sender score is:

```text
sender_score(m,t) = P_score(m,t)^2 / (P_score(m,t) + C_score(m,t))
```

The result is set to zero when the denominator is zero. `E_score` remains an
export-support intermediate and does not enter the formal sender score.

Intermediates include raw `P`, `C`, and `E`, their positive-reference scores
and references, pseudobulk expression, expression fractions, cell counts, cell
fractions, abundance weights, and QC annotations. Obsolete signed contrasts and
duplicate score aliases are not returned.
