# Methodology and relationship to Simoens--Van Overberghe

## Provenance

The fixed switching matrices, irreducible switching-subgraph catalogues, and
paper benchmark examples are based on Robin Simoens and Steven Van
Overberghe, "An algorithm to find cospectral mates", and their SageMath
reference implementation:

https://github.com/robinsimoens/implementing-switching-methods

The matrices and catalogues were adapted with the authors' permission and
re-expressed as rational numerator matrices and compact edge masks.

## Common mathematical core

For the fixed methods, both implementations use the same mathematical data
and tests:

- A method is represented by a rational orthogonal switching matrix `Q` and
  a catalogue of irreducible switching subgraphs `Gamma`.
- The search embeds a catalogue graph `Gamma` as an induced subgraph of the
  input graph.
- Every vertex outside the embedding must have an admissible 0/1 incidence
  block `b`: the transformed block `Q^T b` must again be a 0/1 block.
- A valid embedding is switched by the conjugation `Q^T A Q`. Results for
  which the adjacency matrix is unchanged are not emitted.
- Canonical labelling of the resulting graphs can be applied afterwards to
  retain one representative of each graph isomorphism class.

Thus the catalogue data, admissibility conditions, and switching operation
are common. The implementations differ in how they enumerate embeddings and
remove equivalent search branches.

## Search organization

The Simoens--Van Overberghe implementation uses a general
canonical-construction path. A partial switching set is represented as a
matching between the input graph and a switching template. Canonical deletion
in a colored auxiliary graph determines whether the partial construction is
retained. This gives one systematic search framework for all switching
methods and accounts for symmetries of both the input graph and the template.

This program uses method-specific search instead:

- Ordinary GM and WQH enumeration retains family-specific candidate masks,
  delayed regularity tests, and delayed outside-block feasibility tests.
- Ordinary fixed-method enumeration traverses each catalogue template in a
  fixed labelled order, using pair-compatibility and block-feasibility tests
  before applying the switching matrix.
- With `--sym`, nauty computes `Aut(G)` once. The recursion then extends a
  partial switching set through point- or set-stabilizer orbits rather than
  canonically labelling a new auxiliary graph at every search node.
- The command-line interface groups catalogue instances into families; the
  family is selected by `--method` and its parameter by `--part-size`.

## Consequences

The family-specific recursion has lower per-node overhead on the supplied
benchmarks. Its duplicate guarantee is weaker: fixed template automorphisms
can leave multiple switching operations that the general matching
construction would identify. The `--sym` implementation guarantees that the
set of resulting graph isomorphism classes is preserved, not that exactly
one operation or output graph is emitted per class.

For that reason, `--sym` is an orbit-reduced switching search rather than a
final graph-isomorphism filter. Use a canonical-labelling pipeline such as
`amtog | labelg | sort -u` when exactly one output per isomorphism class is
required. Implementation pseudocode is in `docs/symmetry-search.md`.
