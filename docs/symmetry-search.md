# Symmetry-aware switching search

This note describes the current `--sym` implementation. The search changes
which switching operations are emitted, but preserves the resulting graph
isomorphism classes. It is not a final graph-isomorphism filter: distinct
retained operations can still produce isomorphic graphs.

The mathematical provenance of the catalogue methods and a comparison with
the reference search methodology are documented separately in
`docs/methodology.md`.

## Common setup

`--sym` computes `Aut(G)` once with nauty. A Schreier representation is
kept for point stabilizers. When the group has at most 4,096 elements, its
elements are also cached so small-set orbit keys and setwise stabilizers can
be obtained by filtering without another nauty call.

All modes use the same validity tests and switching transformation. The
ordinary search traverses labelled candidates directly; `--sym` replaces
those loops with orbit-representative loops and adds orbit-duplicate checks.

## GM4

The GM search constructs unordered four-sets. Its tests are scheduled by
cost: a stabilizer is delayed when trying the remaining vertices is cheaper.

```text
GM-SYM(G):
    A := empty set
    SEARCH-GM(A, vertex orbits of Aut(G))

SEARCH-GM(A, orbits):
    if |A| = 3:
        candidates := bit mask of vertices whose adjacency to A
                      makes G[A union {v}] regular
    else:
        candidates := all vertices not in A

    for one v from each current candidate orbit:
        B := sorted(A union {v})

        if |B| = 4:
            reject unless B is a valid, nontrivial GM switching set
            reject unless B is the first valid set in its Aut(G)-orbit
            emit the switched graph
            continue

        if |B| = 1:
            reject a repeated singleton orbit
            next_orbits := point-stabilizer orbits of v
        else if |B| = 2 and Aut(G) is cached:
            reject a repeated set orbit
            next_orbits := orbits of the setwise stabilizer Aut(G)_{B}
        else:
            reject a repeated literal subset
            next_orbits := singleton orbits

        SEARCH-GM(B, next_orbits)
```

The three-set stabilizer is always postponed. For groups too large to cache,
the two-set stabilizer is postponed as well. Complete set orbits are marked
only after GM validity is known; an invalid set cannot hide a valid member of
the same orbit.

Without `--sym`, GM uses the same regular-completion bit mask but visits
every unordered four-set and performs no orbit checks.

## WQH

For part size `p`, WQH constructs `C1` and `C2` separately. It retains
the delayed regularity schedule: internal tests start only after more than
half a part is known.

```text
WQH-SYM(G, p):
    construct C1 through successive Aut(G) stabilizer orbits
    when |C1| = p:
        reject a repeated C1 orbit
        compute the setwise stabilizer Aut(G)_{C1}
        cache every vertex's adjacency signature on C1

    construct C2 through successive Aut(G)_{C1} stabilizer orbits
        use specialized candidate masks for the small supported p
        apply internal and cross-regularity tests only when informative
        after the second C2 vertex, reject if too many outside vertices
            are forced into C1 union C2 than remaining positions

    when |C2| = p:
        reject a repeated C2 orbit under Aut(G)_{C1}
        run the full outside-block validity test
        emit the switched graph
```

Without `--sym`, the same specialized candidate masks and validity
schedule are used, but every labelled partition is visited. The delayed
outside-block feasibility scan is also used in ordinary mode for `p=3,4`,
where it reduces work; it is disabled for `p=5`, where its extra scan costs
more than it saves.

## Fixed catalogue methods

The fixed catalogues are selected by `gm -p 3`, `gm -p 4`, `ah`,
`gm2`, `is3`, `is5`, and `fano`. Generic `wqh` retains its
specialized partition search.

```text
FIXED-SYM(G, Q, catalogue):
    derive every 0/1 block b for which Q^T b is also 0/1

    for each irreducible labelled switching subgraph Gamma:
        MATCH(0)

    MATCH(depth):
        if depth = |Gamma|:
            apply Q^T A Q
            emit the graph if the transformation changes it
            return

        compute point-stabilizer orbits after the selected graph vertices

        for one graph vertex v from each orbit:
            reject if v was already selected
            reject unless its adjacency to selected vertices
                agrees with the next labelled vertex of Gamma
            reject if too many outside vertices cannot extend
                to an admissible block
            select v
            MATCH(depth + 1)
            unselect v
```

Ordinary fixed-method enumeration uses the same template order,
pair-compatibility checks, block-feasibility pruning, and transformation,
but tries every eligible graph vertex instead of one per stabilizer orbit.
