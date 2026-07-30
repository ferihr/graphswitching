# Test graphs

`make check` runs the explicit WQH method in `graphswitching` and the
original `gen_all_srgs_wqh_generic` program for every entry in `cases.txt`,
then requires their output to be byte-for-byte identical. It also compares
the default GM method with `gen_all_srgs` on Sp(6,2), and checks automatic
order detection, `n=...` headers, file input/output, help, and version
reporting. `generate_algebraic_fixtures.py` independently reconstructs the
named symplectic, bilinear forms, and generalized-quadrangle graphs and
checks both their stored matrices and strongly regular parameters.

The fixtures use the following deterministic vertex orderings:

| Fixture | Parameters | Construction | Test purpose |
| --- | --- | --- | --- |
| `cycle4` | 4 vertices, `p=1` | Vertices in cyclic order | Small switching and output smoke test |
| `petersen` | 10 vertices, `p=2` | The 2-subsets of `{0,...,4}` in lexicographic order; disjoint subsets are adjacent | Small nontrivial strongly regular graph |
| `clebsch` | 16 vertices, `p=2` | Four-bit words in numeric order; words at Hamming distance 1 or 4 are adjacent | Medium switching and output test |
| `generalized-quadrangle-gq2-4` | `SRG(27,10,1,5)`, `p=5` | Nonzero singular vectors of the elliptic quadratic form on `GF(2)^6`, in numeric order; orthogonal points are adjacent | Small output-producing test for WQH `5,5,n-10` |
| `symplectic-sp6-2` | `SRG(63,30,13,15)`, `p=3` | Nonzero vectors of `GF(2)^6` in numeric order; distinct vectors are adjacent when their standard symplectic product is 0 | Sp(6,2) comparison and a tuple-enumeration-heavy test |
| `paley73` | 73 vertices, `p=3` | Residues modulo 73 in numeric order; distinct residues are adjacent when their difference is a nonzero square | Larger tuple-enumeration-heavy test |
| `bilinear-forms-bil2-2-3` | `SRG(81,32,13,12)`, `p=3` | The `2` by `2` matrices over `GF(3)` in lexicographic order; matrices are adjacent when their difference has rank 1 | Full benchmark for WQH `3,3,n-6` |
| `symplectic-sp4-4` | `SRG(85,20,3,5)`, `p=4` | Projective points of `GF(4)^4` in normalized lexicographic order; distinct points are adjacent when their standard symplectic product is 0 | Full benchmark for WQH `4,4,n-8` |

The Sp(6,2) and Paley(73) `p=3` cases emit only the `n=...` header.
Consequently, their running time primarily exercises candidate-tuple
enumeration instead of matrix output. The two other large fixtures exercise
the one-step switching cases reported for their respective SRG families.

Regenerate the algebraically defined fixtures after intentionally changing
their construction:

```sh
python3 tests/generate_algebraic_fixtures.py --write
```
