# Test graphs

`make check` runs both `graphswitching` and the original
`gen_all_srgs_wqh_generic` program for every entry in `cases.txt`, then
requires their output to be byte-for-byte identical.

The fixtures use the following deterministic vertex orderings:

| Fixture | Parameters | Construction | Test purpose |
| --- | --- | --- | --- |
| `cycle4` | 4 vertices, `p=1` | Vertices in cyclic order | Small switching and output smoke test |
| `petersen` | 10 vertices, `p=2` | The 2-subsets of `{0,...,4}` in lexicographic order; disjoint subsets are adjacent | Small nontrivial strongly regular graph |
| `clebsch` | 16 vertices, `p=2` | Four-bit words in numeric order; words at Hamming distance 1 or 4 are adjacent | Medium switching and output test |
| `symplectic-sp6-2` | 63 vertices, `p=3` | Nonzero vectors of `GF(2)^6` in numeric order; vectors are adjacent when their standard symplectic product is 1 | Requested Sp(6,2) case and a tuple-enumeration-heavy test |
| `paley73` | 73 vertices, `p=3` | Residues modulo 73 in numeric order; distinct residues are adjacent when their difference is a nonzero square | Larger tuple-enumeration-heavy test |

The two largest cases emit only the `n=...` header for these part sizes.
Consequently, their running time primarily exercises candidate-tuple
enumeration instead of matrix output.
