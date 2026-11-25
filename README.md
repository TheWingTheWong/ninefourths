# ninefourths

On November 24, 2026 OpenAI posted a preprint claiming that the matrix multiplication exponent ω is exactly 9/4. The upper bound comes from October, and the new lower bound has not yet been checked outside OpenAI. Rafael Pass then published an explicit, constructive algorithm for the upper bound ([arXiv 2610.28627](https://arxiv.org/abs/2610.28627)). This repo builds that algorithm, runs it, and works out what it costs.

**The best exact identity found here whose exponent is below Strassen's 2.807 multiplies matrices whose side has 1,078,626 digits (V = 10^1,078,625). At that size it uses 10^3,028,063 multiplications, 22 orders of magnitude fewer than Strassen's 10^3,028,085, but nothing smaller can be fed to it without padding.** That is the best point of a grid search over the free parameters of Pass's lemmas (c ≤ 60, d ≤ 2500, k up to 64× its first feasible value). Pass's own choices of c and d (the smallest the proof allows) give an identity of side 10^(1.6×10^11). Almost all of the size comes from removing the formal parameter ε at the end (Bini's conversion). Before that step, the identity already has exponent 2.71 at side 10^17,977, but its ε-degree is 10^103,314, so Bini's conversion needs 60 tensor powers.

Along the way the build turned up one fact the paper doesn't state: **the amplified packing of Lemma 3.2 is exact.** The ε cancels identically on the balanced targets the lemma keeps, so every ε in the final algorithm comes from the interpolation in Lemma 3.4. It only changes the final size by about 0.8%, because the ε that remains still compounds every round.

![verify output](assets/verify.svg)

Every construction in the proof is a working program here. The primal (c+2)-for-c gadget, role exchanges, tensoring, amplification, the same-left-to-general conversion, interpolation, the round, composition, and Bini's conversion are each built as an explicit bilinear algorithm. Each is checked exactly at random points of F_p, with ε carried symbolically. The whole chain then multiplies real integer matrices and matches schoolbook. Two random 2×2 matrices take 373,248 multiplications, against 8 for schoolbook.

![cost chart](assets/cost.svg)

Every row is produced by [`bench/run.sh`](bench/run.sh) (the `cost` binary) and stored in [`bench/cost_table.jsonl`](bench/cost_table.jsonl). CI recomputes the 2.807 rows and checks them against the file.
- **Pass:** Pass's (c, d), meaning the smallest c whose limit exponent is below the target and the smallest d satisfying his condition (10). The paper leaves k and h as "sufficiently large", so they are chosen the way the search chooses them.
- **Searched:** the best (c, d, k, h) on the grid. The grid is c ≤ 60 and d ≤ 2500 (d ≤ 8000 at 2.4), with k sampled from its first feasible value k0 up to 64·k0.
- **Columns:** h is the number of tensor powers Bini's conversion needs. "Strassen at side V" is 7^⌈log2 V⌉.

| target exponent | parameters | c | d | k | h | mod-ε side v (exponent) | ε-degree | exact side V | multiplications R | Strassen at side V |
|---|---|---|---|---|---|---|---|---|---|---|
| log2 7 = 2.807 | Pass | 3 | 133 | 53,867 | 25,801 | 10^6,128,780 (2.807) | 10^(1.99×10^7) | **10^(1.58×10^11)** | 10^(4.44×10^11) | 10^(4.44×10^11) |
| log2 7 = 2.807 | searched | 5 | 209 | 147 | 60 | 10^17,977 (2.712) | 10^103,314 | **10^1,078,625** | 10^3,028,063 | 10^3,028,085 |
| 2.75 | Pass | 3 | 292 | 261,393 | 121,846 | 10^(6.64×10^7) (2.750) | 10^(2.12×10^8) | **10^(8.10×10^12)** | 10^(2.23×10^13) | 10^(2.27×10^13) |
| 2.75 | searched | 6 | 265 | 194 | 76 | 10^25,704 (2.655) | 10^184,956 | **10^1,953,502** | 10^5,372,090 | 10^5,484,174 |
| 2.7 | Pass | 3 | 2,411 | 7,324,047 | 3.34×10^6 | 10^(1.56×10^10) (2.700) | 10^(4.93×10^10) | **10^(5.21×10^16)** | 10^(1.41×10^17) | 10^(1.46×10^17) |
| 2.7 | searched | 6 | 318 | 237 | 93 | 10^37,621 (2.622) | 10^271,289 | **10^3,498,796** | 10^9,446,684 | 10^9,822,362 |
| 2.6 | Pass | 4 | 654 | 1,010,724 | 407,410 | 10^(4.60×10^8) (2.600) | 10^(2.05×10^9) | **10^(1.87×10^14)** | 10^(4.87×10^14) | 10^(5.26×10^14) |
| 2.6 | searched | 7 | 545 | 428 | 152 | 10^101,568 (2.542) | 10^888,412 | **10^(1.54×10^7)** | 10^(4.01×10^7) | 10^(4.33×10^7) |
| 2.5 | Pass | 5 | 6,879 | 123,477,188 | 4.51×10^7 | 10^(4.96×10^11) (2.500) | 10^(2.87×10^12) | **10^(2.24×10^19)** | 10^(5.59×10^19) | 10^(6.28×10^19) |
| 2.5 | searched | 10 | 1,084 | 934 | 324 | 10^319,673 (2.458) | 10^4,365,716 | **10^(1.04×10^8)** | 10^(2.59×10^8) | 10^(2.91×10^8) |
| 2.45 | Pass | 7 | 1,307 | 10,776,252 | 3.78×10^6 | 10^(6.08×10^9) (2.450) | 10^(5.37×10^10) | **10^(2.30×10^16)** | 10^(5.63×10^16) | 10^(6.45×10^16) |
| 2.45 | searched | 13 | 1,743 | 1,617 | 536 | 10^698,348 (2.415) | 10^(1.32×10^7) | **10^(3.74×10^8)** | 10^(9.17×10^8) | 10^(1.05×10^9) |
| 2.4 | Pass | 8 | 24,771 | 2,341,145,379 | 7.88×10^8 | 10^(2.25×10^13) (2.400) | 10^(2.32×10^14) | **10^(1.77×10^22)** | 10^(4.25×10^22) | 10^(4.97×10^22) |
| 2.4 | searched | 16 | 3,207 | 3,043 | 984 | 10^1,991,526 (2.375) | 10^(4.90×10^7) | **10^(1.96×10^9)** | 10^(4.70×10^9) | 10^(5.50×10^9) |

The "Pass" line jumps around because the smallest admissible c changes with the target. When that c's limit exponent sits just under the target (c = 3 has limit 2.69, so a target of 2.7 is barely above it), condition (10) forces a huge d and then a huge k. The search escapes this by taking a larger c (6 at a target of 2.7) and moving k well past feasibility. That puts the mod-ε identity clearly below the target (2.62 against 2.7) before ε is removed.

## Quickstart

```sh
cargo run --release --bin verify     # 451 exact checks, one line per lemma (~6 min, several GB of memory; --quick skips the 2x2 run)
cargo run --release --bin cost       # faithful + searched cost for target exponent log2 7
TAUS=strassen,2.6 CMAX=60 DMAX=2500 cargo run --release --bin cost   # more targets, wider search
./bench/run.sh                       # regenerate the whole table and chart (~30 min)
cargo test --release -- --include-ignored
```

Rust stable, no unsafe. The only dependency is `num-bigint`, used for the exact integers in the cost model.

## What is implemented

| Paper | Code | Checked |
|---|---|---|
| Def. 2.1 selected products, Full / Mid supports (2), (9) | `src/shape.rs` | Full^{12}_{a,b} = Full_{b,a}, Mid = dual of Full |
| Def. 2.2 packings (SL and general), modulo ε | `src/packing.rs` | exact verifier with ε symbolic |
| Thm 3.1 primal (c+2)-for-c packing | `gadgets::primal` | 304 parameter sets valid; all 64 sets with ct < a+2b−1 rejected |
| OpenAI Lemma 4.2 three-for-one degeneration | `primal(a, b, 1, a+2b−1, ..)` | 16 sets |
| Lemma 2.7 tensoring, 2.8 role exchanges, 2.9 composition (with the ε → ε^N substitution), 2.10 renaming | `Packing::{tensor, ex12, ex23, ex13, compose, rename_*}` | on primal, dual and SL→PP outputs |
| Lemma 3.2 amplification + Claim 3.3 word count | `gadgets::amplify` | ≥ m^{2d}/(d+1) targets, all isomorphic to T_{a,b} |
| Lemma 3.4 SL-PP → PP (points i + εj) | `gadgets::sl_to_pp` | w targets via 2w−1 sources |
| Lemma 2.3 evaluation in eff(P) scalar products | `gadgets::evaluate` | source count equals eff(P) exactly |
| Lemma 4.1 Q_{a,b} | `gadgets::q_shape` | eff = (a+b−1)^{6d}, Q_{a,b} ≅ Q_{b,a}, Q_{1,1} = matmul |
| Lemma 4.2 the round | `gadgets::round` | exact and ε-perturbed inputs |
| Lemma 2.4 Bini | `gadgets::bini` | mod-ε chain → exact identity |
| Thm 4.3 chain | `gadgets::minimal_chain` | 2×2 integer matrices, equal to schoolbook |
| Section 4 parameter rules, ε bookkeeping, Bini's h | `src/cost.rs`, `bin/cost` | ε ranges per lemma equal those measured on materialized packings (`tests/profiles.rs`) |

### How a packing is represented and checked

A packing is three sparse matrices whose entries are Laurent polynomials in ε over F_p, with p = 2^61 − 1. They are the left encoder, the right encoder and the decoder, together with the target and source shapes. The verifier draws random inputs from F_p, encodes them, multiplies the source products, and decodes. It then checks that the result minus the true target output has only positive powers of ε. Because ε stays a symbol, "valid modulo ε" is an exact property here, not a numerical one. All of the paper's identities are rational, so they hold over F_p whenever p divides no denominator. A wrong identity survives one random point with probability below about 10^-15 at these sizes (Schwartz–Zippel). A negative control flips one decoder coefficient and is always caught.

Role exchange is done generically. A packing is a trilinear identity between targets and sources, so exchanging two roles only permutes the three maps (EL, ER, Dᵀ) and re-indexes coordinates with the paper's conventions, Prod^{12} = Prod_{J,I,K}[r,q,p] and Prod^{23} = Prod_{I,K⁻¹,J⁻¹}[q,p,r]. Isomorphisms (Lemma 2.10) are found by searching variable permutations, which takes at most 720 at d = 1.

## The finding: Lemma 3.2's output is exact

Lemma 3.2 tensors the primal packing with its dual d times and keeps the words with exactly d Full and d Mid factors. The paper says validity modulo ε is preserved. In fact every ε cancels.

Each column of the primal right encoder carries a single power of ε: ε^0 for the P_j, ε^1 for Q. Each row of its decoder does too: ε^0 for Full outputs, ε^−1 for the middle product. The dual swaps the two maps. A kept word has α primal Mid letters and β dual-of-Mid letters, and the balance condition is exactly α = β. So every entry of the amplified right encoder carries ε^(α−β) = ε^0, and every decoder entry ε^(β−α) = ε^0. The amplified packing, for example 17 targets from 9 sources at (a, b, c, t, d) = (2, 1, 3, 1, 1), contains no ε at all, and the verifier confirms it is exact (`verify` line "finding: 3.2 is exact").

The effect on cost is small. Without exactness, every round's ε-degree factor would grow from about 9w² to about (6d+3)²w². At the searched optimum that adds 17 bits per round to about 2,300, which makes the final identity about 0.8% larger in digits. The finding is about the structure of the construction: the formal parameter in Pass's chain is only needed for Lemma 3.4's interpolation.

## How the cost is computed

Every quantity follows a formula in Section 4 of the paper. `src/cost.rs` cites the equation next to each one.

- m = c+2, n = c; u = n^{2d}, w = ⌈m^{2d}/(d+1)⌉ (8); Λ = ⌊w/u⌋² (Lemma 4.2); δ_L = ⌈(3L−1)/(n−1)⌉, L_{j+1} = L_j + δ_{L_j}, L_0 = 1.
- After k rounds, the mod-ε identity has side v_k = w²Λ^k and uses r = E_k = 64^k w^9 (2L_k−1)^{6d} scalar multiplications (proof of Thm 4.3). The 64 counts padded sources, as in the paper.
- **ε-degree.** Every map's ε-exponent range is tracked: tensoring adds ranges, role exchange permutes them, and composition first substitutes ε → ε^N into the inner packing, with N = pL + pR + pD + 1 (proof of Lemma 2.9), then adds. One round multiplies the negative range by about 9w², and k rounds are chained, so d_ε ≈ (9w²)^k. Each rule is checked against ranges measured on materialized packings. At sizes too large to materialize, the closed forms are upper bounds on the same rules.
- **Bini** (Lemma 2.4). The exact identity has side V = v^h and uses R = (h·d_ε + 1)·r^h multiplications, with h the smallest power making log_V R < target.
- **"Pass"** uses the smallest c whose threshold is below the target and the smallest d satisfying (10). **"Searched"** applies the same lemmas over c ≤ 60 and d ≤ 2500. Both pick k from k_0 (the first feasible k) to 64·k_0, using 1,500 log-spaced samples plus local refinement, take the smallest feasible h at each k, and minimize R. This is a grid search, not a proof of optimality. Spot checks up to 80·k_0 and across d found the optimum in the interior.
- **Why not the smallest k.** At the very first feasible k, the mod-ε exponent sits within about 10^-15 (relative) of the target. The per-power margin in Bini's conversion is then below f64 resolution, so h cannot be computed reliably. Every reported identity is below its target by at least 64 bits of margin plus 10^-12 relative, which is far above rounding error.

**Comparing with Strassen.** The exact identity multiplies V×V matrices with R multiplications, so its recursion handles sides V, V², V³, … with R, R², R³, … products. Strassen at side V uses 7^⌈log2 V⌉ ≥ V^(log2 7) > R, which is what "exponent below 2.807" buys at the identity's own size. Below V, the identity can only be used by padding the input up to V. A careful implementation could then skip products that vanish on the padding, so no claim is made about smaller sides. Additions and the linear maps are not counted anywhere. The statement is about the identities in the table, not about every algorithm with exponent 9/4.

## What this verifies and what it does not

- **Verified:** finite instances of every lemma, built exactly as the proofs describe and checked exactly, plus the full chain end to end on 2×2 integer matrices. Also, that the cost model's ε bookkeeping matches the materialized constructions.
- **Not verified:** the proof in general. The checks cover all a, b ≤ 3, but they say nothing about all a, b, t, about the counting in Claim 3.3 beyond the instances run, or about the limit argument in Thm 4.3. Pass's paper reports a Lean formalization of its rational identities (about 3,600 lines, generated by Aristotle, which the paper says it has not independently verified). That covers the general statements, which this repo does not try to do. What this repo adds is concrete costs and a construction you can run.
- **The end-to-end run uses minimal retention.** The chain run keeps w = u = 1 (Observation 2.5 allows retaining fewer targets), so Λ = 1 and nothing is saved. It shows that primal → amplification → SL→PP feed the round correctly as separate pieces, and that the round, composition, Lemma 2.3 evaluation and Bini compose into a correct matrix multiplication. The source count matches the formula E_k = 64·3^6·x³. The smallest instance with real retention (c = 2, d = 1, w = 8, u = 4) needs about 10^14 scalar multiplications, so it is computed rather than run. No instance with an actual saving can ever be run.

## Design decisions

- **Exact over F_p, not floating point.** Bini-style algorithms cancel huge negative powers of ε. Any numerical ε would test approximation quality rather than the identity.
- **Materialized maps.** Each construction produces concrete sparse matrices, not a closure, so its ε-exponent ranges can be read off and compared with the cost model. The cost is memory. The 2×2 chain has 373,248 sources and takes a few minutes.
- **Paper's type labels for amplification.** When a + b − 1 = 1, Full and Mid have identical supports. The kept words are chosen by the labels in Claim 3.3 (primal Mid is the last target, the dual swaps types), not by looking at shapes. Choosing by shape keeps unbalanced words, which still carry ε.
- **Log space for the cost.** Exact big integers are used where they fit (u, w, ⌊w/u⌋ below 2^20000; L_k below 2^100). Beyond that, the code uses log2 values with upper-bound rounding.

## Limitations

- The search fixes the retained count at w from (8). It does not explore retaining fewer targets (smaller Λ, smaller ε-degree), the exact balanced count Σ_h C(d,h)²(c+1)^{2(d−h)} (larger than (8)), or the unpadded source count (2w−1)³ instead of 8w³. All three are small effects next to the ε-degree.
- Lemma 3.4 uses the paper's nodes i + εj. Other node choices with fewer negative powers might lower d_ε a lot. That would be a change to the construction, not a parameter, so it is left out.
- OpenAI's claimed lower bound is a separate, non-constructive argument and is not touched here.

## References

- R. Pass, *Matrix Multiplication from Polynomial Packings: A Constructive Proof of ω ≤ 9/4*, arXiv:2610.28627 (Oct 2026).
- OpenAI, *The matrix multiplication exponent is nine fourths* (preprint, Nov 24, 2026, github.com/openai/math), and the October upper-bound preprints.
- Doğan, Lysikov, Sakabe, Walter, *Tensor value and the 9/4 bound for matrix multiplication*, arXiv:2611.09716; *Semidirect sums of tensors and the matrix multiplication exponent*, arXiv:2611.04175.
- D. Bini, *Relations between exact and approximate bilinear algorithms. Applications*, Calcolo 17 (1980).
- P. Bürgisser, M. Clausen, A. Shokrollahi, *Algebraic Complexity Theory* (1997), Ch. 13–15.
- V. Strassen, *Gaussian elimination is not optimal*, Numer. Math. 13 (1969).

## License

MIT, Copyright (c) 2026 Ansor Juraev.
