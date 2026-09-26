# Engine notes: what the RWR learner does and why each test matters

This file is the conceptual companion to `d1/engine/rwr.py` and `tests/test_rwr.py`.
It explains the engine to anyone reading the code for the first time.

---

## 1. What RWR computes, in five sentences

1. The network is a sparse symmetric matrix A; A[i, j] = 1 if genes i and j are
   connected (unweighted, decision D-01).
2. W = D^{-1} A divides each row by the row's degree, so row i is the probability
   distribution of "step from i to one of its neighbours" — every row of W sums to 1.
3. A random walker starts on the seed genes (the training-fold known drivers, each
   column of P0 sums to 1). At every step it either restarts on a random seed
   (probability alpha) or moves to a random neighbour (probability 1 - alpha).
4. Iterating P <- alpha P0 + (1 - alpha) W^T P until the L1 change in every column
   is below 1e-8 gives the long-run probability P[i, j] of finding the walker at
   gene i for seed set j. Genes close to many seeds (in the network) score high.
5. Each column of P0 is one seed set (one fold). The harness propagates all 50
   seed vectors of an arm in one matrix product (batched), which is why the
   6-network x 4-label-set factorial finishes in seconds.

## 2. Why W^T and not W

The walker at gene j should move to gene i with probability proportional to A[i, j]
(which is A[j, i] since A is symmetric) divided by deg(j) (its out-degree under W).

- W[j, i] = A[j, i] / deg(j) — "from j, step to a neighbour including i".
- W^T[i, j] = W[j, i] — "to i, was at j", exactly what we want for the matrix product.

`transition_T` returns W^T (sparse CSR) so the iteration is `WT @ P` (no extra
transposition per step). The code also checks at construction that A is symmetric,
square, has no self-loops, and has no degree-0 nodes; the LCC restriction in
`finalize_edges` guarantees the last point.

## 3. Why alpha = 0.5 (and 0.3 / 0.7)

alpha is the **restart probability**. It controls how much the walker trusts the
seeds vs the network topology.

- alpha -> 0: P approaches the stationary distribution of W, which is the degree
  distribution (a "hub prior"). Useful as a baseline; we explicitly score degree
  separately so we can subtract its effect.
- alpha -> 1: P = P0 exactly, no network contribution.
- alpha = 0.5: balanced. The benchmark's primary analysis uses this value;
  0.3 and 0.7 are sensitivity arms to check that conclusions are not artefacts
  of one specific value.

## 4. Why unweighted (decision D-01)

We discard STRING's edge weights and use only the binary structure. Reason: weights
are evidence-type specific and hard to compare across networks (STRING combines
experiments, text mining, co-expression; IntAct combines affinity methods; Reactome
has curator scores). The benchmark is about **what topology contributes**, not how
good each network's weights are. Weight is kept in the network file (contract C2)
for future extensions but ignored in the analysis.

## 5. Why LCC only

A random walker cannot reach degree-0 nodes (no neighbours to step to), so RWR
returns P=0 for them. Including disconnected components would dilute the score
without changing the ordering within the LCC. `finalize_edges` strips everything
outside the LCC before this function ever runs, so `transition_T`'s check on
degree-0 is a backstop, not the primary restriction.

## 6. Why each test exists

I read every test, ran them, and asked "what bug would this catch?". The matrix
below is the answer.

| Test | Catches |
|---|---|
| `test_transition_is_row_stochastic` | Off-by-one in the row normalisation; column-normalising W instead of row-normalising (Bug B in `docs/a0_log.md`) |
| `test_mass_is_conserved` | Mass leak: W^T used without the right `deg` divisor; alpha/1-alpha swapped with each other so mass grows each iteration |
| `test_alpha_one_returns_seeds` | Off-by-one in the alpha term: `alpha * P0 + (1-alpha) * WT @ P` rewritten to `alpha * P0 + alpha * WT @ P` (Bug A) — at alpha=1 the buggy form adds `WT @ P0` on top of `P0`, this test catches that |
| `test_matches_closed_form[0.3 / 0.5 / 0.7]` | Closed-form is P = alpha (I - (1-alpha) W^T)^{-1} P0, derived by setting P_{k+1}=P_k. Any deviation in the alpha weighting, the matrix transpose direction, or the convergence loop will show up as a numerical mismatch. **Three alpha values, because at 0.5 a common bug (alpha vs 1-alpha confusion) is invisible.** |
| `test_star_centre_ranks_first` | Behavioural sanity: in a star graph with all leaves as seeds, the centre (node 0) is 1 step from all seeds and must outrank every leaf. This catches column-vs-row normalisation (Bug B) because column-normalisation makes leaves accumulate probability. |
| `test_signal_stays_in_seeded_community` | Behavioural sanity: in a barbell graph (two cliques connected by a single bridge), seeding one clique keeps the score higher there than in the other. Catches a bug where the walker leaks across the bridge too quickly. |
| `test_single_vector_input` | API: rwr must accept both (n,) and (n, k) inputs. The squeeze / unsqueeze logic on lines 74-77 must round-trip. |
| `test_rejects_degree_zero_nodes` | Input validation: degree-0 nodes would make W have a divide-by-zero. Fails fast rather than silently returning NaN. |
| `test_rejects_asymmetric_and_self_loops` | Input validation: directed or self-loopy A would invalidate the row-stochastic interpretation. |
| `test_empty_seed_set_rejected` | Input validation: empty P0 would make a column sum to 0 and the convergence check meaningless. |
| `test_warns_when_not_converged` | Robustness: alpha=0.01 with max_iter=3 on a 200-node path won't converge. We want a warning, not a silent wrong answer. |

## 7. Performance characteristics (FunMap-sized)

- 10,525 genes, 196,800 edges -> sparse transition matrix 4.8 MB
- Dense equivalent: 0.89 GB (177x larger)
- Per seed vector: ~7 ms on this PC (16-core, numpy 1.26.4)
- 50 seed vectors in one batched call: ~358 ms (50x cheaper than 50 separate calls)
- Full factorial (6 networks x 4 labels x 50 folds x 3 alphas): ~9 s

This is why we batch — the matrix-vector product is the bottleneck and scipy.sparse
runs it efficiently when all columns are stacked.

## 8. Why "engine-locked" matters for the benchmark

Any learner change between weeks would make the delta-AUROC comparisons
uninterpretable: a positive delta could be "the network helps" or "I tuned the
learner better this week". Locking the learner means the only variable we are
measuring is the network itself.

Person B's labels and Person A's networks change; the engine does not.