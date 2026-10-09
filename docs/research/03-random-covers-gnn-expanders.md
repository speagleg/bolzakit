# 03 — Random covers of the Bolza surface as GNN expanders: prior art and feasibility

Status: read-only research note, 2026-10-09. Every arXiv ID below was found in a search result or fetched page during this session, except the ones marked *(from memory, not fetched)*. A quoted phrase means I read that text in the fetched page or abstract. A few numerical experiments were run in a scratch directory and their results are reported in §3. No bolzakit code was changed.

---

## TL;DR

- **The size-flexibility argument is already taken.** Exphormer (arXiv 2303.06147) builds its expander overlay by picking "d/2 permutations independently and uniformly among all possible permutations of n elements". That is exactly a random lift of a bouquet of loops, the Friedman/Bordenave–Collins model. It gives a 2d-regular, near-Ramanujan graph on **exactly N nodes for every N**. "Any size, provably near-optimal gap" is therefore not new relative to Exphormer.
- **The Bolza cover adds a relation, and the relation cannot improve the graph gap.** Take an n-sheeted cover of the one-vertex {8,8} cellulation of Bolza. Its 1-skeleton is the Schreier graph of 4 permutations satisfying the surface relation, which is 8-regular on exactly n vertices. By Magee–Puder–van Handel (arXiv 2504.08988, Thm 1.1), its nontrivial adjacency norm tends to ‖λ(x)‖ on the surface group. That equals 8·μ₂, where μ₂ ∈ [0.6623, 0.6629] (Nagnibeda; math/0109069). The Ramanujan value is 2√7/8 ≈ 0.6614. So asymptotically the gap is **about 0.1–0.2% worse than the free (Exphormer) model**, which is essentially equal and never better.
- **The genuinely new structure is the faces.** Each cover is a cellulation of a genus-(n+1) surface by n octagons. That makes it an **expander 2-complex of any size**, with a hyperbolic geometry and a deterministic congruence variant. That is the strongest honest claim. It matters for cellular, topological or sheaf message passing, not for plain MPNN rewiring.
- **Uniform sampling is solved for practical n.** I derived an exact rejection sampler and verified it is exactly uniform at n = 3. Its acceptance rate is ≈ 2/p(n), where p(n) is the partition count. That is fine up to n ≈ 40–60 in pure Python and needs a smarter sampler (character-weighted or MCMC) for n in the hundreds (§3).
- **Recommendation:** treat this as a **cheap negative-control-style experiment plus a niche positive angle** (2-complexes), not a headline result. Expect parity with Exphormer-style random regular overlays on OGB/LRGB.

---

## 1. GNN side: expander propagation and rewiring

### 1.1 Core papers

| Paper | arXiv | What it does | Expander and size handling |
|---|---|---|---|
| Deac, Lackenby, Veličković, *Expander Graph Propagation* (EGP) | 2210.02997 | Interleaves message passing over the input graph and over Cay(SL(2,ℤ_n); S_n). The graphs are "sparse and 4-regular", have "logarithmic diameter", and are precomputed independently of the input | Truncates the Cayley graph to the input size. Best paper at GLFrontiers 2022 |
| Wilson, Bechler-Speicher, Veličković, *Cayley Graph Propagation* (CGP) | 2410.03424 (LoG 2024) | Argues that truncated Cayley graphs "lose expansion properties", so it uses the **full** Cayley graph | Fills the size mismatch with **virtual nodes** (zeros, ones, N(0,1) or learnable). The most virtual nodes added range from 23 (MUTAG) to 1809 (REDDIT-BINARY). \|V\| = n³∏_{p\|n}(1−1/p²) |
| Mishra, Martinez, Lin, *Schreier Coset Graph Rewiring* (SCGR) | 2607.27479 | Schreier coset graphs of SL(2,ℤ_n) modulo the diagonal subgroup. These are 4-regular and have fewer vertices than the Cayley graph | Its "FindN" algorithm "searches for smallest n such that … \|V_Γ\| ≥ n_in". Overshoot is still possible. No mention of covers or lifts |
| Shirzad et al., *Exphormer* | 2303.06147 | Sparse graph transformer combining local edges, an expander overlay and virtual nodes | Expander built by "pick[ing] d/2 permutations independently and uniformly". The repo exposes `exp_algorithm ∈ {Random-d, Random-d2, Hamiltonian}`, with true degree 2·`exp_deg`. **Exactly N nodes, any N** |
| Christie, *Higher-Order Expander Graph Propagation* | OpenReview `0lZjoxTZFb` | Higher-order (hypergraph or bipartite) expander propagation | Could not read it, because OpenReview showed a bot-check page. It is relevant to the 2-complex angle and should be read before claiming novelty there |
| Topping et al., curvature rewiring (SDRF) | 2111.14522 | Proves "negatively curved edges are responsible for the over-squashing" and rewires using curvature | Input-dependent rewiring |
| Karhadkar et al., FoSR | 2210.11790 | Adds edges greedily to increase the spectral gap. Labels added edges as a separate relation | Input-dependent |
| Arnaiz-Rodríguez et al., DiffWire | 2206.07369 | Differentiable layers for commute times (CT-Layer) and spectral gap (GAP-Layer) | Learned |
| Qian et al., PR-MPNN | 2310.02156 | Learns which edges to add via k-subset sampling. Identifies "conditions under which they outperform purely randomized approaches" | Learned |
| Black et al., effective-resistance view | 2302.06835 | Over-squashing measured by total effective resistance | Theory |
| Surveys | 2308.15568, 2411.17429, 2605.00951 | Rewiring for over-squashing and over-smoothing | — |
| *Local Virtual Nodes* | 2508.20597 | Learnable local virtual nodes placed by centrality | — |
| *Understanding Virtual Nodes* (ICLR 2025) | 2405.13526 | Virtual nodes reduce commute time | Shows how strong the virtual-node baseline is |

### 1.2 Known pain points

1. **Size mismatch.** SL(2,ℤ_n) sizes are sparse: n = 2,3,4,5,6,7,… gives 6, 24, 48, 120, 144, 336, … nodes. EGP truncates, which hurts expansion according to CGP. CGP pads with up to ~1.8k virtual nodes, which costs compute and creates feature-initialisation choices. SCGR reduces the overshoot with coset graphs. **Random permutation models (Exphormer) have no mismatch at all.**
2. **Determinism and reproducibility.** Cayley graphs are deterministic, which is EGP's selling point. Random overlays need fixed seeds, or else vary per epoch, which can act as a regulariser. Congruence covers of Bolza would be deterministic, but they are size-restricted again (sizes are |PSL₂-type quotient| of ℤ[ζ₈]/𝔭).
3. **How much provable gap is needed?** Nothing in the literature I read ties downstream accuracy to a quantitative gap. EGP's theory is qualitative: there is a Cheeger lower bound, and negatively curved edges are necessary but not "sufficiently negatively curved" to cause over-squashing. Every 2d-regular random expander already sits at the Ramanujan bound up to o(1). **There is no accuracy headroom to gain from a "better" gap.**
4. **Benchmark fragility.** Tönshoff et al. (2309.00367) find that on LRGB "the performance gap completely vanishes after basic hyperparameter optimization". Any claimed improvement on LRGB needs tuned baselines.
5. **Node correspondence.** All these methods map input nodes to expander nodes arbitrarily (by index). The expander's geometry is therefore unrelated to the input graph's geometry. This weakens any "hyperbolic geometric structure" claim unless the embedding is learned or chosen.

---

## 2. Maths side: spectral gaps of random covers

### 2.1 Surface results

| Result | arXiv | Statement (as read) | Applies to Bolza? |
|---|---|---|---|
| Magee–Naud–Puder (GAFA 2022) | 2003.10911 | For a compact hyperbolic X and a uniform degree-n cover, "no new eigenvalues of X_n below 3/16 − ε" whp. They conjecture 1/4 | Yes (compact) |
| Hide–Magee (Annals 2023) | 2107.05292 | Finite-area **non-compact** X: no new eigenvalues in [0, 1/4 − ε) | **No.** It covers cusped surfaces only. Bolza is closed |
| **Magee–Puder–van Handel (2025)** | **2504.08988** | Thm 1.1: "For all x∈ℂ[Γ], for φₙ∈Hom(Γ,Sₙ) uniformly random, … ‖std∘φₙ(x)‖→‖λ(x)‖ in probability." Thm 1.2: for closed hyperbolic X, spec(Δ_{Xₙ}) ∩ [0, 1/4 − o(1)) = spec(Δ_X) ∩ [same] whp | **Yes. This is the result to cite for Bolza,** not Hide–Magee |
| Polynomial rate | 2505.08479 | "Spectral gap with polynomial rate for random covering surfaces" (title seen, not read) | Probably |
| Anantharaman–Monk | 2304.02678, 2502.12268 | Weil–Petersson random surfaces: λ₁ ≥ 1/4 − ε whp as g → ∞ | A different model (random moduli, not covers). Context only |
| Surveys | 2507.00346, 2510.12520 | The strong convergence phenomenon | — |

### 2.2 Graph results

- **Bordenave–Collins** (**1801.00876**, Annals 2019). Nontrivial eigenvalues of random n-lifts of a fixed base graph "approximately achieve the Alon-Boppana bound" whp. This generalises Friedman's theorem for Schreier graphs of independent random permutations.
- Friedman, Alon's conjecture: cs/0405020 *(from memory, not fetched)*. Bilu–Linial, lifts and near-optimal gap: math/0312022 *(from memory, not fetched)*.

### 2.3 Transferring the surface gap to graphs: you don't need to

The direct route is cleaner than any surface-to-graph comparison theorem.

- **{8,8} 1-skeleton (n vertices, 8-regular).** Its adjacency on the nontrivial part is std∘φₙ(x) with x = Σ(gᵢ + gᵢ⁻¹) ∈ ℂ[Γ]. MPvH Thm 1.1 applies verbatim. The norm tends to ‖λ(x)‖ = 8·μ₂(Γ₂), where μ₂ is the spectral radius of simple random walk on the genus-2 surface group with standard generators. Search results report the upper bound μ₂ ≤ 0.6629 (Nagnibeda) and the lower bound μ₂ ≥ 0.6623 (math/0109069; also Gouëzel 1310.4265). Ramanujan for 8-regular is 2√7/8 = 0.66144. So surface covers are **asymptotically near-Ramanujan, but strictly above the free bound.** The relation is a constraint, so it costs a little.
- **{8,3} tiling graph (16n vertices, cubic).** Its adjacency is a matrix-coefficient element of M₁₆(ℂ[Γ]). Strong convergence for scalar ℂ[Γ] is standard to upgrade to matrix coefficients; this is a standard fact, but I did not verify the statement in 2507.00346. The limit is the spectral radius of the infinite {8,3} tiling graph, which exceeds 2√2 (the 3-regular tree value).
- **Surface-to-graph comparison theorems (only needed for the Laplacian-on-surface story).**
  - Brooks (1986, *J. Diff. Geom.*, "spectral geometry of a tower of coverings") *(from memory, not fetched)*: λ₁ of covers is bounded below uniformly iff the Schreier graphs are uniform expanders. This is qualitative, with poor constants.
  - Mantuano's weighted discretisation, per the search summary, proves c·λ_k(M) ≤ λ_k(X) ≤ C·λ_k(M) under bounded geometry, with constants depending on curvature bounds and the discretisation parameter. This is also qualitative.
  - **Conclusion:** the "1/4 − ε" surface statement does **not** translate into a sharp graph statement via comparison. Use MPvH Thm 1.1 directly on the graph instead.

---

## 3. Generation

### 3.1 Sampling uniform φ ∈ Hom(Γ₂, S_n)

MPvH (2504.08988) "does not discuss algorithms for sampling" (fetched summary). I derived the following sampler and tested it.

**Exact rejection sampler (relator [a₁,b₁][a₂,b₂] = 1):**
1. Draw a₁ and b₁ uniformly from S_n and set c = [a₁,b₁].
2. Draw a **partition λ ⊢ n uniformly**, then a₂ uniformly in the conjugacy class λ. This makes P(a₂) ∝ |Z(a₂)|, which is the correct marginal weight.
3. A solution b₂ to [a₂,b₂] = c⁻¹ exists iff a₂⁻¹c⁻¹ is conjugate to a₂⁻¹. If not, reject and go back to step 1.
4. If a solution exists, the solution set is a coset b₀·Z(a₂⁻¹). Take b₀ from matching cycle structures and multiply by a uniform centralizer element (a wreath product, easy to sample).

Why it is exact: the target marginal on (a₁,b₁,a₂) is ∝ |Z(a₂)|·1[solvable], and the proposal density is ∝ |Z(a₂)|. So acceptance is just the indicator.

**Tests** (scratch scripts, Python):
- n = 3, 200k samples. All 486 homomorphisms were hit. 486 matches |Hom| = n!³·Σ_χ χ(1)⁻² = 216·2.25, and a brute-force enumeration. χ² = 523.6 on 485 df, consistent with uniform.
- Acceptance rates were 0.75 (n=3), 0.097 (n=8) and 0.0029 (n=20). These match ≈ 2/p(n), since |Hom| ≈ 2·n!³ (trivial and sign characters dominate).
- **Cost:** about p(n)/2 tries at O(n) each. That gives p(40) ≈ 3.7·10⁴ (seconds), p(60) ≈ 10⁶ (minutes, pure Python), and p(100) ≈ 1.9·10⁸ (too slow).
- **Fix for large n:** draw λ from its exact marginal using character values of S_n (Murnaghan–Nakayama), or use a Gibbs chain. In the Gibbs chain each move resamples one generator uniformly in its solution coset, plus mapping-class-group moves. The chain preserves the uniform measure. **Its irreducibility is unproven.**
- **bolzakit caveat:** bolzakit's side-pairing generators satisfy a different octagon relator (opposite-side pairings), not the standard commutator form. You need an explicit isomorphism to the standard presentation, or a sampler rederived for that relator.

**Small-n spectra.** I averaged 6 samples each, normalised by 8. "|max|" is the largest nontrivial eigenvalue in absolute value; λ₂ is the second-largest eigenvalue.

| Model | n = 30, \|max\| | n = 30, λ₂ | n = 40, \|max\| | n = 40, λ₂ |
|---|---|---|---|---|
| Surface cover | 0.654 | 0.621 | 0.638 | 0.622 |
| Free 4-permutation | 0.620 | 0.608 | 0.644 | 0.627 |

These are indistinguishable at this scale. Both are below the asymptotic ~0.662, as expected for small n.

### 3.2 Random lift of the tiling graph instead

- A **random lift of the {8,8} one-vertex graph without the relation is exactly the Exphormer/Friedman model.** It is trivial to sample, works for any n, and has the Bordenave–Collins guarantee (best possible).
- A random lift of the {8,3} graph (16 base vertices, 24 edges, one permutation per edge, ignoring faces) gives 16n vertices. It is near-Ramanujan relative to the base graph's universal cover, the 3-regular tree, which is better than the {8,3} tiling value.
- **What you lose:** the faces. A random lift is not a surface cover, octagons do not close up, and there is no hyperbolic geometry, no genus and no Bolza symmetry. **You lose nothing on the spectral gap and in fact gain slightly.**

### 3.3 Size flexibility

- {8,8} surface covers come in **every n ≥ 1**, as 8-regular multigraphs on exactly n vertices.
- {8,3} covers come in multiples of 16. Covers of the {3,8} map come in multiples of 6.
- Congruence covers are deterministic but restricted to group-order sizes, like SL(2,ℤ_n).
- So size flexibility is real versus EGP/CGP/SCGR, but **not versus Exphormer**.

---

## 4. Existing software

- **HyperCells** (GAP) and **HyperBloch** (Mathematica), by Lenggenhager, Maciejko and Bzdušek (PRL 131, 226401; 2305.04945). They build unit cells and supercells of {p,q} lattices as quotients by torsion-free normal subgroups, which are regular (normal) covers and therefore congruence-style. They are useful for cross-checking bolzakit congruence covers. They do not do random non-normal covers or GNNs.
- **networkx:** `random_regular_graph` exists. I did not verify a dedicated random-lift function; it is a ~20-line function either way.
- **PyG:** `transforms.VirtualNode` exists. I found no expander transform in the transforms docs.
- **DGL:** nothing found.
- **Reference implementations:** Exphormer at github.com/hamed1375/exphormer (GraphGPS-based, with the `exp_algorithm` config). EGP and CGP code are referenced from their papers; I did not inspect the repos.

---

## 5. A concrete, credible experiment

**Setup (copy CGP's protocol, 2410.03424):**
- GCN and GIN backbones.
- OGB: 5 layers, hidden 300, dropout 0.5, batch norm. LRGB: 5 layers, hidden 64. TUDataset: 4 layers, hidden 64.
- Datasets: ogbg-molhiv, Peptides-func, Peptides-struct, plus a few TU sets.
- **Skip ogbg-ppa.** CGP needed 8× Quadro RTX 8000 for it.

**Arms (all with the same interleaved-layer schedule):**
1. Input graph only, with tuned hyperparameters (per 2309.00367).
2. One global virtual node.
3. CGP: full Cay(SL(2,ℤ_n)) plus virtual-node padding.
4. **Free random 8-regular on exactly N nodes** (4 permutations, the Exphormer-style control).
5. **Bolza {8,8} random cover on exactly N nodes** (the treatment).
6. Bolza congruence cover (deterministic) plus padding.
7. Optional: an arm with the cover's octagonal faces as 2-cells in a cellular message-passing layer.

Use ≥ 10 seeds per arm and a fresh or fixed overlay per graph as an ablation.

**Success criteria (pre-register):**
- *Weak, likely achievable:* arm 5 ≈ arm 4 ≥ arm 3. This shows size-matched covers remove CGP's padding overhead at no accuracy cost. It is publishable only as a workshop note, since arm 4 already achieves it.
- *Real success:* arm 5 > arm 4 by more than seed std on ≥ 2 datasets, **or** arm 7 > arm 5 (the faces help). **Prior: low for the first, unknown for the second.**
- *Negative result also useful:* it confirms that the relation and geometry are irrelevant to MPNN expander overlays.

**Compute (estimates, not measured):**
- Overlay generation: N ≤ ~150 for molhiv and peptides, so the exact sampler costs seconds to minutes per distinct N in Python. Cache one overlay per N, or use a vectorised or C sampler.
- Training: ogbg-molhiv (41k graphs) with GIN-300 runs at roughly a few minutes per epoch on 12 CPU cores. That is ~100 epochs × 10 seeds × 6 arms, i.e. **days on the laptop**.
- Peptides runs are similar or larger.
- **Recommendation:** do the TU sets and a molhiv pilot (3 seeds) on CPU. Then move the full grid to a single cloud GPU (one A10 or 4090-class card), for roughly 1–2 GPU-days total.

---

## 6. Recommendation

- **Novelty.** No paper I found uses random covers of hyperbolic surfaces, or surface-group Schreier graphs, as GNN overlays. SCGR (2607.27479) is the closest: a group-theoretic Schreier graph, but of SL(2,ℤ_n). Using random permutations for size-flexible expanders is **not** novel (Exphormer).
- **Strongest honest claim:** "Size-matched expander overlays on exactly N nodes for every N. Each is the 1-skeleton of a hyperbolic surface cellulated by N octagons. Their nontrivial spectrum converges to the genus-2 surface-group bound 8·μ₂ ≈ 5.30 (vs 2√7 ≈ 5.29 Ramanujan) by Magee–Puder–van Handel. A deterministic congruence family is also available."
- **Do not claim** "better gap than random regular", and do not claim "1/4 − ε gap for graphs". The 1/4 is for the surface Laplacian, and for closed Bolza it comes from 2504.08988, not Hide–Magee.
- **Risks:**
  - (a) No accuracy gain over the Exphormer control. This is the most likely outcome.
  - (b) The sampler scales as p(n), so N in the hundreds needs the character-weighted or MCMC sampler, and MCMC irreducibility is unproven.
  - (c) bolzakit's relator differs from the standard commutator form.
  - (d) LRGB baseline tuning can erase small effects.
  - (e) The higher-order EGP work (Christie) may already cover the 2-complex angle; read it first.
- **Verdict:** worth a one- to two-week pilot because the infrastructure (bolzakit generators plus the sampler) is cheap. Frame it as a **controlled test of whether surface or hyperbolic structure matters beyond the spectral gap**, with arm 4 as the critical control. The 2-complex (faces) angle is the only part with a plausible novelty upside.
